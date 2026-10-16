"""Synthetic regression tests; these are not real AF3 results."""
import argparse
import csv
import json
import sys
import uuid
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import af3_postprocess as p


def cif_text(shift=(0, 0, 0), ligand_shift=0):
    columns = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
               "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num"]
    text = "data_synthetic\n#\nloop_\n" + "".join(f"_atom_site.{c}\n" for c in columns)
    atoms = [("A", 1, "CA", [0, 0, 0]), ("A", 2, "CA", [3, 0, 0]), ("A", 3, "CA", [0, 3, 0]),
             ("B", 1, "C1", [1, 1, 2.5]), ("B", 1, "C2", [1, 2, 2.5])]
    for i, (chain, residue, name, xyz) in enumerate(atoms, 1):
        xyz = np.array(xyz) + shift + ([0, 0, ligand_shift] if chain == "B" else [0, 0, 0])
        text += f"ATOM {i} C '{name}' . {'ALA' if chain == 'A' else 'LIG'} {chain} {residue} "
        text += " ".join(map(str, xyz)) + " 1\n"
    return text + "#\n"


class PostprocessTest(unittest.TestCase):
    def test_rigid_transform_and_displaced_ligand(self):
        atoms = [{"chain": c, "resid": r, "resname": "X", "name": n, "xyz": xyz}
                 for c, r, n, xyz in [("A", 1, "CA", [0, 0, 0]), ("A", 2, "CA", [3, 0, 0]),
                                      ("A", 3, "CA", [0, 3, 0]), ("B", 1, "C1", [1, 1, 2.5])]]
        ref = p.geometry(atoms, "A", "B")
        rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        rotated = [dict(a, xyz=np.array(a["xyz"]) @ rotation + [12, 4, 6]) for a in atoms]
        self.assertAlmostEqual(p.aligned_pose_rmsd(ref, p.geometry(rotated, "A", "B")), 0, places=8)
        displaced = [dict(a, xyz=np.array(a["xyz"]) + ([0, 0, 4] if a["chain"] == "B" else [0, 0, 0])) for a in atoms]
        self.assertAlmostEqual(p.aligned_pose_rmsd(ref, p.geometry(displaced, "A", "B")), 4, places=8)

    def test_audit_duplicates_partial_mapping_geometry(self):
        root = Path(__file__).parent / ("af3_test_" + uuid.uuid4().hex[:8])
        root.mkdir()
        self.addCleanup(__import__("shutil").rmtree, root)
        results = root / "results"
        for ligand, attempt, count, iptm in [(1, 100, 5, .85), (1, 200, 5, .99), (2, 300, 2, .95)]:
            for sample in range(count):
                folder = (results / "000" / f"lig_{ligand:06d}" / "seed_1" / f"attempt_{attempt}_0" /
                          f"bcat_lig_{ligand:06d}_seed_1" / f"seed-1_sample-{sample}")
                folder.mkdir(parents=True)
                prefix = "synthetic"
                summary = dict(ranking_score=iptm, has_clash=False, chain_ids=["B", "A"],
                               chain_pair_iptm=[[.1, iptm], [iptm, .9]], chain_pair_pae_min=[[0, 3], [6, 0]])
                (folder / f"{prefix}_summary_confidences.json").write_text(json.dumps(summary))
                (folder / f"{prefix}_model.cif").write_text(cif_text(ligand_shift=4 if sample == 4 else 0))
                confidence = dict(atom_chain_ids=["A"] * 3 + ["B"] * 2,
                                  atom_plddts=[90, 90, 90, 60, 80], token_chain_ids=["A"] * 3 + ["B"] * 2,
                                  pae=np.full((5, 5), 4.0).tolist())
                (folder / f"{prefix}_confidences.json").write_text(json.dumps(confidence))
        # The root-level best file must not count as a sixth sample.
        (results / "000/lig_000001/seed_1/best_summary_confidences.json").write_text('{"ranking_score":1.0}')
        manifest = root / "ligands.csv"
        manifest.write_text("id,inchikey,smiles,source_name\n1,K1,CC,F1\n2,K2,CN,F2\n3,K3,CO,F3\n")
        p.summarize(argparse.Namespace(results=results, jobs=None, ligands=manifest, out=root / "summary",
                    chain_order="A,B", samples_per_seed=5, protein_chain="A", ligand_chain="B"))
        rows = p.read_csv(root / "summary/ligand_metrics.csv")
        self.assertEqual(len(rows), 1)
        self.assertEqual(float(rows[0]["interface_iptm"]), .85)
        self.assertEqual(float(rows[0]["pae_min_A_to_B"]), 6)
        self.assertEqual(int(rows[0]["n_samples"]), 5)
        self.assertEqual(len(p.read_csv(root / "summary/duplicate_attempts.csv")), 1)
        status = {r["ligand_id"]: r["status"] for r in p.read_csv(root / "summary/coverage.csv")}
        self.assertEqual(status["lig_000002"], "partial_or_invalid")
        self.assertEqual(status["lig_000003"], "missing")
        p.inspect(argparse.Namespace(analysis=root / "summary", out=root / "geometry", top=100, export_top=100,
                  pockets=None, include_ids=None, residue_offset=130, contact_cutoff=4.0, pose_cutoff=2.0))
        detail = p.read_csv(root / "geometry/ligand_geometry.csv")[0]
        self.assertEqual(float(detail["ligand_plddt_mean"]), 70)
        self.assertEqual(float(detail["pose_support_fraction"]), .8)
        self.assertEqual(detail["contact_residues_reference"], "131;132;133")
        self.assertEqual(detail["geometry_triage_pass"], "True")
        self.assertEqual(len(list((root / "geometry/top_structures").glob("*/*_model.cif"))), 1)

        # Older AF3 releases write bare sample filenames. Both plural and
        # singular summary names must find model.cif and confidences.json.
        for path in list(results.rglob("synthetic_*")):
            path.rename(path.with_name(path.name.removeprefix("synthetic_")))
        for suffix in ("plural", "singular"):
            if suffix == "singular":
                for path in list(results.rglob("summary_confidences.json")):
                    path.rename(path.with_name("summary_confidence.json"))
            p.summarize(argparse.Namespace(results=results, jobs=None, ligands=manifest,
                out=root / ("legacy_" + suffix), chain_order="A,B", samples_per_seed=5,
                protein_chain="A", ligand_chain="B"))
            rows = p.read_csv(root / ("legacy_" + suffix) / "ligand_metrics.csv")
            self.assertEqual(len(rows), 1)
            self.assertEqual(float(rows[0]["interface_iptm"]), .85)

        empty = root / "unrecognized_results"
        empty.mkdir()
        (empty / "best_summary_confidences.json").write_text('{"ranking_score":0.9}')
        with self.assertRaisesRegex(SystemExit, "No recognized sample summaries"):
            p.summarize(argparse.Namespace(results=empty, jobs=None, ligands=manifest,
                out=root / "unrecognized_audit", chain_order="A,B", samples_per_seed=5,
                protein_chain="A", ligand_chain="B"))
        unresolved = p.read_csv(root / "unrecognized_audit/coverage.csv")
        self.assertTrue(all(r["status"] == "not_resolved_by_scan" for r in unresolved))


if __name__ == "__main__":
    unittest.main()
