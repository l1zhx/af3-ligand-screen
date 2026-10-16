#!/usr/bin/env python3
"""Join AF3 to the existing Glide parent/Butina cluster table (run with RDKit)."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r)) or ["ligand_id"]
    with Path(path).open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def score(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


def parent_smiles(smiles):
    # Same primary standardization as cluster_glide_top1w.py; retain stereochemistry.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("cannot parse SMILES")
    mol = Chem.RemoveHs(mol, sanitize=True)
    mol = rdMolStandardize.FragmentParent(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def average_ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = np.empty(len(values), dtype=float)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        ranks[order[i:j]] = (i + 1 + j) / 2
        i = j
    return ranks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--af3-csv", required=True, help="AF3 ligand_metrics.csv; use all ligands for rank comparison")
    parser.add_argument("--geometry-csv", help="Optional ligand_geometry.csv for contact/pose annotations")
    parser.add_argument("--glide-csv", required=True, help="cluster_assignments.csv from cluster_glide_top1w.py")
    parser.add_argument("--glide-selection-csv", help="Optional diverse_top30/60 CSV to retain previous docking selections")
    parser.add_argument("--out", required=True)
    parser.add_argument("--top-fraction", type=float, default=0.10)
    parser.add_argument("--diverse-n", type=int, default=30)
    parser.add_argument("--plots", action="store_true")
    args = parser.parse_args()
    if not 0 < args.top_fraction <= 1 or args.diverse_n < 1:
        parser.error("top-fraction must be in (0,1], diverse-n positive")
    out = Path(args.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit("Output must be new or empty")
    out.mkdir(parents=True, exist_ok=True)
    af3 = read_csv(args.af3_csv)
    glide = read_csv(args.glide_csv)
    geometry = {r["ligand_id"]: r for r in read_csv(args.geometry_csv)} if args.geometry_csv else {}
    selected_parents = set()
    if args.glide_selection_csv:
        for r in read_csv(args.glide_selection_csv):
            selected_parents.add(parent_smiles(r["parent_smiles"]))
    if not glide or not {"parent_smiles", "docking_score", "cluster_id"} <= glide[0].keys():
        raise SystemExit("Need cluster_assignments.csv with parent_smiles, docking_score, cluster_id")
    if not af3 or not {"smiles", "interface_iptm", "ligand_id"} <= af3[0].keys():
        raise SystemExit("AF3 summary needs ligand manifest SMILES and interface_iptm")

    by_parent, issues = {}, []
    for r in glide:
        try:
            key = parent_smiles(r["parent_smiles"])
            ds = score(r["docking_score"])
            if ds is None:
                raise ValueError("invalid docking score")
            if key not in by_parent or ds < score(by_parent[key]["docking_score"]):
                by_parent[key] = r
        except ValueError as exc:
            issues.append(dict(source="glide", id=r.get("title", ""), reason=str(exc)))

    matched, unmatched = [], []
    for i, r in enumerate(af3, 1):
        try:
            key = parent_smiles(r["smiles"])
        except ValueError as exc:
            unmatched.append(dict(ligand_id=r["ligand_id"], reason=str(exc)))
            continue
        if key not in by_parent:
            unmatched.append(dict(ligand_id=r["ligand_id"], reason="absent_from_docking_parent_table"))
            continue
        d = by_parent[key]
        item = dict(r)
        item.update({k: v for k, v in geometry.get(r["ligand_id"], {}).items() if k not in item})
        item.update(matched_parent_smiles=key, identity_level="standardized_parent_not_exact_protonation_or_tautomer",
                    glide_title=d.get("title", ""), docking_score=d["docking_score"],
                    docking_ligand_efficiency=d.get("ligand_efficiency", ""), docking_hac=d.get("hac", ""),
                    docking_cluster_id=d["cluster_id"], docking_cluster_size=d.get("cluster_size", ""),
                    docking_pains=d.get("pains", ""), docking_review_flags=d.get("review_flags", ""),
                    docking_previously_selected=key in selected_parents,
                    docking_selection_priority=d.get("selection_priority", ""))
        matched.append(item)
        if i % 2000 == 0:
            print(f"Matched/standardized {i}/{len(af3)} AF3 states", flush=True)
    parent_states = defaultdict(list)
    for r in matched:
        parent_states[r["matched_parent_smiles"]].append(r)
    parents = []
    for key, states in parent_states.items():
        # One representative per standardized parent; keep the all-state table separately.
        states.sort(key=lambda r: (r.get("has_clash") != "False", score(r["interface_iptm"]) is None,
                                  -(score(r["interface_iptm"]) or 0), int(r["af3_rank"])))
        row = dict(states[0], af3_parent_state_count=len(states),
                   af3_parent_ligand_ids=";".join(r["ligand_id"] for r in states))
        parents.append(row)
    usable = [r for r in parents if score(r["interface_iptm"]) is not None and r.get("has_clash") == "False"]
    if usable:
        iptm = [score(r["interface_iptm"]) for r in usable]
        ds = [score(r["docking_score"]) for r in usable]
        rank_a = average_ranks([-v for v in iptm])
        rank_d = average_ranks(ds)
        n = len(usable)
        for r, ra, rd in zip(usable, rank_a, rank_d):
            r.update(af3_rank_in_matched_parents=float(ra), docking_rank_in_matched_parents=float(rd),
                     af3_rank_fraction=float(ra / n), docking_rank_fraction=float(rd / n),
                     consensus_top_fraction=ra / n <= args.top_fraction and rd / n <= args.top_fraction)
        # Both ranks: low is better. Positive rho means agreement of ranking.
        rho = float(np.corrcoef(rank_a, rank_d)[0, 1]) if n > 1 and np.std(rank_a) and np.std(rank_d) else None
    else:
        rho = None
    usable.sort(key=lambda r: (not r["consensus_top_fraction"],
                              max(r["af3_rank_fraction"], r["docking_rank_fraction"]),
                              r["af3_rank_fraction"], r["ligand_id"]))
    representatives, seen = [], set()
    for r in usable:
        # Require geometric inspection before choosing final review representatives.
        if r.get("geometry_triage_pass") != "True":
            continue
        pocket = r.get("pocket", "unassigned")
        cluster = r["docking_cluster_id"]
        k = (pocket, cluster)
        if k in seen:
            continue
        seen.add(k)
        representatives.append(dict(r, review_selection_rank=len(representatives) + 1))
        if len(representatives) >= args.diverse_n:
            break

    write_csv(out / "all_matched_states.csv", matched)
    write_csv(out / "matched_parent_comparison.csv", usable)
    write_csv(out / "consensus_parents.csv", [r for r in usable if r["consensus_top_fraction"]])
    write_csv(out / "docking_selected_af3_ids.csv", [dict(ligand_id=r["ligand_id"], glide_title=r["glide_title"])
              for r in matched if r["docking_previously_selected"]])
    write_csv(out / "parent_mapping.csv", parents)
    write_csv(out / "unmatched_af3.csv", unmatched)
    write_csv(out / "diverse_review.csv", representatives)
    write_csv(out / "mapping_errors.csv", issues)
    report = dict(af3_states=len(af3), glide_input_rows=len(glide), glide_unique_parents=len(by_parent),
                  matched_af3_states=len(matched), matched_parents=len(parents), usable_nonclashing_parents=len(usable),
                  spearman_rank_agreement=rho, top_fraction=args.top_fraction,
                  consensus_parents=sum(bool(r.get("consensus_top_fraction")) for r in usable),
                  diverse_review_count=len(representatives), unmatched_af3_states=len(unmatched),
                  limitations=["If docking input is top10000 only, comparison is conditional on that selected cohort.",
                               "Matching collapses protonation/tautomer states to parents; exact-state pose comparison still needs structures.",
                               "Geometry is only available for inspected candidates; unknown geometry is excluded from diverse_review.",
                               "Diverse review uses existing docking cluster IDs, not re-clustering the full AF3 library.",
                               "Rank consensus is an exploratory selection rule, not calibrated binding probability."])
    (out / "comparison.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    if args.plots and usable:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        axes[0].scatter([score(r["docking_score"]) for r in usable], [score(r["interface_iptm"]) for r in usable],
                        s=7, alpha=0.35, rasterized=True)
        axes[0].set(xlabel="Glide DockingScore (lower is better)", ylabel="AF3 A-B ipTM", title=f"Matched parents (n={len(usable)})")
        axes[1].scatter([r["docking_rank_fraction"] for r in usable], [r["af3_rank_fraction"] for r in usable],
                        s=7, alpha=0.35, rasterized=True)
        axes[1].set(xlabel="Glide rank fraction (lower is better)", ylabel="AF3 rank fraction (lower is better)",
                    title=f"Spearman rank agreement: {rho:.3f}" if rho is not None else "Insufficient variation")
        fig.savefig(out / "rank_comparison.png", dpi=300)
        fig.savefig(out / "rank_comparison.pdf")
        plt.close(fig)
        full_interface = [score(r["interface_iptm"]) for r in af3 if r.get("has_clash") == "False" and score(r["interface_iptm"]) is not None]
        pose_support = [score(r.get("pose_support_fraction")) for r in geometry.values() if score(r.get("pose_support_fraction")) is not None]
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
        axes[0].hist(full_interface, bins=50, color="#377eb8")
        axes[0].set(xlabel="AF3 A-B ipTM", ylabel="Ligand-state count",
                    title=f"Nonclashing AF3 library (n={len(full_interface)})")
        axes[1].hist(pose_support, bins=np.linspace(0, 1, 11), color="#4daf4a")
        axes[1].set(xlabel="Pose agreement fraction (analysis cutoff)", ylabel="Inspected ligand count",
                    title=f"Inspected subset only (n={len(pose_support)})")
        fig.savefig(out / "af3_confidence_and_pose_support.png", dpi=300)
        fig.savefig(out / "af3_confidence_and_pose_support.pdf")
        plt.close(fig)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
