#!/usr/bin/env python3
"""Audit AF3 samples, rank interfaces, and inspect a shortlist. NumPy only."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shlex
import shutil
import stat
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SAMPLE = re.compile(r"seed-(\d+)_sample-(\d+)$")
LIGAND = re.compile(r"lig_\d{6}$")
JOB = re.compile(r"bcat_(lig_\d{6})_(.+)\.json$")
SUMMARY_ENDINGS = ("summary_confidences.json", "summary_confidence.json")


def summary_prefix(name):
    """AF3 releases use bare or job-prefixed sample output filenames."""
    for ending in SUMMARY_ENDINGS:
        if name == ending:
            return ""
        if name.endswith("_" + ending):
            return name[:-len(ending)]  # Retain the separator underscore.
    return None


def companion_files(summary):
    prefix = summary_prefix(summary.name)
    if prefix is None:
        raise ValueError(f"Unsupported summary filename: {summary.name}")
    return (summary.with_name(prefix + "model.cif"),
            summary.with_name(prefix + "confidences.json"))


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def new_dir(path):
    path = Path(path).resolve()
    if path.exists() and any(path.iterdir()):
        raise SystemExit(f"Output must be a new or empty directory: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def truth(value):
    if value in (True, 1, "1", "True", "true"):
        return True
    if value in (False, 0, "0", "False", "false"):
        return False
    return None


def pair(data, key, order, a, b):
    try:
        return number(data[key][order.index(a)][order.index(b)])
    except (ValueError, KeyError, IndexError, TypeError):
        return None


def priority(row):
    # Only known non-clashing structures lead the queue; missing metrics remain missing.
    clash = truth(row.get("has_clash"))
    interface = number(row.get("interface_iptm"))
    ranking = number(row.get("ranking_score"))
    return (0 if clash is False else 1 if clash is None else 2,
            0 if interface is not None else 1,
            -(interface if interface is not None else -math.inf),
            -(ranking if ranking is not None else -math.inf), row.get("summary_path", ""))


def attempt_order(path):
    match = re.match(r"attempt_(\d+)(?:_\d+)?$", Path(path).name)
    return (int(match.group(1)) if match else math.inf, path)


def summarize(args):
    root = Path(args.results).resolve()
    if not root.is_dir():
        raise SystemExit(f"Results directory not found: {root}")
    out = new_dir(args.out)
    groups = defaultdict(list)
    errors = []
    found = set()
    n_files = 0
    started = time.monotonic()
    last_report = started
    n_directories = 0
    discovered_summaries = 0
    skipped_directory_examples = []
    print(f"Starting summary scan: {root}; sample summaries are read serially", flush=True)

    def paths():
        nonlocal last_report, n_directories, discovered_summaries
        def scan_error(exc):
            errors.append(dict(path=str(exc.filename or root), error=f"directory scan failed: {exc}"))
        for current, dirs, files in os.walk(root, onerror=scan_error):
            n_directories += 1
            now = time.monotonic()
            if now - last_report >= 30:
                elapsed = now - started
                print(f"Scan elapsed={elapsed / 60:.1f} min; directories={n_directories}; "
                      f"sample summaries={n_files}; ligands seen={len(found)}; "
                      f"rate={n_files / elapsed:.1f} samples/s", flush=True)
                last_report = now
            folder = Path(current)
            for name in files:
                if summary_prefix(name) is not None:
                    discovered_summaries += 1
                    if SAMPLE.fullmatch(folder.name):
                        yield folder / name
                    elif len(skipped_directory_examples) < 10:
                        skipped_directory_examples.append(str(folder / name))

    for summary in paths():
        match = SAMPLE.fullmatch(summary.parent.name)
        if not match:
            continue  # Root-level best outputs are copies, never additional samples.
        ligand_dir = next((p for p in summary.parents if LIGAND.fullmatch(p.name)), None)
        if ligand_dir is None:
            continue
        ligand = ligand_dir.name
        found.add(ligand)
        relative = summary.relative_to(ligand_dir).parts
        round_name = relative[0]
        attempt = next((p for p in summary.parents if p.name.startswith("attempt_")), summary.parent.parent)
        n_files += 1
        try:
            data = json.loads(summary.read_text())
            score = number(data.get("ranking_score"))
            if score is None:
                raise ValueError("missing or nonfinite ranking_score")
            cif, confidence = companion_files(summary)
            cif_stat, confidence_stat = cif.stat(), confidence.stat()
            if not stat.S_ISREG(cif_stat.st_mode) or not stat.S_ISREG(confidence_stat.st_mode) or not cif_stat.st_size or not confidence_stat.st_size:
                raise ValueError("matching model/confidence missing or empty")
            order = data.get("chain_ids") or args.chain_order.split(",")
            if args.protein_chain not in order or args.ligand_chain not in order:
                raise ValueError(f"chain IDs absent from order: {order}")
            row = dict(ligand_id=ligand, round=round_name, seed=int(match[1]), sample=int(match[2]),
                       ranking_score=score, ptm=number(data.get("ptm")), iptm=number(data.get("iptm")),
                       interface_iptm=pair(data, "chain_pair_iptm", order, args.protein_chain, args.ligand_chain),
                       pae_min_A_to_B=pair(data, "chain_pair_pae_min", order, args.protein_chain, args.ligand_chain),
                       pae_min_B_to_A=pair(data, "chain_pair_pae_min", order, args.ligand_chain, args.protein_chain),
                       has_clash=truth(data.get("has_clash")), fraction_disordered=number(data.get("fraction_disordered")),
                       chain_order_source="summary" if data.get("chain_ids") else "explicit_cli_order",
                       attempt_path=str(attempt), summary_path=str(summary), cif_path=str(cif),
                       confidence_path=str(confidence))
            groups[(ligand, round_name, int(match[1]), str(attempt))].append(row)
        except (ValueError, OSError, TypeError) as exc:
            errors.append(dict(path=str(summary), error=str(exc)))
        if n_files % 10000 == 0:
            print(f"Read {n_files} sample summaries", flush=True)

    seed_groups = defaultdict(list)
    print(f"Scan finished: {n_files} sample summaries; validating groups and duplicate attempts", flush=True)
    for (ligand, rnd, seed, attempt), rows in groups.items():
        samples = {r["sample"] for r in rows}
        valid = len(rows) == len(samples) == args.samples_per_seed and samples == set(range(args.samples_per_seed))
        if valid:
            seed_groups[(ligand, rnd, seed)].append((attempt, rows))
        else:
            errors.append(dict(path=attempt, error=f"partial/duplicate seed {seed}: sample indices {sorted(samples)}"))
    selected, duplicates = [], []
    by_ligand = defaultdict(list)
    for (ligand, rnd, seed), attempts in sorted(seed_groups.items()):
        attempts.sort(key=lambda x: attempt_order(x[0]))
        chosen, rows = attempts[0]
        selected.extend(rows)
        by_ligand[ligand].extend(rows)
        for discarded, other in attempts[1:]:
            duplicates.append(dict(ligand_id=ligand, round=rnd, seed=seed, kept_attempt=chosen,
                                   excluded_attempt=discarded, excluded_samples=len(other)))

    manifest = {}
    if args.ligands:
        for r in read_csv(args.ligands):
            ligand = r.get("ligand_id") or f"lig_{int(r['id']):06d}"
            if ligand in manifest:
                raise SystemExit(f"Duplicate ligand ID in manifest: {ligand}")
            manifest[ligand] = r
    expected = set(manifest)
    if args.jobs:
        for path in Path(args.jobs).glob("*/*.json"):
            match = JOB.fullmatch(path.name)
            if match:
                expected.add(match[1])
    expected = expected or found
    coverage = [dict(ligand_id=ligand, status="not_resolved_by_scan" if n_files == 0 else
                     "complete_sample_group" if ligand in by_ligand else
                     "partial_or_invalid" if ligand in found else "missing", selected_samples=len(by_ligand.get(ligand, [])))
                for ligand in sorted(expected)]
    candidates = []
    for ligand, rows in by_ligand.items():
        best = min(rows, key=priority)
        rank_best = max(rows, key=lambda r: r["ranking_score"])
        item = dict(best)
        for key in ("inchikey", "smiles", "source_name", "source_file", "source_record"):
            item[key] = manifest.get(ligand, {}).get(key, "")
        metrics = [r["interface_iptm"] for r in rows if r["interface_iptm"] is not None]
        item.update(n_samples=len(rows), n_seeds=len({r["seed"] for r in rows}),
                    interface_iptm_mean=float(np.mean(metrics)) if metrics else None,
                    interface_iptm_sd=float(np.std(metrics)) if metrics else None,
                    max_ranking_score=rank_best["ranking_score"], ranking_best_cif=rank_best["cif_path"],
                    known_nonclashing_samples=sum(r["has_clash"] is False for r in rows))
        candidates.append(item)
    candidates.sort(key=priority)
    for rank, row in enumerate(candidates, 1):
        row["af3_rank"] = rank
    write_csv(out / "sample_metrics.csv", selected, list(selected[0]) if selected else ["ligand_id"])
    write_csv(out / "ligand_metrics.csv", candidates, list(candidates[0]) if candidates else ["ligand_id"])
    write_csv(out / "coverage.csv", coverage, ["ligand_id", "status", "selected_samples"])
    write_csv(out / "duplicate_attempts.csv", duplicates,
              ["ligand_id", "round", "seed", "kept_attempt", "excluded_attempt", "excluded_samples"])
    write_csv(out / "errors.csv", errors, ["path", "error"])
    audit = dict(results=str(root), expected_ligands=len(expected), ranked_ligands=len(candidates),
                 coverage=dict(Counter(r["status"] for r in coverage)), selected_samples=len(selected),
                 duplicate_seed_groups_excluded=len(duplicates), error_records=len(errors),
                 ligands_missing_interface_metric=sum(r["interface_iptm"] is None for r in candidates),
                 unexpected_result_ligands=sorted(found - expected),
                 discovered_summary_files=discovered_summaries, recognized_sample_summary_files=n_files,
                 visited_directories=n_directories, skipped_non_sample_summary_examples=skipped_directory_examples,
                 samples_per_seed=args.samples_per_seed, protein_chain=args.protein_chain, ligand_chain=args.ligand_chain,
                 ranking_method="known nonclashing, then A-B ipTM descending, then ranking_score",
                 duplicate_policy="earliest complete attempt per ligand/round/seed; no file deletion")
    (out / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2), flush=True)
    if n_files == 0:
        raise SystemExit("No recognized sample summaries. Do not interpret this as failed predictions. "
                         "Check results path, directory names, and audit skipped_non_sample_summary_examples; "
                         "run: python af3_postprocess.py probe --results <results>")


def probe(args):
    root = Path(args.results).resolve()
    if not root.is_dir():
        raise SystemExit(f"Results directory not found: {root}")
    # Look at a few ligands, never perform a full-library scan for diagnosis.
    ligands = []
    if LIGAND.fullmatch(root.name):
        ligands = [root]
    else:
        for shard in sorted(root.iterdir()):
            if not shard.is_dir():
                continue
            if LIGAND.fullmatch(shard.name):
                ligands.append(shard)
            else:
                for ligand in sorted(shard.iterdir()):
                    if ligand.is_dir() and LIGAND.fullmatch(ligand.name):
                        ligands.append(ligand)
                        if len(ligands) >= args.max_ligands:
                            break
            if len(ligands) >= args.max_ligands:
                break
    print(f"Probe results root: {root}; examining {len(ligands)} ligand directories", flush=True)
    if not ligands:
        print("No lig_NNNNNN directories found at the expected depth. Root entries:")
        for item in list(root.iterdir())[:20]:
            print(item.name)
    for ligand in ligands:
        print(f"LIGAND {ligand}", flush=True)
        listed = 0
        for current, dirs, files in os.walk(ligand):
            for name in sorted(files):
                if "confiden" not in name and not name.endswith(".cif"):
                    continue
                path = Path(current) / name
                print(f"  {path.relative_to(root)}", flush=True)
                if summary_prefix(name) is not None:
                    cif, confidence = companion_files(path)
                    print(f"    sample_directory_match={bool(SAMPLE.fullmatch(path.parent.name))}; "
                          f"model_exists={cif.is_file()}; confidence_exists={confidence.is_file()}")
                listed += 1
                if listed >= 18:
                    break
            if listed >= 18:
                break


def read_atoms(path):
    """Read the AF3 _atom_site loop; quoted names and wrapped rows supported.

    This intentionally reads AF3 atom loops, not a general CIF document library.
    Ambiguous alt locations/multiple models raise instead of being counted twice.
    """
    rows, headers, pending = [], [], []
    collecting = False
    with Path(path).open() as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("_atom_site."):
                headers.append(line.split()[0].split(".", 1)[1])
                collecting = True
                continue
            if not collecting or not headers:
                continue
            if not line or line.startswith("#"):
                if rows:
                    break
                continue
            if line.startswith(("loop_", "data_", "_", ";")):
                if rows:
                    break
                raise ValueError("unsupported/empty atom_site loop")
            pending.extend(shlex.split(line, comments=True, posix=True))
            while len(pending) >= len(headers):
                values, pending = pending[:len(headers)], pending[len(headers):]
                r = dict(zip(headers, values))
                if r.get("type_symbol", "").upper() in ("H", "D"):
                    continue
                if r.get("pdbx_PDB_model_num", "1") != "1" or r.get("label_alt_id", ".") not in (".", "?", "A"):
                    raise ValueError("multiple models/alternate positions require explicit handling")
                chain = r["label_asym_id"]
                resid = r.get("label_seq_id", ".")
                if resid in (".", "?"):
                    resid = r.get("auth_seq_id", "0")
                xyz = [float(r[c]) for c in ("Cartn_x", "Cartn_y", "Cartn_z")]
                if not all(math.isfinite(x) for x in xyz):
                    raise ValueError("nonfinite atom coordinate")
                rows.append(dict(chain=chain, resid=int(resid), resname=r["label_comp_id"], name=r["label_atom_id"], xyz=xyz))
    if pending or not rows:
        raise ValueError("truncated or empty atom_site loop")
    return rows


def geometry(atoms, protein, ligand, contact_cutoff=4.0):
    a = [r for r in atoms if r["chain"] == protein]
    b = [r for r in atoms if r["chain"] == ligand]
    if not a or not b:
        raise ValueError("protein/ligand atoms missing")
    ca = {r["resid"]: np.array(r["xyz"]) for r in a if r["name"] == "CA"}
    names = [r["name"] for r in b]
    if len(names) != len(set(names)):
        raise ValueError("ligand atom names not unique")
    distances = np.linalg.norm(np.array([r["xyz"] for r in a])[:, None, :] - np.array([r["xyz"] for r in b])[None, :, :], axis=2)
    residue_min = {}
    for atom, d in zip(a, distances.min(axis=1)):
        residue_min[atom["resid"]] = min(residue_min.get(atom["resid"], math.inf), float(d))
    contacts = {resid for resid, d in residue_min.items() if d <= contact_cutoff}
    return dict(min_distance=float(distances.min()), short_contacts=int(np.sum(distances < 1.5)),
                contact_pairs=int(np.sum(distances <= contact_cutoff)), ligand_hac=len(b), contacts=contacts,
                residue_min=residue_min, ca=ca, ligand={r["name"]: np.array(r["xyz"]) for r in b})


def aligned_pose_rmsd(reference, moving, cutoff=10.0):
    # Protein-only fit. Never fit the ligand to itself, which hides different sites.
    residues = [r for r, d in reference["residue_min"].items() if d <= cutoff and r in reference["ca"] and r in moving["ca"]]
    if len(residues) < 3 or set(reference["ligand"]) != set(moving["ligand"]):
        return None
    x, y = np.array([moving["ca"][r] for r in residues]), np.array([reference["ca"][r] for r in residues])
    x0, y0 = x.mean(axis=0), y.mean(axis=0)
    if np.linalg.matrix_rank(x - x0) < 2:
        return None
    u, _, vt = np.linalg.svd((x - x0).T @ (y - y0))
    fix = np.diag([1.0, 1.0, np.linalg.det(u @ vt)])
    rotation = u @ fix @ vt
    names = sorted(reference["ligand"])
    ligand_x = (np.array([moving["ligand"][n] for n in names]) - x0) @ rotation + y0
    ligand_y = np.array([reference["ligand"][n] for n in names])
    return float(np.sqrt(np.mean(np.sum((ligand_x - ligand_y) ** 2, axis=1))))


def full_confidence(path, protein, ligand):
    data = json.loads(Path(path).read_text())
    labels = data.get("atom_chain_ids", [])
    plddt = data.get("atom_plddts", [])
    if len(labels) != len(plddt):
        raise ValueError("atom_chain_ids / atom_plddts lengths differ")
    vals = [number(v) for c, v in zip(labels, plddt) if c == ligand]
    vals = [v for v in vals if v is not None]
    out = dict(ligand_plddt_mean=float(np.mean(vals)) if vals else None)
    tokens = data.get("token_chain_ids", [])
    ia = [i for i, c in enumerate(tokens) if c == protein]
    ib = [i for i, c in enumerate(tokens) if c == ligand]
    pae = data.get("pae")
    if ia and ib and pae is not None:
        pae = np.asarray(pae, dtype=float)
        if pae.shape != (len(tokens), len(tokens)):
            raise ValueError("PAE shape does not match token_chain_ids")
        if not np.isfinite(pae).all():
            raise ValueError("nonfinite PAE values")
        for name, x, y in (("A_to_B", ia, ib), ("B_to_A", ib, ia)):
            block = pae[np.ix_(x, y)]
            out[f"pae_mean_{name}"] = float(np.mean(block))
            out[f"pae_min_{name}"] = float(np.min(block))
    return out


def inspect(args):
    analysis = Path(args.analysis)
    all_candidates = read_csv(analysis / "ligand_metrics.csv")
    extra = set()
    if args.include_ids:
        extra = {r["ligand_id"] for r in read_csv(args.include_ids)}
    candidates = [r for i, r in enumerate(all_candidates) if i < args.top or r["ligand_id"] in extra]
    if not candidates:
        raise SystemExit("No complete ligand sample groups in ligand_metrics.csv")
    audit = json.loads((analysis / "audit.json").read_text())
    protein, ligand = audit["protein_chain"], audit["ligand_chain"]
    by_ligand = defaultdict(list)
    chosen_ids = {r["ligand_id"] for r in candidates}
    for r in read_csv(analysis / "sample_metrics.csv"):
        if r["ligand_id"] in chosen_ids:
            by_ligand[r["ligand_id"]].append(r)
    pockets = json.loads(Path(args.pockets).read_text()) if args.pockets else {}
    out = new_dir(args.out)
    detail, sample_detail, contacts_out, errors = [], [], [], []
    hotspot = Counter()
    started = last_report = time.monotonic()
    print(f"Starting structural inspection: {len(candidates)} candidates (not the full library)", flush=True)
    for i, candidate in enumerate(candidates, 1):
        lid = candidate["ligand_id"]
        models = []
        for r in by_ligand[lid]:
            now = time.monotonic()
            if now - last_report >= 30:
                print(f"Inspection elapsed={(now - started) / 60:.1f} min; "
                      f"ligand {i}/{len(candidates)}: {lid}; completed={len(detail)}", flush=True)
                last_report = now
            try:
                g = geometry(read_atoms(r["cif_path"]), protein, ligand, args.contact_cutoff)
                metrics = full_confidence(r["confidence_path"], protein, ligand)
                metrics.update(min_heavy_atom_distance=g["min_distance"], short_contact_pairs_lt_1_5A=g["short_contacts"],
                               contact_pairs=g["contact_pairs"], ligand_hac=g["ligand_hac"],
                               contact_residues_af3=";".join(map(str, sorted(g["contacts"]))),
                               contact_residues_reference=";".join(str(x + args.residue_offset) for x in sorted(g["contacts"])))
                models.append((r, g, metrics))
            except (ValueError, KeyError, OSError, TypeError) as exc:
                errors.append(dict(ligand_id=lid, path=r["cif_path"], error=str(exc)))
        reference = next((m for m in models if m[0]["summary_path"] == candidate["summary_path"]), None)
        if reference is None:
            continue
        r0, g0, m0 = reference
        comparable, support = 0, 0
        for r, g, metrics in models:
            rmsd = aligned_pose_rmsd(g0, g)
            metrics["pose_rmsd_to_best_A"] = rmsd
            metrics["pose_rmsd_method"] = "same-name ligand atoms; reference-pocket CA alignment; no symmetry correction"
            if rmsd is not None:
                comparable += 1
                support += rmsd <= args.pose_cutoff
            sample_detail.append(dict(r, **metrics))
        contact_counts = Counter(x for _, g, _ in models for x in g["contacts"])
        for res, n in sorted(contact_counts.items()):
            contacts_out.append(dict(ligand_id=lid, af3_residue=res, reference_residue=res + args.residue_offset,
                                     contacting_samples=n, analyzed_samples=len(models), frequency=n / len(models)))
        hotspot.update(g0["contacts"])
        overlaps = {p: len(g0["contacts"] & {int(x) - args.residue_offset for x in residues}) for p, residues in pockets.items()}
        best_pocket = max(overlaps, key=overlaps.get) if overlaps and max(overlaps.values()) > 0 else "unassigned"
        flags = []
        if truth(candidate["has_clash"]) is not False:
            flags.append("global_clash_or_unknown")
        if g0["short_contacts"]:
            flags.append("very_short_atom_contacts_review")
        if not g0["contacts"]:
            flags.append("no_contacts_within_cutoff")
        if number(candidate.get("interface_iptm")) is None:
            flags.append("interface_metric_missing")
        if m0["ligand_plddt_mean"] is None:
            flags.append("ligand_confidence_missing")
        if len(models) != len(by_ligand[lid]):
            flags.append("some_samples_unreadable")
        if comparable < len(models):
            flags.append("some_pose_comparisons_unavailable")
        row = dict(candidate, **m0, pose_comparable_samples=comparable,
                   pose_support_fraction=support / comparable if comparable else None,
                   n_geometry_samples=len(models), pocket=best_pocket,
                   pocket_contact_counts=json.dumps(overlaps, separators=(",", ":")), review_flags=";".join(flags),
                   # Geometric triage only; chemical/functional validation remains external.
                   geometry_triage_pass=truth(candidate["has_clash"]) is False and not g0["short_contacts"] and bool(g0["contacts"]))
        detail.append(row)
        if i % 50 == 0:
            print(f"Inspected {i}/{len(candidates)} ligands", flush=True)
    write_csv(out / "ligand_geometry.csv", detail, None if detail else ["ligand_id"])
    write_csv(out / "sample_geometry.csv", sample_detail, None if sample_detail else ["ligand_id"])
    write_csv(out / "residue_contact_frequencies.csv", contacts_out,
              ["ligand_id", "af3_residue", "reference_residue", "contacting_samples", "analyzed_samples", "frequency"])
    write_csv(out / "candidate_hotspots.csv", [dict(af3_residue=r, reference_residue=r + args.residue_offset,
              contacting_best_models=n, analyzed_ligands=len(detail)) for r, n in hotspot.most_common()],
              ["af3_residue", "reference_residue", "contacting_best_models", "analyzed_ligands"])
    write_csv(out / "errors.csv", errors, ["ligand_id", "path", "error"])
    shortlist = [r for r in detail if r["geometry_triage_pass"]][:args.export_top]
    for i, r in enumerate(shortlist, 1):
        dest = out / "top_structures" / f"{i:03d}_{r['ligand_id']}"
        dest.mkdir(parents=True, exist_ok=True)
        for key in ("cif_path", "summary_path", "confidence_path"):
            path = Path(r[key])
            shutil.copy2(path, dest / path.name)
    write_csv(out / "shortlist.csv", shortlist, list(detail[0]) if detail else ["ligand_id"])
    (out / "settings.json").write_text(json.dumps(dict(residue_offset=args.residue_offset,
        inspected_limit=args.top, inspected_ligands=len(detail), contact_cutoff=args.contact_cutoff,
        pose_cutoff=args.pose_cutoff, exported_structures=len(shortlist), pockets=pockets,
        note="Geometry triage is not PoseBusters or affinity prediction; pose RMSD has no symmetry correction."), indent=2))
    print(f"Inspected {len(detail)} ligands, exported {len(shortlist)} structures, errors={len(errors)}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    p = subs.add_parser("probe")
    p.add_argument("--results", required=True)
    p.add_argument("--max-ligands", type=int, default=2)
    p.set_defaults(func=probe)
    p = subs.add_parser("summarize")
    p.add_argument("--results", required=True)
    p.add_argument("--jobs")
    p.add_argument("--ligands")
    p.add_argument("--out", required=True)
    p.add_argument("--samples-per-seed", type=int, default=5)
    p.add_argument("--chain-order", required=True, help="Explicit summary matrix order for older AF3; here A,B")
    p.add_argument("--protein-chain", default="A")
    p.add_argument("--ligand-chain", default="B")
    p.set_defaults(func=summarize)
    p = subs.add_parser("inspect")
    p.add_argument("--analysis", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--top", type=int, default=1000)
    p.add_argument("--export-top", type=int, default=100)
    p.add_argument("--include-ids", help="CSV with ligand_id to inspect additional docking-selected candidates")
    p.add_argument("--residue-offset", type=int, default=0, help="reference residue = AF3 label_seq_id + offset")
    p.add_argument("--pockets", help="JSON: pocket name -> reference residue integer list")
    p.add_argument("--contact-cutoff", type=float, default=4.0)
    p.add_argument("--pose-cutoff", type=float, default=2.0)
    p.set_defaults(func=inspect)
    args = parser.parse_args()
    if getattr(args, "samples_per_seed", 1) < 1 or getattr(args, "top", 1) < 1:
        parser.error("sample count and top must be positive")
    args.func(args)


if __name__ == "__main__":
    main()
