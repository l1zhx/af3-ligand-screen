#!/usr/bin/env python3
"""Build compact AF3 protein-ligand jobs from precomputed monomer data and CSV.

AF3 input format version 2 refers to one shared A3M and shared template CIFs.
Run on the preparation machine or cluster login node. No RDKit dependency is needed because SMILES are in CSV.
"""
import argparse
import csv
import hashlib
import json
import re
import zipfile
from pathlib import Path


def json_bytes(obj):
    return (json.dumps(obj, ensure_ascii=False, separators=(",", ":")) + "\n").encode("utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--msa-json", type=Path, required=True)
    parser.add_argument("--ligands-csv", type=Path, required=True)
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--out-dir", type=Path, help="write extracted jobs/shared files")
    output.add_argument("--archive", type=Path, help="write a ZIP bundle")
    parser.add_argument("--seeds", default="1", help="comma-separated model seeds, default 1")
    parser.add_argument("--round", default="seed_1", help="label appended to each unique job name")
    parser.add_argument("--example", type=Path, help="also save the first ligand JSON for inspection")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_]+", args.round):
        parser.error("--round must contain only letters, digits, underscore")
    try:
        seeds = [int(x) for x in args.seeds.split(",")]
    except ValueError:
        parser.error("--seeds must be comma-separated integers")
    if not seeds or len(set(seeds)) != len(seeds) or any(x < 0 for x in seeds):
        parser.error("--seeds must be distinct nonnegative integers")

    original = json.loads(args.msa_json.read_text(encoding="utf-8"))
    proteins = [entry["protein"] for entry in original["sequences"] if "protein" in entry]
    if len(proteins) != 1:
        raise SystemExit("MSA JSON must contain exactly one protein")
    p = proteins[0]
    if p.get("id") != "A" or not p.get("sequence") or not p.get("unpairedMsa"):
        raise SystemExit("MSA JSON must contain protein chain A and nonempty unpairedMsa")
    lines = p["unpairedMsa"].splitlines()
    if not lines or not lines[0].startswith(">"):
        raise SystemExit("unpairedMsa must start with a FASTA/A3M header")
    query_lines = []
    for line in lines[1:]:
        if line.startswith(">"):
            break
        query_lines.append(line.strip())
    query = "".join(query_lines)
    if query != p["sequence"]:
        raise SystemExit("unpairedMsa query differs from protein sequence")
    templates = p.get("templates")
    if not isinstance(templates, list):
        raise SystemExit("Precomputed template list is missing")

    shared = {"shared/beta_catenin_unpaired.a3m": p["unpairedMsa"].encode("utf-8")}
    template_refs = []
    for i, t in enumerate(templates, 1):
        cif = t.get("mmcif")
        if not cif or not isinstance(t.get("queryIndices"), list) or not isinstance(t.get("templateIndices"), list):
            raise SystemExit(f"Template {i} lacks mmcif or residue mappings")
        name = f"shared/template_{i:02d}.cif"
        shared[name] = cif.encode("utf-8")
        template_refs.append({
            "mmcifPath": f"../../{name}",
            "queryIndices": t["queryIndices"],
            "templateIndices": t["templateIndices"],
        })

    base_protein = {
        "id": "A", "sequence": p["sequence"],
        "unpairedMsaPath": "../../shared/beta_catenin_unpaired.a3m",
        "pairedMsa": "",
        "templates": template_refs,
    }
    if p.get("modifications"):
        base_protein["modifications"] = p["modifications"]

    if args.out_dir:
        root = args.out_dir
        if any((root / name).exists() for name in ("jobs", "shared", "bundle_meta.json", "ligands.csv")):
            raise SystemExit(f"Job bundle already exists under {root}; choose a fresh directory")
        root.mkdir(parents=True, exist_ok=True)

        def write(name, data):
            dest = root / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
    else:
        args.archive.parent.mkdir(parents=True, exist_ok=True)
        if args.archive.exists():
            raise SystemExit(f"Archive already exists: {args.archive}")
        archive = zipfile.ZipFile(args.archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True)

        def write(name, data):
            archive.writestr(name, data)

    try:
        for name, data in shared.items():
            write(name, data)
        write("ligands.csv", args.ligands_csv.read_bytes())
        count = 0
        seen_ids = set()
        seen_keys = set()
        example = None
        with args.ligands_csv.open(newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            required = {"id", "inchikey", "smiles"}
            if not required.issubset(reader.fieldnames or []):
                raise SystemExit(f"CSV must contain {sorted(required)}")
            for row in reader:
                ligand_id = int(row["id"])
                if not 1 <= ligand_id <= 999999 or ligand_id in seen_ids or row["inchikey"] in seen_keys or not row["inchikey"] or not row["smiles"]:
                    raise SystemExit(f"Invalid/duplicate ligand ID, InChIKey, or SMILES at row {count+2}")
                seen_ids.add(ligand_id)
                seen_keys.add(row["inchikey"])
                count += 1
                name = f"bcat_lig_{ligand_id:06d}_{args.round}"
                job = {
                    "name": name,
                    "modelSeeds": seeds,
                    "sequences": [
                        {"protein": base_protein},
                        {"ligand": {"id": "B", "smiles": row["smiles"]}},
                    ],
                    "dialect": "alphafold3", "version": 2,
                }
                payload = json_bytes(job)
                shard = (ligand_id - 1) // 1000
                write(f"jobs/{shard:03d}/{name}.json", payload)
                if example is None:
                    example = payload
        if count == 0 or seen_ids != set(range(1, count + 1)):
            raise SystemExit("Ligand IDs must be contiguous from 1")
        meta = {
            "count": count, "round": args.round, "modelSeeds": seeds,
            "sourceMsaSha256": hashlib.sha256(args.msa_json.read_bytes()).hexdigest(),
            "sourceCsvSha256": hashlib.sha256(args.ligands_csv.read_bytes()).hexdigest(),
        }
        write("bundle_meta.json", json_bytes(meta))
    finally:
        if args.archive:
            archive.close()

    if args.example:
        args.example.parent.mkdir(parents=True, exist_ok=True)
        args.example.write_bytes(example)
    print(f"Generated {count} AF3 jobs; round={args.round}; seeds={seeds}; shared templates={len(template_refs)}")
    print(f"Output: {args.out_dir or args.archive}")


if __name__ == "__main__":
    main()
