#!/usr/bin/env python3
"""Build a traceable ligand manifest or a single-protein AF3 input."""
import argparse
import csv
import hashlib
import json
import sqlite3
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def msa_input(args):
    lines = args.fasta.read_text().splitlines()
    if sum(line.startswith('>') for line in lines) != 1:
        raise SystemExit('Provide exactly one FASTA record')
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>')).upper()
    if not sequence or set(sequence) - set('ACDEFGHIKLMNPQRSTVWYX'):
        raise SystemExit('Empty or unsupported protein sequence')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x', encoding='utf-8') as f:
        json.dump(dict(name='beta_catenin_msa', modelSeeds=[1], dialect='alphafold3', version=1,
                       sequences=[dict(protein=dict(id='A', sequence=sequence))]), f, indent=2)


def library(args):
    from rdkit import Chem, rdBase
    from rdkit.Chem import inchi
    files = sorted(args.sdf_dir.glob('*.sdf'))
    if not files:
        raise SystemExit('No *.sdf files found')
    args.out.mkdir(parents=True, exist_ok=False)
    counts = dict(records=0, unique=0, duplicates=0, rejected=0)
    with sqlite3.connect(args.out / 'ligands.sqlite') as db:
        db.execute('CREATE TABLE ligands (id INTEGER PRIMARY KEY, inchikey TEXT UNIQUE, smiles TEXT, source_file TEXT, source_record INTEGER, source_name TEXT)')
        db.execute('CREATE TABLE records (source_file TEXT, source_record INTEGER, ligand_id INTEGER, status TEXT, reason TEXT)')
        for path in files:
            for index, mol in enumerate(Chem.SDMolSupplier(str(path), sanitize=True, removeHs=True), 1):
                counts['records'] += 1
                try:
                    if mol is None:
                        raise ValueError('RDKit parse/sanitize failed')
                    smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                    key = inchi.MolToInchiKey(mol)
                    if not smiles or not key:
                        raise ValueError('Empty SMILES/InChIKey')
                    old = db.execute('SELECT id FROM ligands WHERE inchikey=?', (key,)).fetchone()
                    if old:
                        lid, status = old[0], 'duplicates'
                    else:
                        name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
                        lid = db.execute('INSERT INTO ligands (inchikey,smiles,source_file,source_record,source_name) VALUES (?,?,?,?,?)',
                                         (key, smiles, path.name, index, name)).lastrowid
                        status = 'unique'
                    db.execute('INSERT INTO records VALUES (?,?,?,?,?)', (path.name, index, lid, status, ''))
                    counts[status] += 1
                except (ValueError, RuntimeError) as exc:
                    counts['rejected'] += 1
                    db.execute('INSERT INTO records VALUES (?,?,?,?,?)', (path.name, index, None, 'rejected', str(exc)))
            db.commit()
            print(f'{path.name}: {counts}', flush=True)
        with (args.out / 'ligands.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['id', 'inchikey', 'smiles', 'source_file', 'source_record', 'source_name'])
            writer.writerows(db.execute('SELECT * FROM ligands ORDER BY id'))
        with (args.out / 'record_audit.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['source_file', 'source_record', 'ligand_id', 'status', 'reason'])
            writer.writerows(db.execute('SELECT * FROM records'))
    report = dict(counts, rdkit_version=rdBase.rdkitVersion,
                  inputs=[dict(file=p.name, sha256=sha(p)) for p in files],
                  policy='Standard InChIKey equality; first occurrence retained; canonical isomeric SMILES. No salt stripping or neutralization.')
    (args.out / 'preparation_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('msa-input')
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.set_defaults(func=msa_input)
    p = sub.add_parser('library')
    p.add_argument('--sdf-dir', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.set_defaults(func=library)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
