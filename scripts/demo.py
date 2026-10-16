#!/usr/bin/env python3
"""Create clearly synthetic AF3-like outputs and run the CPU analysis pipeline."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys


def cif_text(ligand_shift=0):
    columns = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id',
               'label_asym_id', 'label_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'pdbx_PDB_model_num']
    text = 'data_synthetic_not_a_prediction\n#\nloop_\n' + ''.join(f'_atom_site.{c}\n' for c in columns)
    atoms = [('A', 1, 'CA', (0, 0, 0)), ('A', 2, 'CA', (3, 0, 0)), ('A', 3, 'CA', (0, 3, 0)),
             ('B', 1, 'C1', (1, 1, 2.5)), ('B', 1, 'C2', (1, 2, 2.5))]
    for i, (chain, residue, name, (x, y, z)) in enumerate(atoms, 1):
        z += ligand_shift if chain == 'B' else 0
        text += f"ATOM {i} C {name} . {'ALA' if chain == 'A' else 'LIG'} {chain} {residue} {x} {y} {z} 1\n"
    return text + '#\n'


def make_fixture(root):
    root.mkdir(parents=True, exist_ok=False)
    (root / 'SYNTHETIC_DATA.txt').write_text('Artificial coordinates and scores for software validation only. Not AF3 predictions.\n')
    with (root / 'ligands.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['id', 'inchikey', 'smiles', 'source_name'])
        for i in range(1, 5):
            w.writerow([i, f'SYNTHETIC_KEY_{i}', 'CC', f'artificial_case_{i}'])
    # Complete case + higher scoring duplicate; clashing case; partial case; missing case.
    for lid, attempt, samples, iptm, clash in [(1, 100, 5, .85, False), (1, 200, 5, .99, False),
                                             (2, 300, 5, .91, True), (3, 400, 2, .8, False)]:
        folder = root / 'results/000' / f'lig_{lid:06d}' / 'seed_1' / f'attempt_{attempt}_0' / f'bcat_lig_{lid:06d}_seed_1'
        for s in range(samples):
            sample = folder / f'seed-1_sample-{s}'
            sample.mkdir(parents=True)
            data = dict(ranking_score=iptm, has_clash=clash, chain_pair_iptm=[[.9, iptm], [iptm, .7]],
                        chain_pair_pae_min=[[0, 3], [4, 0]])
            (sample / 'summary_confidences.json').write_text(json.dumps(data))
            (sample / 'model.cif').write_text(cif_text(4 if s == 4 else 0))
            confidence = dict(atom_chain_ids=['A']*3 + ['B']*2, atom_plddts=[90, 90, 90, 60, 80],
                              token_chain_ids=['A']*3 + ['B']*2, pae=[[4.0]*5 for _ in range(5)])
            (sample / 'confidences.json').write_text(json.dumps(confidence))
        (folder / f'bcat_lig_{lid:06d}_seed_1_summary_confidences.json').write_text(json.dumps(data))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    root = args.out.resolve()
    make_fixture(root)
    script = Path(__file__).with_name('af3_postprocess.py')
    subprocess.run([sys.executable, str(script), 'summarize', '--results', str(root/'results'),
                    '--ligands', str(root/'ligands.csv'), '--chain-order', 'A,B', '--out', str(root/'summary')], check=True)
    subprocess.run([sys.executable, str(script), 'inspect', '--analysis', str(root/'summary'),
                    '--out', str(root/'geometry'), '--top', '4', '--export-top', '100'], check=True)
    audit = json.loads((root/'summary/audit.json').read_text())
    assert audit['ranked_ligands'] == 2 and audit['selected_samples'] == 10
    assert audit['duplicate_seed_groups_excluded'] == 1
    assert audit['coverage'] == dict(complete_sample_group=2, partial_or_invalid=1, missing=1)
    settings = json.loads((root/'geometry/settings.json').read_text())
    assert settings['exported_structures'] == 1
    print(f'DEMO PASS: 4 artificial cases; 2 ranked; 1 duplicate excluded; 1 structure exported. {root}')


if __name__ == '__main__':
    main()
