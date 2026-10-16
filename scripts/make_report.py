#!/usr/bin/env python3
"""Export aggregate analysis statistics without source paths or molecule identities."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary', type=Path, required=True)
    p.add_argument('--geometry', type=Path)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--synthetic', action='store_true')
    p.add_argument('--plots', action='store_true')
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    audit = json.loads((args.summary/'audit.json').read_text())
    with (args.summary/'ligand_metrics.csv').open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    values = [float(r['interface_iptm']) for r in rows if r.get('interface_iptm') and r.get('has_clash') == 'False']
    report = {key: audit[key] for key in ('expected_ligands', 'ranked_ligands', 'coverage', 'selected_samples',
              'duplicate_seed_groups_excluded', 'error_records', 'ligands_missing_interface_metric')}
    report['synthetic'] = args.synthetic
    report['nonclashing_interface_count'] = len(values)
    report['nonclashing_interface_quantiles'] = dict(zip(['p10', 'p50', 'p90'], map(float, np.quantile(values, [.1, .5, .9])))) if values else {}
    if args.geometry:
        settings = json.loads((args.geometry/'settings.json').read_text())
        report['inspected_ligands'] = settings['inspected_ligands']
        report['exported_structures'] = settings['exported_structures']
    (args.out/'aggregate.json').write_text(json.dumps(report, indent=2) + '\n')
    title = 'Synthetic demonstration — not scientific results' if args.synthetic else 'AF3 screening audit'
    lines = [f'# {title}', '', '| Metric | Value |', '| --- | --- |']
    for key, value in report.items():
        if key != 'synthetic':
            lines.append(f'| {key} | {json.dumps(value)} |')
    lines += ['', 'Coverage means at least one complete seed group per ligand. Inspect planned seed coverage separately for multi-seed studies.',
              '', 'AF3 confidence is not a measured binding affinity. Structural inspection covers the chosen subset; counts are not validated hits.', '']
    (args.out/'report.md').write_text('\n'.join(lines), encoding='utf-8')
    if args.plots and values:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        ax.hist(values, bins=40)
        ax.set(xlabel='A–B ipTM', ylabel='Nonclashing ligand states', title=title)
        fig.savefig(args.out/'interface_distribution.png', dpi=300)
        fig.savefig(args.out/'interface_distribution.pdf')
        plt.close(fig)
    print(f'Aggregate report written: {args.out}')


if __name__ == '__main__':
    main()
