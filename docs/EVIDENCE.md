# Evidence register

This page distinguishes campaign observations from software demonstration and future validation. It does not contain private structures, ligand identities, internal paths or cluster account information.

| Claim | Available evidence | What remains unverified here |
| --- | --- | --- |
| Campaign scale: 59,635 ligand entries | Manifest size reported by the original analysis audit | Full production manifest is private and not distributed |
| Inference campaign finished | Operator reported completion across result shards 000–059 | A directory or job completion alone does not establish valid samples for all entries |
| Real AF3 output parsing works for a representative ligand | User-provided command output on 2026-10-15: 1 expected and ranked ligand, 5 selected samples, 0 errors, 0 missing interface metrics | Does not establish whole-library success rate or structural correctness |
| Root best copy is excluded | Real probe discovered 6 summaries and recognized 5 sample summaries | Full-library duplicate statistics remain pending |
| Portable workflow functionality | Synthetic tests and demo in this repository | Refactored Slurm scripts need a four-ligand test on the target cluster |
| Whole-library hit rates, performance and docking agreement | No audited public result supplied yet | Must be computed from the actual full-library outputs and accounting |
| Biological efficacy | No experimental validation supplied | Binding, inhibition, specificity and cellular activity remain open |

## Real one-ligand parser check

The following scalar values are transcribed from the operator's successful probe; paths have been omitted:

```json
{
  "expected_ligands": 1,
  "ranked_ligands": 1,
  "coverage": {"complete_sample_group": 1},
  "selected_samples": 5,
  "duplicate_seed_groups_excluded": 0,
  "error_records": 0,
  "ligands_missing_interface_metric": 0,
  "discovered_summary_files": 6,
  "recognized_sample_summary_files": 5,
  "samples_per_seed": 5
}
```

## Reproducing software evidence

```bash
python -m unittest discover -s tests -v
python scripts/demo.py --out demo_output
python scripts/make_report.py --summary demo_output/summary \
  --geometry demo_output/geometry --synthetic --out demo_output/report
```

The demo intentionally includes missing, partial, clashing and duplicate cases. Its output counts are software expectations, not a retrospective estimate of the real screen. [Validation record](VALIDATION.md) describes the local checks and environmental limitations. GitHub Actions runs Linux tests after this repository is uploaded; a configured workflow is not evidence that CI has already passed.

## Promote actual analysis to a portfolio result

After the real CPU analysis completes, review coverage, error reasons, repeated attempts and planned seed coverage. Generate an aggregate report using `make_report.py` without `--synthetic`. Add only data approved for public sharing and update this register with the date, workflow commit, AF3 release, input hashes, number of valid molecules and sampled structures, geometry-review denominator and unresolved failures.

For performance, record GPU model, resource requests, ligand/token size distribution, seeds, samples, cold/warm compilation, median and p90 per-ligand time, completed structures per GPU-hour and failure/retry rate. Slurm elapsed time is wall time for a task; summing GPU allocations gives allocated GPU-hours, not measured utilization. Do not extrapolate a single 13-minute observation into a claimed whole-campaign benchmark.
