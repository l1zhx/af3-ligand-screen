# Cluster runbook

Run commands from the repository root on a Linux cluster. AF3 and its model parameters must already be installed for your GPU architecture. Preparation and analysis use the dependencies here; the AF3 runtime remains in its own environment. Verify the installed AF3 version against its own documentation, since upstream flags can change.

## 1. Configure the cluster

```bash
cp config/cluster.example.sh config/cluster.local.sh
# Edit the local file: AF3_DIR, AF3_PYTHON, MODEL_DIR, DB_DIR, WORK_DIR and partitions.
source config/cluster.local.sh
export CLUSTER_CONFIG="$PWD/config/cluster.local.sh"
mkdir -p "$WORK_DIR" "$WORK_DIR/logs"
test -f "$AF3_DIR/run_alphafold.py"
test -x "$AF3_PYTHON"
test -d "$MODEL_DIR"
"$WORKFLOW_PYTHON" -c 'import numpy; print(numpy.__version__)'
sinfo -o '%P %G %C'
```

All configured paths must be visible on compute nodes. `WORK_DIR` is your data directory; `AF3_DIR` is the software installation. Slurm opens log files **before** running the script, so create the log directory before submission. A partition listing GPUs does not necessarily require requesting a GPU; use the site's CPU queue for postprocessing.

Record deployment information privately:

```bash
mkdir -p "$WORK_DIR/provenance"
"$AF3_PYTHON" -m pip freeze > "$WORK_DIR/provenance/af3_environment.txt"
"$WORKFLOW_PYTHON" -m pip freeze > "$WORK_DIR/provenance/analysis_environment.txt"
git rev-parse HEAD > "$WORK_DIR/provenance/workflow_commit.txt"
git -C "$AF3_DIR" rev-parse HEAD > "$WORK_DIR/provenance/af3_commit.txt"
# If AF3 was not installed from Git, record the installed release manually instead.
"$AF3_PYTHON" "$AF3_DIR/run_alphafold.py" --helpfull > "$WORK_DIR/provenance/af3_help.txt" 2>&1
```

Some flag libraries return a nonzero exit code after printing help. Inspect the saved text. This workflow uses `run_data_pipeline=false` and `flash_attention_implementation=xla`; it does not pass the previously incompatible `resolve_msa_overlaps` flag. `SAMPLES_PER_SEED` is the **expected output count**, not a flag changing the model configuration. Confirm the runtime's diffusion sample setting; default examples expect 5.

## 2. Prepare the library and compute the receptor MSA once

```bash
python -m pip install -r requirements-chem.txt
python scripts/prepare.py library --sdf-dir /private/input_sdf --out data/prepared
python scripts/prepare.py msa-input --fasta examples/beta_catenin.fasta --out data/receptor.json

export MSA_INPUT="$PWD/data/receptor.json"
export MSA_OUTPUT="$WORK_DIR/msa_output"
sbatch --partition="$CPU_PARTITION" \
  --output="$WORK_DIR/logs/msa-%j.out" slurm/msa.slurm
```

The MSA job runs `--run_inference=false`. Use the resulting `*_data.json` containing the expanded protein `unpairedMsa` and template list. If MSA computation occurs on another machine, transfer that JSON and the deduplicated `ligands.csv` to the cluster. Model parameters and sequence databases remain in the cluster's existing AF3 installation.

Generate the bundle (adjust the MSA JSON filename to the actual output):

```bash
python scripts/build_jobs.py \
  --msa-json "$MSA_OUTPUT/beta_catenin_msa/beta_catenin_msa_data.json" \
  --ligands-csv data/prepared/ligands.csv \
  --out-dir "$WORK_DIR" --seeds 1 --round seed_1
```

The generator accepts an existing work directory but refuses to overwrite `jobs/`, `shared/`, `ligands.csv` or bundle metadata. It requires one chain-A protein, a query-matching nonempty unpaired MSA, an explicit template list (possibly empty), and contiguous CSV IDs from 1. Generated input version 2 requires AF3 support for `unpairedMsaPath` and template `mmcifPath`. Relative paths are resolved from each JSON, so transfer the entire `jobs/` and `shared/` tree together.

```text
WORK_DIR/
  bundle_meta.json
  ligands.csv
  shared/beta_catenin_unpaired.a3m
  shared/template_01.cif            # if templates are present
  jobs/000/bcat_lig_000001_seed_1.json
  jobs/001/bcat_lig_001001_seed_1.json
  results/                         # created by inference
  logs/
```

## 3. Four-ligand smoke test

```bash
"$WORKFLOW_PYTHON" scripts/workflow.py plan \
  --work "$WORK_DIR" --start-id 1 --end-id 4 \
  --chunk-size 1 --batch-tasks 4 --out "$WORK_DIR/plans/test4"
"$WORKFLOW_PYTHON" scripts/workflow.py submit \
  --plan "$WORK_DIR/plans/test4/plan.json" --max-active 4 --max-batches 1 --dry-run
"$WORKFLOW_PYTHON" scripts/workflow.py submit \
  --plan "$WORK_DIR/plans/test4/plan.json" --max-active 4 --max-batches 1
```

`--dry-run` prints arguments and does not submit. The real command calls `sbatch slurm/infer.slurm`, with one GPU per task. Four tasks can run on different nodes. The submission limit is a maximum, not a reservation or a guarantee of four simultaneous GPUs. CPU, memory, account limits, priority and available GPUs must all permit scheduling.

Before the full screen, inspect the logs and `complete.json` for all four inputs. Audit the four-result directory with `summarize` and inspect representative structures. The new runner validates every requested seed and sample, not merely AF3's process exit code or a root-level best structure.

## 4. Submit the library and continue without duplicating arrays

```bash
"$WORKFLOW_PYTHON" scripts/workflow.py plan \
  --work "$WORK_DIR" --round seed_1 \
  --chunk-size 4 --batch-tasks 500 --out "$WORK_DIR/plans/full"
"$WORKFLOW_PYTHON" scripts/workflow.py submit \
  --plan "$WORK_DIR/plans/full/plan.json" --max-active 12 --max-batches 5
```

This is the replacement for the legacy `MAX_ACTIVE=12 MAX_BATCHES=5 bash submit_gpu_c.sh full`. One array covers up to 2,000 ligands; five arrays cover up to 10,000. Arrays use `afterany` dependencies, so only one array in this chain is eligible at a time, keeping the chain's concurrency bounded. Failures remain visible while later arrays can proceed. Near the end of an array, the chain may temporarily underutilize GPUs while its last tasks finish; this is a deliberate simple scheduling tradeoff.

To submit the next five arrays, repeat the **same submit command using the same plan**. `submission.json` records accepted IDs; already submitted arrays are not resubmitted. `--max-batches 0` requests all remaining arrays, subject to site submission limits. The helper stops on the first error and records its uncertain submission state rather than retrying blindly.

To increase concurrency for arrays already submitted, use the actual IDs in `submission.json`:

```bash
scontrol update JobId=YOUR_ARRAY_ID ArrayTaskThrottle=24
```

Change the IDs individually for queued successor arrays as needed. Changing `--max-active` affects only future submissions. If your partition has 40 GPUs, a throttle of 40 permits up to 40 of your array tasks, but Slurm still respects other jobs, quotas and CPU/memory availability. Do not manually set GPU device IDs. A chunk of 4 ligands runs sequentially on one GPU; it is not four concurrent GPU jobs.

With one GPU per ligand, an idealized runtime is `remaining_ligands × measured_minutes_per_ligand / (concurrent_GPUs × 60)` hours. This excludes queue time, failures, compilation changes and array-tail idle time. Use representative measurements rather than a single ligand's timing.

## 5. Monitor, recover and retry

```bash
squeue -u "$USER" -n bcat-af3 -o '%.25i %.10T %.12M %.40R'
sacct -j YOUR_ARRAY_ID -X --format=JobID,State,Elapsed,ExitCode,AllocTRES
tail -n 50 "$WORK_DIR/logs/af3-YOUR_ARRAY_ID_0.out"
```

`JobArrayTaskLimit` means the current array throttle is reached. `Dependency` means a predecessor array is unfinished. `Resources` means scheduling requirements are currently unavailable. `No error` from `scontrol` indicates a successful update, not an inference failure.

For an **uncertain submission**, inspect the last `inflight.batch` in `submission.json`. Use `squeue`/`sacct` around its timestamp and `scontrol show job JOB_ID`; `scontrol write batch_script JOB_ID -` can help check the submitted script. Correlate the output log path, submission time, task range and exported plan before deciding. After determining whether Slurm accepted it, back up `submission.json`. If accepted, append the confirmed `{ "batch": N, "job_id": "ID" }` to `accepted` in order; if definitely rejected, leave `accepted` unchanged. Then set `inflight` to `null` and resume. Do not clear the state while acceptance is unknown.

Dependencies on very old job IDs can expire from Slurm's retained history. If a resume reports an invalid dependency, verify the predecessor's terminal state with accounting and create a fresh plan covering only the unsubmitted ligand range. Use `--start-id`/`--end-id` from the old plan, not guessed shard counts.

For a failed ligand, wait for its old task to exit, inspect the error, correct the environment, then create a new small retry plan:

```bash
"$WORKFLOW_PYTHON" scripts/workflow.py plan --work "$WORK_DIR" \
  --start-id 1708 --end-id 1708 --out "$WORK_DIR/plans/retry_1708"
"$WORKFLOW_PYTHON" scripts/workflow.py submit \
  --plan "$WORK_DIR/plans/retry_1708/plan.json" --max-active 1 --max-batches 1
```

Attempt folders are immutable. A requeue of the same Slurm task will refuse an existing attempt directory; submit a new retry task instead. A completed record is reused only when input hashes match and its files still pass the sample check. Old `complete` sentinel files from earlier scripts are not trusted by this runner. For an old campaign, use the analysis entry point rather than resubmitting its completed inputs.

Per-ligand Linux `flock` prevents simultaneous workers in the same work directory. Verify advisory lock support on your shared filesystem. It does not coordinate different copies of WORK_DIR or old scripts that do not acquire this lock. Do not run both submission systems for the same inputs. Completion reuse assumes the same scientific protocol; for changed model parameters, sample settings or intended replicate experiments, use a new round and bundle.

## 6. Audit existing or new results on the CPU queue

```bash
"$WORKFLOW_PYTHON" scripts/af3_postprocess.py probe \
  --results "$WORK_DIR/results/000/lig_000001"
"$WORKFLOW_PYTHON" scripts/af3_postprocess.py summarize \
  --results "$WORK_DIR/results/000/lig_000001" --chain-order A,B \
  --out "$WORK_DIR/probe_summary_v1"

export ANALYSIS_DIR="$WORK_DIR/analysis_v1"
# Set 130 only for the verified 556-aa beta-catenin construct; see METHODS.md.
export RESIDUE_OFFSET=130
sbatch --partition="$CPU_PARTITION" \
  --output="$WORK_DIR/logs/af3-analysis-%j.out" slurm/postprocess.slurm
```

Always choose a fresh analysis directory. For a completed analysis job, check elapsed time and exit status with `sacct`, not the log filename. Startup and progress messages are unbuffered. Summary scanning is serial and mostly small-file I/O; a higher CPU allocation alone will not speed that loop. A directory's existence is not proof of successful output.

The real one-ligand smoke test found 6 summaries: 5 sample summaries plus one root-level best copy. Correct counts were `ranked_ligands=1`, `selected_samples=5`, `error_records=0`, and no missing A–B interface metric. For the full campaign check `audit.json`, `errors.csv`, `coverage.csv` and all expected groups before interpreting rankings.

## 7. Compare with existing docking analyses

Move only the needed AF3 tables to the machine with RDKit and your docking cluster table:

```bash
python scripts/merge_af3_glide.py \
  --af3-csv /private/analysis/summary/ligand_metrics.csv \
  --geometry-csv /private/analysis/geometry/ligand_geometry.csv \
  --glide-csv /private/docking/cluster_assignments.csv \
  --glide-selection-csv /private/docking/diverse_top60_for_review.csv \
  --out /private/comparison_v1 --plots
```

`--glide-selection-csv` is optional. Its entries must have `parent_smiles`. Required cluster columns are `parent_smiles`, `docking_score`, `cluster_id`; optional annotations such as cluster size and PAINS are retained. The clustering itself is an external input, not rerun by this repository.

The comparison emits `docking_selected_af3_ids.csv`. If important docking candidates were outside the AF3 Top 1,000, transfer that ID file back to the AF3 storage machine and run a new inspection:

```bash
"$WORKFLOW_PYTHON" scripts/af3_postprocess.py inspect \
  --analysis "$ANALYSIS_DIR/summary" --include-ids /private/docking_selected_af3_ids.csv \
  --top 1000 --export-top 100 --residue-offset 130 \
  --out "$ANALYSIS_DIR/geometry_with_docking"
```

The inspection needs the original CIF and confidence paths recorded in the summary. Run it where those paths exist. Rerun the comparison with the expanded geometry table and a new comparison output directory. Report unmapped AF3 molecules separately rather than discarding them.

## 8. Generate a shareable aggregate report

```bash
python scripts/make_report.py --summary "$ANALYSIS_DIR/summary" \
  --geometry "$ANALYSIS_DIR/geometry" --out /private/public_report_v1 --plots
```

This produces aggregate JSON/Markdown and an optional interface distribution figure; source paths and molecular identities are omitted. Only add reviewed aggregate artifacts to GitHub. Large confidence matrices, MSA databases, weights and private ligand files stay in private storage.
