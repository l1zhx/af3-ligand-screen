#!/usr/bin/env bash
#SBATCH --job-name=bcat-analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
set -euo pipefail
trap 'rc=$?; echo "$(date -Is) finished exit=$rc elapsed=${SECONDS}s"' EXIT
REPO="${SLURM_SUBMIT_DIR:?Submit from repository root}"
source "${CLUSTER_CONFIG:-$REPO/config/cluster.local.sh}"
: "${ANALYSIS_DIR:?Set a new absolute ANALYSIS_DIR}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
echo "$(date -Is) starting analysis on $(hostname)"
"$WORKFLOW_PYTHON" -u "$REPO/scripts/af3_postprocess.py" summarize \
  --results "$WORK_DIR/results" --jobs "$WORK_DIR/jobs" --ligands "$WORK_DIR/ligands.csv" \
  --samples-per-seed "${SAMPLES_PER_SEED:-5}" --chain-order A,B --out "$ANALYSIS_DIR/summary"
"$WORKFLOW_PYTHON" -u "$REPO/scripts/af3_postprocess.py" inspect \
  --analysis "$ANALYSIS_DIR/summary" --out "$ANALYSIS_DIR/geometry" \
  --top "${TOP_CANDIDATES:-1000}" --export-top "${EXPORT_TOP:-100}" \
  --residue-offset "${RESIDUE_OFFSET:-0}"
