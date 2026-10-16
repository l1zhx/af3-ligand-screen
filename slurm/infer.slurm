#!/usr/bin/env bash
#SBATCH --job-name=bcat-af3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
set -euo pipefail
REPO="${SLURM_SUBMIT_DIR:?Submit from repository root}"
source "${CLUSTER_CONFIG:-$REPO/config/cluster.local.sh}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export PYTHONUNBUFFERED=1
# Slurm assigns exactly one GPU per array task. This task executes a chunk serially.
exec "$WORKFLOW_PYTHON" -u "$REPO/scripts/workflow.py" run \
  --plan "${PLAN:?Missing PLAN}" --batch "${BATCH_INDEX:?Missing BATCH_INDEX}" \
  --task "${SLURM_ARRAY_TASK_ID:?Requires Slurm array}"
