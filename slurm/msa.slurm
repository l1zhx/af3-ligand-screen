#!/usr/bin/env bash
#SBATCH --job-name=bcat-msa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
set -euo pipefail
REPO="${SLURM_SUBMIT_DIR:?Submit from repository root}"
source "${CLUSTER_CONFIG:-$REPO/config/cluster.local.sh}"
: "${MSA_INPUT:?Set MSA_INPUT to an absolute protein-only JSON path}"
: "${MSA_OUTPUT:?Set MSA_OUTPUT to an absolute output directory}"
cd "$AF3_DIR"
exec "$AF3_PYTHON" -u run_alphafold.py --json_path="$MSA_INPUT" \
  --db_dir="$DB_DIR" --output_dir="$MSA_OUTPUT" --run_inference=false
