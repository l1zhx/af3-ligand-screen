#!/usr/bin/env bash
# Copy to config/cluster.local.sh; set absolute paths on your cluster.
export AF3_DIR="/path/to/alphafold3"
export AF3_PYTHON="/path/to/af3/environment/bin/python"
export MODEL_DIR="/path/to/af3/model_parameters"
export DB_DIR="/path/to/af3/databases"
# Workflow Python needs NumPy; AF3 Python can be used if NumPy is installed.
export WORKFLOW_PYTHON="$AF3_PYTHON"
# WORK_DIR is the job bundle root containing jobs/, shared/, ligands.csv.
export WORK_DIR="/path/to/private/screen"
export GPU_PARTITION="gpu-c"
export CPU_PARTITION="cpu"
export CPUS_PER_TASK=4
export MEMORY="64G"
export TIME_LIMIT="2-00:00:00"
export SAMPLES_PER_SEED=5
export FLASH_ATTENTION_IMPLEMENTATION="xla"
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE="true"
export XLA_CLIENT_MEM_FRACTION="0.90"
# Load site-required modules here. Do not set CUDA_VISIBLE_DEVICES manually.
