#!/bin/bash -l
#SBATCH --job-name=msEdgeMerge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_exact_time_ne50k}"
module load conda
conda activate dnaging
python "$PROJECT/msprime_exact_time_validation.py" merge \
    --input-dir "$ROOT/edge_interval_uniform_inference" \
    --output-dir "$ROOT/edge_interval_uniform_results" --n-replicates 100
