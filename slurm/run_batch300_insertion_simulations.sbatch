#!/bin/bash -l
#SBATCH --job-name=b300sim
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-299%30
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_batch300_insertion}"
module load conda
conda activate dnaging
python "$PROJECT/simulate_batch300_insertion.py" --output-root "$ROOT/simulations" --replicate "$((SLURM_ARRAY_TASK_ID + 1))"
