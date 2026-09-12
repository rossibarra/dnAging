#!/bin/bash -l
#SBATCH --job-name=msInsert
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --array=0-99%25

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_exact_time_ne50k}"
module load conda
conda activate dnaging
python "$PROJECT/msprime_insertion_validation.py" infer \
  --input-root "$ROOT/simulations" --output-dir "$ROOT/insertion_inference" \
  --replicate "$((SLURM_ARRAY_TASK_ID + 1))"
