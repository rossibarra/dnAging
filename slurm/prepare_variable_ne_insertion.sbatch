#!/bin/bash -l
#SBATCH --job-name=vneInsPrep
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-10%5
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
module load conda
conda activate dnaging
python "$PROJECT/validation/prepare_variable_ne_insertion.py" --root "$PROJECT/msprime_variable_ne_error" \
  --replicate "$SLURM_ARRAY_TASK_ID" --keep-last "${KEEP_LAST:-50}" \
  --output-root "$PROJECT/msprime_variable_ne_error/insertion_posterior"
