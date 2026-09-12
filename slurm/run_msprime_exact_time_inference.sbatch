#!/bin/bash -l
#SBATCH --job-name=msExactInf
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-99%25

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_exact_time_ne50k}"
module load conda
conda activate dnaging
python "$PROJECT/msprime_exact_time_validation.py" infer \
    --input-root "$ROOT/simulations" \
    --frequency-table "$ROOT/constant_ne_50000.freq_table.npz" \
    --output-dir "$ROOT/exact_time_inference" \
    --replicate "$((SLURM_ARRAY_TASK_ID + 1))" \
    --n-modern 26 --epsilon 0 --mutation-age-max-tau 3
