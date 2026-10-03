#!/bin/bash -l
#SBATCH --job-name=vneSim
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=2-100%25

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_variable_ne_error}"
module load conda
conda activate dnaging
# Replicate 1 was produced as the smoke test; the script refuses to overwrite an
# existing replicate directory, so reruns are safe.
python "$PROJECT/simulate_variable_ne_error.py" \
    --output-root "$ROOT/simulations" \
    --replicate "$SLURM_ARRAY_TASK_ID" \
    --base-seed "${BASE_SEED:-1730000000}" \
    --n-epochs 10 --ne-min 5000 --ne-max 150000 --ne-prior loguniform \
    --n-modern 26 --n-ancient 10 --age-min 0 --age-max 10000 \
    --length 10000000 --mutation-rate 1e-8 --recombination-rate 1e-8 \
    --epsilon 0.01
