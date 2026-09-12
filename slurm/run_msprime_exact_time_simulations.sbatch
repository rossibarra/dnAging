#!/bin/bash -l
#SBATCH --job-name=msExactSim
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-99%25

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PROJECT/msprime_exact_time_ne50k/simulations}"
module load conda
conda activate dnaging
python "$PROJECT/msprime_exact_time_validation.py" simulate \
    --output-root "$OUTPUT_ROOT" \
    --replicate "$((SLURM_ARRAY_TASK_ID + 1))" \
    --base-seed "${BASE_SEED:-1700000000}" \
    --ne 50000 --n-modern 26 --length 10000000 \
    --mutation-rate 1e-8 --recombination-rate 1e-8 \
    --age-min 0 --age-max 10000
