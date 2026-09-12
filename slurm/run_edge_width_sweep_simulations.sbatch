#!/bin/bash -l
#SBATCH --job-name=ewSim
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=0-99%25

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/edge_width_ne_sweep}"
NE="${NE:?set NE}"
module load conda
conda activate dnaging
# Identical to the Ne=50K control in every respect but Ne, so the contrast
# across Ne varies one thing.
python "$PROJECT/msprime_exact_time_validation.py" simulate \
    --output-root "$ROOT/simulations_ne${NE}" \
    --replicate "$((SLURM_ARRAY_TASK_ID + 1))" \
    --base-seed "${BASE_SEED:?set BASE_SEED}" \
    --ne "$NE" --n-modern 26 --length 10000000 \
    --mutation-rate 1e-8 --recombination-rate 1e-8 \
    --age-min 0 --age-max 10000
