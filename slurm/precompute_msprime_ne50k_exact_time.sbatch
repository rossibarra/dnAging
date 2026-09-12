#!/bin/bash -l
#SBATCH --job-name=msExactTab
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_exact_time_ne50k}"
mkdir -p "$ROOT"
module load conda
conda activate dnaging
NE_FILE="$ROOT/constant_ne_50000.tsv"
if [[ ! -f "$NE_FILE" ]]; then
    printf 'series\ttime_left\ttime_right\teffective_population_size\nconstant\t0\t301000\t50000\n' > "$NE_FILE"
fi
python "$PROJECT/precompute_freq_trajectory_moments.py" \
    --ne "$NE_FILE" --ne-series constant \
    --n-sample 26 --min-n 26 \
    --t-min 0 --t-max 10000 --n-t 401 \
    --age-min 1 --age-max 301000 --n-age 250 \
    --output "$ROOT/constant_ne_50000.freq_table.npz"
