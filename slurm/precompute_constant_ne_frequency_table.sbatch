#!/bin/bash -l
#SBATCH --job-name=constNeFreq
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00

set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PROJECT/slim_single_site_diffusion_test_A}"
mkdir -p "$OUTPUT_ROOT" "$PROJECT/logs"

module load conda
conda activate dnaging

NE_FILE="$OUTPUT_ROOT/constant_ne_10000.tsv"
if [[ ! -f "$NE_FILE" ]]; then
    printf 'series\ttime_left\ttime_right\teffective_population_size\nconstant\t0\t61000\t10000\n' > "$NE_FILE"
fi

python "$PROJECT/precompute_freq_trajectory_moments.py" \
    --ne "$NE_FILE" --ne-series constant \
    --n-sample 26 --min-n 26 \
    --t-min 0 --t-max 7000 --n-t 701 \
    --age-min 1 --age-max 61000 --n-age 250 \
    --output "$OUTPUT_ROOT/constant_ne_10000.freq_table.npz"
