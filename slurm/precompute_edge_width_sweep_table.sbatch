#!/bin/bash -l
#SBATCH --job-name=ewTable
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=20:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/edge_width_ne_sweep}"
NE="${NE:?set NE}"
# tau = 3 in generations, plus a margin so the cutoff is inside the grid.
AGE_MAX=$(( NE * 6 + 1000 ))
mkdir -p "$ROOT/tables"
module load conda
conda activate dnaging
NE_FILE="$ROOT/tables/constant_ne_${NE}.tsv"
printf 'series\ttime_left\ttime_right\teffective_population_size\nconstant\t0\t%s\t%s\n' \
    "$AGE_MAX" "$NE" > "$NE_FILE"
# T grid widened to 25,000: the Ne=50K run censored 22%% of the 30-100K stratum
# at its 10,000 ceiling, so the measured shift there was a lower bound.
python "$PROJECT/precompute_freq_trajectory_moments.py" \
    --ne "$NE_FILE" --ne-series constant \
    --n-sample 26 --min-n 26 \
    --t-min 0 --t-max 25000 --n-t 501 \
    --age-min 1 --age-max "$AGE_MAX" --n-age 250 \
    --output "$ROOT/tables/constant_ne_${NE}.freq_table.npz"
