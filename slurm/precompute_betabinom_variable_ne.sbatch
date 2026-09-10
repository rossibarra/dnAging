#!/bin/bash
#SBATCH --job-name=TbetaNe
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
python -c 'import numpy, mpmath' 2>/dev/null || {
    echo "ERROR: activate the dnaging conda environment before submission" >&2
    exit 90
}
python "$PROJECT/precompute_betabinom_variable_ne.py" \
    --ne "$NE_FILE" --ne-series posterior_mean --n-panel 26 \
    --u-min 1e-9 --u-max 3 --n-u 80 --dps 360 --output "$OUTPUT"
