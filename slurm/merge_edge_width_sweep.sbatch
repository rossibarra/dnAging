#!/bin/bash -l
#SBATCH --job-name=ewMerge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/edge_width_ne_sweep}"
NE="${NE:?set NE}"
module load conda
conda activate dnaging
python "$PROJECT/scripts/edge_width_scaling.py" merge \
    --input-dir "$ROOT/inference_ne${NE}" \
    --output-dir "$ROOT/results_ne${NE}" \
    --width-edges 3000 10000 30000 100000 300000
