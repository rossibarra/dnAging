#!/bin/bash -l
#SBATCH --job-name=insertMerge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
set -euo pipefail
: "${MANIFEST:?set MANIFEST to the draw manifest}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTROOT="${OUTROOT:-$PROJECT/insertion_arg_draws}"
module load conda
conda activate dnaging
python "$PROJECT/insertion_real_data.py" merge --manifest "$MANIFEST" \
  --draw-root "$OUTROOT" --output "$OUTROOT/merged"
