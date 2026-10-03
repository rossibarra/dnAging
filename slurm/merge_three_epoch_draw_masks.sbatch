#!/bin/bash -l
#SBATCH --job-name=maskMerge3
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
cd "${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
module load conda
conda activate dnaging
python validation/merge_three_epoch_draw_masks.py
