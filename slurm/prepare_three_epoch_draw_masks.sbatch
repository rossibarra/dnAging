#!/bin/bash -l
#SBATCH --job-name=mask3chr
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-49%20
set -euo pipefail
cd "${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
module load conda
conda activate dnaging
python prepare_three_epoch_draw_mask.py --draw "$((SLURM_ARRAY_TASK_ID+50))"
