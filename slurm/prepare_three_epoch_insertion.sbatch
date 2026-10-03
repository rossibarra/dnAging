#!/bin/bash -l
#SBATCH --job-name=prep3chr
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=08:00:00
set -euo pipefail
cd "${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
module load conda
conda activate dnaging
python prepare_three_epoch_insertion.py
