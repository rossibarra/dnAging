#!/bin/bash -l
#SBATCH --job-name=true3chr
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=24:00:00
set -euo pipefail
cd "${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
module load conda
conda activate dnaging
python three_epoch_true_arg_insertion.py
