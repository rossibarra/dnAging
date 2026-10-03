#!/bin/bash -l
#SBATCH --job-name=vneCompare
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
module load conda
conda activate dnaging
cd "$PROJECT"
python scripts/compare_arg_ne_estimation_effects.py
