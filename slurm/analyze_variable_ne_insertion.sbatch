#!/bin/bash -l
#SBATCH --job-name=varNeAnalyze
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_variable_ne_error}"
module load conda
conda activate dnaging
python "$PROJECT/analyze_variable_ne_insertion.py" --root "$ROOT" --output-dir "$ROOT/insertion_results"
