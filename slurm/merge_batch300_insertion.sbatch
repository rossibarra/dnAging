#!/bin/bash -l
#SBATCH --job-name=b300merge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_batch300_insertion}"
module load conda
conda activate dnaging
python "$PROJECT/msprime_insertion_validation.py" merge --input-dir "$ROOT/insertion_inference" --output-dir "$ROOT/insertion_results" --n-replicates 300
