#!/bin/bash -l
#SBATCH --job-name=linvinf
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-3

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="$PROJECT/msprime_lineage_variance"
case "$SLURM_ARRAY_TASK_ID" in
  1) AGE=2000 ;; 2) AGE=5000 ;; 3) AGE=8000 ;;
esac
module load conda
conda activate dnaging
python "$PROJECT/lineage_variance_infer.py" \
    --replicate-dir "$ROOT/simulations_age${AGE}/replicate_001" \
    --output "$ROOT/results/age${AGE}.json" \
    --n-blocks 50 --n-boot 1000 --age-min 0 --age-max 15000 --age-step 50
