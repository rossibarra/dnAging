#!/bin/bash -l
#SBATCH --job-name=b300ins
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --array=0-299%30
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_batch300_insertion}"
module load conda
conda activate dnaging
NMODERN=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["n_modern_haploid"])' "$ROOT/simulations/replicate_$(printf '%03d' "$((SLURM_ARRAY_TASK_ID + 1))")/metadata.json")
python "$PROJECT/msprime_insertion_validation.py" infer --input-root "$ROOT/simulations" --output-dir "$ROOT/insertion_inference" --replicate "$((SLURM_ARRAY_TASK_ID + 1))" --n-modern "$NMODERN" --age-min 0 --age-max 15000 --age-step 50
