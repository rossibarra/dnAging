#!/bin/bash -l
#SBATCH --job-name=varNeIns
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --array=0-999%30

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_variable_ne_error}"
REPLICATE="$((SLURM_ARRAY_TASK_ID / 10 + 1))"
ANCIENT_INDEX="$((SLURM_ARRAY_TASK_ID % 10 + 1))"
ANCIENT_SAMPLE="$(printf 'ancient_%02d' "$ANCIENT_INDEX")"
module load conda
conda activate dnaging
python "$PROJECT/msprime_insertion_validation.py" infer \
  --input-root "$ROOT/simulations" \
  --output-dir "$ROOT/insertion_inference/$ANCIENT_SAMPLE" \
  --replicate "$REPLICATE" --n-modern 26 \
  --ancient-sample "$ANCIENT_SAMPLE" \
  --ne-table "$ROOT/simulations/$(printf 'replicate_%03d' "$REPLICATE")/constant_ne_epochs.tsv" \
  --ne-series variable --ne-extrapolate-last --epsilon "${EPSILON:-0.01}" \
  --age-min 0 --age-max 10000 --age-step 20
