#!/bin/bash -l
#SBATCH --job-name=vneInsDraw
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-499%50
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
SCENARIO="${SCENARIO:?set SCENARIO to true_ne or estimated_ne}"
OUTPUT_TAG="${OUTPUT_TAG:-$SCENARIO}"
KEEP="${KEEP_LAST:-50}"
REP=$((SLURM_ARRAY_TASK_ID / KEEP + 1))
ROW=$((SLURM_ARRAY_TASK_ID % KEEP + 2))
R=$(printf '%03d' "$REP")
PREP_BASE="${PREP_BASE:-$PROJECT/msprime_variable_ne_error/insertion_posterior/prepared}"
PREP="$PREP_BASE/replicate_$R"
MANIFEST="$PREP/draw_manifest.tsv"
DRAW_ID=$(awk -F '\t' -v row="$ROW" 'NR==row{print $1}' "$MANIFEST")
TREE=$(awk -F '\t' -v row="$ROW" 'NR==row{print $2}' "$MANIFEST")
case "$SCENARIO" in
  true_ne)
    NE="$PROJECT/msprime_variable_ne_error/simulations/replicate_$R/constant_ne_epochs.tsv"
    NE_SERIES=variable
    ;;
  estimated_ne)
    NE="$PROJECT/msprime_variable_ne_error/argtest/replicate_$R/ne/estimated_ne.tsv"
    NE_SERIES=estimated
    ;;
  *) echo "ERROR: unknown SCENARIO=$SCENARIO" >&2; exit 2 ;;
esac
test -n "$DRAW_ID" && test -s "$TREE" && test -s "$NE"
OUT="$PROJECT/msprime_variable_ne_error/insertion_posterior/$OUTPUT_TAG/draws/replicate_$R/draw_$DRAW_ID"
module load conda
conda activate dnaging
python "$PROJECT/insertion_real_data.py" run-draw --tree "$TREE" --draw-id "$DRAW_ID" \
  --calls "$PREP/calls" --ne "$NE" --ne-series "$NE_SERIES" \
  --epsilon "${EPSILON:-0.01}" --age-min "${AGE_MIN:-0}" \
  --age-max "${AGE_MAX:-15000}" --age-step "${AGE_STEP:-50}" --output "$OUT"
