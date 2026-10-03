#!/bin/bash -l
#SBATCH --job-name=insertDraw
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-9%10
set -euo pipefail
: "${MANIFEST:?set MANIFEST to a TSV with draw_id and tsz_path columns}"
: "${CALLS:?set CALLS to the prepared-call directory}"
: "${NE_TABLE:?set NE_TABLE to the piecewise-Ne TSV}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTROOT="${OUTROOT:-$PROJECT/insertion_arg_draws}"
ROW="$((SLURM_ARRAY_TASK_ID + 2))"
DRAW_ID="$(awk -F '\t' -v row="$ROW" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR==row{print $h["draw_id"]}' "$MANIFEST")"
TREE="$(awk -F '\t' -v row="$ROW" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR==row{print $h["tsz_path"]}' "$MANIFEST")"
test -n "$DRAW_ID" && test -n "$TREE"
module load conda
conda activate dnaging
python "$PROJECT/pipeline/insertion_real_data.py" run-draw --tree "$TREE" --draw-id "$DRAW_ID" \
  --calls "$CALLS" --ne "$NE_TABLE" --ne-series "${NE_SERIES:-combined}" \
  --epsilon "${EPSILON:-0.01}" --output "$OUTROOT/draw_$DRAW_ID"
