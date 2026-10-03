#!/bin/bash -l
#SBATCH --job-name=ins3chr
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=36:00:00
#SBATCH --array=0-49%20
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"; CONDITION="${CONDITION:?}"; NE_TABLE="${NE_TABLE:?}"; NE_SERIES="${NE_SERIES:?}"
PREP="$PROJECT/msprime_three_epoch/insertion_three_chrom/prepared"; ROW="$((SLURM_ARRAY_TASK_ID+2))"
DRAW_ID="$(awk -F '\t' -v n="$ROW" 'NR==n{print $1}' "$PREP/draw_manifest.tsv")"
TREE="$(awk -F '\t' -v n="$ROW" 'NR==n{print $2}' "$PREP/draw_manifest.tsv")"
module load conda; conda activate dnaging
python "$PROJECT/insertion_real_data.py" run-draw --tree "$TREE" --draw-id "$DRAW_ID" --calls "$PREP/calls" \
 --ne "$NE_TABLE" --ne-series "$NE_SERIES" --epsilon "${EPSILON:-0}" --age-min 0 --age-max 10000 --age-step 50 \
 --output "$PROJECT/msprime_three_epoch/insertion_three_chrom/$CONDITION/draw_$DRAW_ID"
