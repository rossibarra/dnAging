#!/bin/bash -l
# Array index maps to 11 epsilon candidates x 10 retained inferred-ARG draws.
#SBATCH --job-name=rcInfer2
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=16:00:00
#SBATCH --array=0-109%20
#SBATCH --output=rcInfer2-%A_%a.out
#SBATCH --error=rcInfer2-%A_%a.err
set -euo pipefail

: "${RUN_ROOT:?set RUN_ROOT to one nine-demography run directory}"
: "${OUTROOT:?set OUTROOT to a new Stage 2 directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
DRAW_INDEX="$((SLURM_ARRAY_TASK_ID % 10))"
EPSILON_INDEX="$((SLURM_ARRAY_TASK_ID / 10))"
FIT_EPSILON="$(awk -v i="$EPSILON_INDEX" 'BEGIN{printf "%.2f",i*0.01}')"
PREP="$RUN_ROOT/insertion_three_chrom/prepared"
ROW="$((DRAW_INDEX+2))"
DRAW_ID="$(awk -F '\t' -v n="$ROW" 'NR==n{print $1}' "$PREP/draw_manifest.tsv")"
TREE="$(awk -F '\t' -v n="$ROW" 'NR==n{print $2}' "$PREP/draw_manifest.tsv")"
test -n "$DRAW_ID" && test -f "$TREE"
FINAL="$OUTROOT/parts/epsilon_$FIT_EPSILON/draw_$DRAW_ID"
if [[ -f "$FINAL/run.json" && -f "$FINAL/ll_by_block.npy" ]]; then
  echo "already complete: $FINAL"
  exit 0
fi
mkdir -p "$OUTROOT/parts/epsilon_$FIT_EPSILON"
TEMP="$OUTROOT/parts/epsilon_$FIT_EPSILON/.draw_${DRAW_ID}.tmp.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
module load conda
conda activate dnaging
python "$PROJECT/pipeline/insertion_real_data.py" run-draw \
  --tree "$TREE" --draw-id "$DRAW_ID" --calls "$PREP/calls" \
  --ne "$RUN_ROOT/demography.tsv" --ne-series constant \
  --epsilon "$FIT_EPSILON" --age-min 0 --age-max 10000 --age-step 50 \
  --block-size 5000000 --output "$TEMP"
mv "$TEMP" "$FINAL"
