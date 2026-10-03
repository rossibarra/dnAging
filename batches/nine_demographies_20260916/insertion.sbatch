#!/bin/bash -l
#SBATCH --job-name=nine_insertion
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=1-12:00:00
#SBATCH --array=0-49

set -euo pipefail
export PROJECT=/quobyte/jrigrp/beil/logan_collab/dnAging
export BATCH_ROOT=/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nine_demographies_20260916
export PYTHONPATH="$BATCH_ROOT/code:$PROJECT:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load conda
conda activate dnaging
RUN_INDEX="${RUN_INDEX:?}"
RUN_LINE=$(awk -F '\t' -v n="$((RUN_INDEX+2))" 'NR==n' "$BATCH_ROOT/runs.tsv")
IFS=$'\t' read -r INDEX RUN_NAME MODEL REPLICATE RUN_ROOT BASE_SEED SINGER_SEED <<< "$RUN_LINE"
export RUN_ROOT
test "$INDEX" = "$RUN_INDEX"
cd "$RUN_ROOT"

PREP="$RUN_ROOT/insertion_three_chrom/prepared"
ROW="$((SLURM_ARRAY_TASK_ID+2))"
DRAW_ID=$(awk -F '\t' -v n="$ROW" 'NR==n{print $1}' "$PREP/draw_manifest.tsv")
TREE=$(awk -F '\t' -v n="$ROW" 'NR==n{print $2}' "$PREP/draw_manifest.tsv")
test -n "$DRAW_ID" && test -f "$TREE"
python "$BATCH_ROOT/code/insertion_real_data.py" run-draw --tree "$TREE" --draw-id "$DRAW_ID" --calls "$PREP/calls" \
 --ne "$RUN_ROOT/argtest/ne/estimated_ne.tsv" --ne-series estimated --epsilon 0 \
 --age-min 0 --age-max 10000.0 --age-step 50.0 \
 --output "$RUN_ROOT/insertion_three_chrom/estimated_arg_estimated_ne/draw_$DRAW_ID"
