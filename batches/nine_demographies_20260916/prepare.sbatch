#!/bin/bash -l
#SBATCH --job-name=nine_prepare
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

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

python "$BATCH_ROOT/code/merge_three_epoch_draw_masks.py"
