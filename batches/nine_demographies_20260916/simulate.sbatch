#!/bin/bash -l
#SBATCH --job-name=nine_simulate
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=1-3

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

python "$BATCH_ROOT/code/simulate_grid.py" --output-root "$RUN_ROOT/simulations" \
 --replicate "$SLURM_ARRAY_TASK_ID" --base-seed "$BASE_SEED" --ne-file "$RUN_ROOT/demography.tsv" \
 --n-modern 100 --n-ancient 21 --length 50000000 --mutation-rate 1e-8 --recombination-rate 1e-8 \
 --epsilon 0 --years-per-generation 1.0
python "$BATCH_ROOT/code/prepare_singer_input.py" --simulations "$RUN_ROOT/simulations" \
 --replicate "$SLURM_ARRAY_TASK_ID" --chrom "chr$SLURM_ARRAY_TASK_ID" --input-dir "$RUN_ROOT/singer/input"
