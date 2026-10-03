#!/bin/bash -l
#SBATCH --job-name=nine_eps01
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=nine-eps01-%A_%a.out
set -euo pipefail
PROJECT=/quobyte/jrigrp/beil/logan_collab/dnAging
BATCH_ROOT="$PROJECT/batches/nine_demographies_20260916"
export PYTHONPATH="$BATCH_ROOT/code:$PROJECT:${PYTHONPATH:-}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK" MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
module load conda
conda activate dnaging
RUN_LINE=$(awk -F '\t' -v n="$((RUN_INDEX+2))" 'NR==n' "$BATCH_ROOT/runs.tsv")
IFS=$'\t' read -r INDEX RUN_NAME MODEL REPLICATE RUN_ROOT BASE_SEED SINGER_SEED <<< "$RUN_LINE"
test "$INDEX" = "$RUN_INDEX"
export RUN_ROOT
cd "$RUN_ROOT"
PREP="$RUN_ROOT/insertion_three_chrom/prepared"
RESULT="$RUN_ROOT/insertion_three_chrom/estimated_arg_estimated_ne_eps_0.01"
if [[ "${STAGE:-draw}" == merge ]]; then
    python "$BATCH_ROOT/code/insertion_real_data.py" merge --manifest "$PREP/draw_manifest.tsv" \
        --draw-root "$RESULT" --skip-d0 --output "$RESULT/merged"
else
    ROW=$((SLURM_ARRAY_TASK_ID+2))
    DRAW_ID=$(awk -F '\t' -v n="$ROW" 'NR==n{print $1}' "$PREP/draw_manifest.tsv")
    TREE=$(awk -F '\t' -v n="$ROW" 'NR==n{print $2}' "$PREP/draw_manifest.tsv")
    test -n "$DRAW_ID" && test -f "$TREE"
    python "$BATCH_ROOT/code/insertion_real_data.py" run-draw --tree "$TREE" --draw-id "$DRAW_ID" \
        --calls "$PREP/calls" --ne "$RUN_ROOT/argtest/ne/estimated_ne.tsv" --ne-series estimated \
        --epsilon 0.01 --age-min 0 --age-max 10000 --age-step 50 \
        --output "$RESULT/draw_$DRAW_ID"
fi
