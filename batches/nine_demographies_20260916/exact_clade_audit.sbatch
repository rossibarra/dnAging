#!/bin/bash -l
#SBATCH --job-name=nine_exact_clade
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=nine-exact-clade-%A_%a.out

set -euo pipefail

PROJECT=/quobyte/jrigrp/beil/logan_collab/dnAging
BATCH_ROOT="$PROJECT/batches/nine_demographies_20260916"
export PYTHONPATH="$BATCH_ROOT/code:$PROJECT:${PYTHONPATH:-}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"

module load conda
conda activate dnaging

RUN_LINE=$(awk -F '\t' -v n="$((RUN_INDEX+2))" 'NR==n' "$BATCH_ROOT/runs.tsv")
IFS=$'\t' read -r INDEX RUN_NAME MODEL REPLICATE RUN_ROOT BASE_SEED SINGER_SEED <<< "$RUN_LINE"
test "$INDEX" = "$RUN_INDEX"
cd "$RUN_ROOT"

PREP="$RUN_ROOT/insertion_three_chrom/prepared"
AUDIT="$RUN_ROOT/insertion_three_chrom/exact_clade_audit"
STRICT="$RUN_ROOT/insertion_three_chrom/prepared_exact_clade_all_50"

if [[ "${STAGE:-audit}" == merge ]]; then
    python "$BATCH_ROOT/code/audit_exact_clade_matches.py" merge \
        --manifest "$PREP/draw_manifest.tsv" \
        --calls "$PREP/calls" \
        --audit-root "$AUDIT" \
        --output "$STRICT"
else
    ROW=$((SLURM_ARRAY_TASK_ID+2))
    DRAW_ID=$(awk -F '\t' -v n="$ROW" 'NR==n{print $1}' "$PREP/draw_manifest.tsv")
    TREE=$(awk -F '\t' -v n="$ROW" 'NR==n{print $2}' "$PREP/draw_manifest.tsv")
    test -n "$DRAW_ID" && test -f "$TREE"
    python "$BATCH_ROOT/code/audit_exact_clade_matches.py" audit-draw \
        --tree "$TREE" \
        --draw-id "$DRAW_ID" \
        --calls "$PREP/calls" \
        --simulations "$RUN_ROOT/simulations" \
        --singer-input "$RUN_ROOT/singer/input" \
        --n-modern 100 \
        --output "$AUDIT/draw_$DRAW_ID"
fi
