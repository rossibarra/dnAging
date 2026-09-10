#!/bin/bash
#SBATCH --job-name=Tbetabinom
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
python -c 'import numpy, mpmath, tszip, tskit' 2>/dev/null || {
    echo "ERROR: activate the dnaging conda environment before submission" >&2
    exit 90
}

mapfile -t DRAWS < <(awk 'NR > 1 && NF && $1 !~ /^#/' "$DRAW_SELECTION")
entry="${DRAWS[${SLURM_ARRAY_TASK_ID:-0}]}"
tree_id="$(awk '{print $1}' <<< "$entry")"
tree="${TREE_ROOT}/run.combined.${tree_id}.tsz"

python "$PROJECT/betabinom_real_data.py" \
    --tree "$tree" --draw-id "$tree_id" --calls "$CALLS" \
    --moment-table "$MOMENT_TABLE" --grid-table "$GRID_TABLE" \
    --epsilon "$EPSILON" --output "$OUTROOT/draw_${tree_id}"
