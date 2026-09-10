#!/bin/bash
#SBATCH --job-name=TbetaMerge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
python -c 'import numpy' 2>/dev/null || {
    echo "ERROR: activate the dnaging conda environment before submission" >&2
    exit 90
}

parts=()
while IFS=$'\t' read -r tree_id _rest; do
    [[ "$tree_id" == "tree_number" || -z "$tree_id" ]] && continue
    parts+=("$DRAW_ROOT/draw_$tree_id")
done < "$DRAW_SELECTION"

python "$PROJECT/betabinom_real_data.py" --merge "${parts[@]}" --output "$OUTPUT"
