#!/bin/bash
#SBATCH --job-name=TbetaCalls
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
python -c 'import numpy' 2>/dev/null || {
    echo "ERROR: activate the dnaging conda environment before submission" >&2
    exit 90
}
python "$PROJECT/prepare_betabinom_calls.py" \
    --epsilon-data-root "$EPSILON_DATA_ROOT" --samples "$SAMPLES_FILE" \
    --ancient-vcf "$ANCIENT_VCF" --output "$OUTPUT"
