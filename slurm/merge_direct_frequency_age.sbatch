#!/bin/bash -l
#SBATCH --job-name=mergeDirectAge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PROJECT/slim_single_site_direct_age}"
N_CHUNKS="${N_CHUNKS:-100}"

module load conda
conda activate dnaging

cache_root="${SLURM_TMPDIR:-/tmp}/direct-age-cache-${SLURM_JOB_ID:-$$}"
export MPLCONFIGDIR="$cache_root/matplotlib"
export XDG_CACHE_HOME="$cache_root/xdg"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

python "$PROJECT/direct_frequency_age_infer.py" merge \
    --chunk-dir "$OUTPUT_ROOT/chunks" \
    --output-dir "$OUTPUT_ROOT/results" \
    --n-chunks "$N_CHUNKS"
