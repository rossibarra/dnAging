#!/bin/bash -l
#SBATCH --job-name=mergeDiffA
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PROJECT/slim_single_site_diffusion_test_A}"
N_CHUNKS="${N_CHUNKS:-100}"

module load conda
conda activate dnaging

python "$PROJECT/diffusion_frequency_validation.py" merge \
    --chunk-dir "$OUTPUT_ROOT/chunks" \
    --output-dir "$OUTPUT_ROOT/results" \
    --n-chunks "$N_CHUNKS"
