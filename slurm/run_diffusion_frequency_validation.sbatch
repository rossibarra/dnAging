#!/bin/bash -l
#SBATCH --job-name=diffFreqA
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=08:00:00
#SBATCH --array=0-99%25

set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
INPUT_DIR="${INPUT_DIR:-$PROJECT/slim_single_site_output}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PROJECT/slim_single_site_diffusion_test_A}"
FREQ_TABLE="${FREQ_TABLE:-$OUTPUT_ROOT/constant_ne_10000.freq_table.npz}"
N_CHUNKS="${N_CHUNKS:-100}"

module load conda
conda activate dnaging

python "$PROJECT/diffusion_frequency_validation.py" chunk \
    --input-dir "$INPUT_DIR" \
    --frequency-table "$FREQ_TABLE" \
    --output-dir "$OUTPUT_ROOT/chunks" \
    --chunk-id "${SLURM_ARRAY_TASK_ID:?This must be an array job}" \
    --n-chunks "$N_CHUNKS" \
    --n-panel "${N_PANEL:-26}" \
    --end-generation "${END_GENERATION:-40000}"
