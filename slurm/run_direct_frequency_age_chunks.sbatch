#!/bin/bash -l
#SBATCH --job-name=directAge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --array=0-99%25

set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
INPUT_DIR="${INPUT_DIR:-$PROJECT/slim_single_site_output}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PROJECT/slim_single_site_direct_age}"
N_CHUNKS="${N_CHUNKS:-100}"
END_GENERATION="${END_GENERATION:-40000}"

module load conda
conda activate dnaging

python "$PROJECT/direct_frequency_age_infer.py" chunk \
    --input-dir "$INPUT_DIR" \
    --output-dir "$OUTPUT_ROOT/chunks" \
    --chunk-id "${SLURM_ARRAY_TASK_ID:?This must be an array job}" \
    --n-chunks "$N_CHUNKS" \
    --end-generation "$END_GENERATION"
