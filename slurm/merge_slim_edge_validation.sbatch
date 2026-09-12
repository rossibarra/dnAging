#!/bin/bash -l
#SBATCH --job-name=slim-edge-merge
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/slim_single_site_edge_validation}"
CHUNK_DIR="${CHUNK_DIR:-$ROOT/chunks}"
RESULT_DIR="${RESULT_DIR:-$ROOT/results}"
MUTATION_AGE_SOURCE="${MUTATION_AGE_SOURCE:-edge}"
module load conda
conda activate dnaging
python "$PROJECT/slim_edge_interval_validation.py" merge \
  --chunk-dir "$CHUNK_DIR" --output-dir "$RESULT_DIR" --n-chunks 100 \
  --mutation-age-source "$MUTATION_AGE_SOURCE"
