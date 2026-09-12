#!/bin/bash -l
#SBATCH --job-name=slim-edge-inf
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-99%100

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/slim_single_site_edge_validation}"
CHUNK_DIR="${CHUNK_DIR:-$ROOT/chunks}"
MUTATION_AGE_SOURCE="${MUTATION_AGE_SOURCE:-edge}"
CHUNK="${SLURM_ARRAY_TASK_ID:?array task required}"
START=$((CHUNK * 100 + 1))
STOP=$((START + 100))
module load conda
conda activate dnaging
python "$PROJECT/slim_edge_interval_validation.py" infer \
  --input-dir "$ROOT" \
  --frequency-table "$PROJECT/slim_single_site_diffusion_test_A/constant_ne_10000.freq_table.npz" \
  --output-dir "$CHUNK_DIR" --chunk-id "$CHUNK" --start "$START" --stop "$STOP" \
  --mutation-age-source "$MUTATION_AGE_SOURCE"
