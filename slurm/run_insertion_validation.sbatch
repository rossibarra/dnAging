#!/bin/bash -l
#SBATCH --job-name=insertion-t9
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=0-99%100

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="$PROJECT/slim_single_site_edge_validation"
module load conda
conda activate dnaging
python "$PROJECT/insertion_validation.py" run-chunk \
  --input-dir "$ROOT" --diagnostic-dir "$ROOT/chunks_random_focal" \
  --output-dir "$ROOT/insertion_chunks" --chunk-id "$SLURM_ARRAY_TASK_ID"
