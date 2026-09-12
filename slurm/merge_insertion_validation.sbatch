#!/bin/bash -l
#SBATCH --job-name=insertion-t9-merge
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="$PROJECT/slim_single_site_edge_validation"
module load conda
conda activate dnaging
python "$PROJECT/insertion_validation.py" merge \
  --chunk-dir "$ROOT/insertion_chunks" --output-dir "$ROOT/insertion_results"
