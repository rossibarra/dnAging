#!/bin/bash -l
#SBATCH --job-name=slim-edge-uniform
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
module load conda
conda activate dnaging
python "$PROJECT/edge_uniformity_validation.py" \
  --input-dir "$PROJECT/slim_single_site_edge_validation" \
  --chunk-dir "$PROJECT/slim_single_site_edge_validation/chunks_random_focal" \
  --output-dir "$PROJECT/slim_single_site_edge_validation/edge_uniformity"
