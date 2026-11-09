#!/bin/bash -l
# Reuse the completed inferred-ARG likelihood tensor for per-chromosome calibration.
#SBATCH --job-name=inferPerChrom
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=inferPerChrom-%j.out
#SBATCH --error=inferPerChrom-%j.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTROOT="${OUTROOT:-$PROJECT/batches/per_chromosome_calibration_20261103}"
LIKELIHOOD="$PROJECT/batches/joint_recovery_matrix_20261029/inferred/recent_growth_r1/draw_50"
RADIOCARBON="$PROJECT/batches/joint_recovery_matrix_20261029/oracle/baseline/recent_growth_r1/synthetic_radiocarbon.tsv"
OUTPUT="$OUTROOT/inferred_topology_polegon/per_chromosome_calibration"

module load conda
conda activate dnaging
python "$PROJECT/validation/per_chromosome_calibration.py" \
  --likelihood-dir "$LIKELIHOOD" --radiocarbon "$RADIOCARBON" \
  --n-folds 5 --output "$OUTPUT"
