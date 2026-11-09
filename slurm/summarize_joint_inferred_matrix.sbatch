#!/bin/bash -l
#SBATCH --job-name=jointInfSum
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=jointInfSum-%j.out
#SBATCH --error=jointInfSum-%j.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTROOT="${OUTROOT:-$PROJECT/batches/joint_recovery_matrix_20261029}"
module load conda
conda activate dnaging
python "$PROJECT/validation/joint_scale_epsilon_recovery.py" summarize-inferred \
  --results "$OUTROOT/inferred" --expected-results 9 \
  --output "$OUTROOT/inferred_summary"
