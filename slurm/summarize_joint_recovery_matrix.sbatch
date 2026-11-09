#!/bin/bash -l
#SBATCH --job-name=jointRecSum
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=jointRecSum-%j.out
#SBATCH --error=jointRecSum-%j.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTROOT="${OUTROOT:-$PROJECT/batches/joint_recovery_matrix_20261029}"
module load conda
conda activate dnaging
python "$PROJECT/validation/joint_scale_epsilon_recovery.py" summarize \
  --results "$OUTROOT/oracle" --expected-results 45 \
  --fail-on-unrecovered --output "$OUTROOT/oracle_summary"
