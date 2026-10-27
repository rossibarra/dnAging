#!/bin/bash -l
# This dependent job assembles the epsilon grid and performs five-fold validation.
#SBATCH --job-name=rcOracleMerge
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=rcOracleMerge-%j.out
#SBATCH --error=rcOracleMerge-%j.err
set -euo pipefail

: "${OUTROOT:?set OUTROOT to the Stage 1 scenario directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
module load conda
conda activate dnaging
python "$PROJECT/validation/radiocarbon_oracle_stage1.py" merge \
  --parts "$OUTROOT/parts" \
  --rc-sd "${RC_SD:-200}" \
  --n-folds "${N_FOLDS:-5}" \
  --bootstrap-replicates "${BOOTSTRAP_REPLICATES:-2000}" \
  --bootstrap-seed "${BOOTSTRAP_SEED:-20261011}" \
  --output "$OUTROOT/result"
