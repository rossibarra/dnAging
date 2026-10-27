#!/bin/bash -l
# Merge draw uncertainty, refit epsilon, and run nested 5 Mb block bootstraps.
#SBATCH --job-name=rcInferMerge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=rcInferMerge-%j.out
#SBATCH --error=rcInferMerge-%j.err
set -euo pipefail

: "${RUN_ROOT:?set RUN_ROOT to one nine-demography run directory}"
: "${OUTROOT:?set OUTROOT to the Stage 2 directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
module load conda
conda activate dnaging
python "$PROJECT/validation/radiocarbon_inferred_stage2.py" \
  --parts "$OUTROOT/parts" --run-root "$RUN_ROOT" \
  --expected-epsilons 11 --expected-draws 10 \
  --rc-sd "${RC_SD:-200}" --n-folds "${N_FOLDS:-5}" \
  --block-bootstrap-replicates "${BLOCK_BOOTSTRAP_REPLICATES:-200}" \
  --individual-bootstrap-replicates "${INDIVIDUAL_BOOTSTRAP_REPLICATES:-2000}" \
  --output "$OUTROOT/result"
