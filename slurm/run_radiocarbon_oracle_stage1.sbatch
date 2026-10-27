#!/bin/bash -l
# One restartable array task evaluates one candidate epsilon on all three true ARGs.
#SBATCH --job-name=rcOracle1
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-20%10
#SBATCH --output=rcOracle1-%A_%a.out
#SBATCH --error=rcOracle1-%A_%a.err
set -euo pipefail

: "${RUN_ROOT:?set RUN_ROOT to one nine-demography run directory}"
: "${OUTROOT:?set OUTROOT to a new Stage 1 scenario directory}"
: "${INJECTED_EPSILON:?set INJECTED_EPSILON, for example 0.01}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
FIT_EPSILON="$(awk -v i="$SLURM_ARRAY_TASK_ID" 'BEGIN{printf "%.3f",i*0.005}')"
OUTPUT="$OUTROOT/parts/epsilon_${FIT_EPSILON}.npz"
if [[ -s "$OUTPUT" ]]; then
  echo "already complete: $OUTPUT"
  exit 0
fi

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
module load conda
conda activate dnaging
python "$PROJECT/validation/radiocarbon_oracle_stage1.py" likelihood \
  --run-root "$RUN_ROOT" \
  --injected-epsilon "$INJECTED_EPSILON" \
  --fit-epsilon "$FIT_EPSILON" \
  --injection-seed "${INJECTION_SEED:-2026101101}" \
  --output "$OUTPUT"
