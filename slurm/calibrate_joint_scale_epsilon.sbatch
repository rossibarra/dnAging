#!/bin/bash -l
# One task performs leakage-safe joint scale/epsilon cross-validation.
#SBATCH --job-name=jointCal
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=jointCal-%j.out
#SBATCH --error=jointCal-%j.err
set -euo pipefail

: "${LIKELIHOOD_DIR:?set LIKELIHOOD_DIR to the merged joint-grid result}"
: "${RADIOCARBON:?set RADIOCARBON to the long calibrated-density TSV}"
: "${OUTPUT:?set OUTPUT to a new result directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"

EXTRA=(--n-folds "${N_FOLDS:-0}" --bootstrap-replicates "${BOOTSTRAP_REPLICATES:-2000}" \
       --bootstrap-seed "${BOOTSTRAP_SEED:-20261028}")
[[ -n "${FOLDS:-}" ]] && EXTRA+=(--folds "$FOLDS")
[[ -n "${UNKNOWN_SAMPLES:-}" ]] && EXTRA+=(--unknown-samples "$UNKNOWN_SAMPLES")
[[ -n "${AGE_PRIOR:-}" ]] && EXTRA+=(--age-prior "$AGE_PRIOR")

module load conda
conda activate dnaging
python "$PROJECT/pipeline/joint_scale_epsilon_calibration.py" \
  --likelihood-dir "$LIKELIHOOD_DIR" \
  --radiocarbon "$RADIOCARBON" \
  --output "$OUTPUT" \
  "${EXTRA[@]}"
