#!/bin/bash -l
# Compute and cross-validate one chunk-matched true-topology timing arm.
#SBATCH --job-name=chunkTrueCal
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=18:00:00
#SBATCH --output=chunkTrueCal-%j.out
#SBATCH --error=chunkTrueCal-%j.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUN_ROOT="$PROJECT/batches/nine_demographies_20260916/runs/recent_growth_r1"
OUTROOT="${OUTROOT:-$PROJECT/batches/per_chromosome_calibration_20261103}"
BASE="$OUTROOT/true_topology_polegon_chunk_matched"
RADIOCARBON="$PROJECT/batches/joint_recovery_matrix_20261029/oracle/baseline/recent_growth_r1/synthetic_radiocarbon.tsv"
ARM="${ARM:?Set ARM to true_times or polegon}"
case "$ARM" in
  true_times) PATTERN='chr{chromosome}.chunked_true_times.trees' ;;
  polegon) PATTERN='chr{chromosome}.chunked_polegon_true_topology.trees' ;;
  *) echo "invalid ARM: $ARM" >&2; exit 2 ;;
esac
LIKELIHOOD="$BASE/$ARM/likelihood"
CALIBRATION="$BASE/$ARM/per_chromosome_calibration"

module load conda
conda activate dnaging
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
python "$PROJECT/validation/true_topology_polegon.py" likelihood \
  --run-root "$RUN_ROOT" --dated-trees "$BASE/trees" --tree-pattern "$PATTERN" \
  --time-scale-grid 0.7,0.8,0.9,1.0,1.1,1.2,1.25,1.3,1.4 \
  --epsilon-grid 0:0.06:0.01 --age-min 0 --age-max 10000 --age-step 100 \
  --output "$LIKELIHOOD"
python "$PROJECT/validation/per_chromosome_calibration.py" \
  --likelihood-dir "$LIKELIHOOD" --radiocarbon "$RADIOCARBON" \
  --n-folds 5 --output "$CALIBRATION"
