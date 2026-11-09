#!/bin/bash -l
# Build true-topology+POLEGON likelihoods, calibrate per chromosome, then combine ages.
#SBATCH --job-name=truePerChrom
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=18:00:00
#SBATCH --output=truePerChrom-%j.out
#SBATCH --error=truePerChrom-%j.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUN_ROOT="$PROJECT/batches/nine_demographies_20260916/runs/recent_growth_r1"
OUTROOT="${OUTROOT:-$PROJECT/batches/per_chromosome_calibration_20261103}"
TREE_DIR="$OUTROOT/true_topology_polegon/trees"
LIKELIHOOD="$OUTROOT/true_topology_polegon/likelihood"
CALIBRATION="$OUTROOT/true_topology_polegon/per_chromosome_calibration"
RADIOCARBON="$PROJECT/batches/joint_recovery_matrix_20261029/oracle/baseline/recent_growth_r1/synthetic_radiocarbon.tsv"

module load conda
conda activate dnaging
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
for chromosome in 1 2 3; do
  test -s "$TREE_DIR/chr$chromosome.polegon_true_topology.trees"
done
python "$PROJECT/validation/true_topology_polegon.py" likelihood \
  --run-root "$RUN_ROOT" --dated-trees "$TREE_DIR" \
  --time-scale-grid 0.7,0.8,0.9,1.0,1.1,1.2,1.25,1.3,1.4 \
  --epsilon-grid 0:0.06:0.01 --age-min 0 --age-max 10000 --age-step 100 \
  --output "$LIKELIHOOD"
python "$PROJECT/validation/per_chromosome_calibration.py" \
  --likelihood-dir "$LIKELIHOOD" --radiocarbon "$RADIOCARBON" \
  --n-folds 5 --output "$CALIBRATION"
