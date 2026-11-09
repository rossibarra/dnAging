#!/bin/bash -l
# Fixed-topology raw-SINGER versus POLEGON dating pilot, recent-growth chr1 draw 50.
#SBATCH --job-name=singerPolegon
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-1%2
#SBATCH --output=singerPolegon-%A_%a.out
#SBATCH --error=singerPolegon-%A_%a.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUN_ROOT="$PROJECT/batches/nine_demographies_20260916/runs/recent_growth_r1"
OUTROOT="${OUTROOT:-$PROJECT/batches/singer_polegon_pilot_20261102}"
SINGER_OUTPUT="$RUN_ROOT/singer/output"
POLEGON_TREE="$SINGER_OUTPUT/chr1/trees/chr1.50.tsz"
RAW_TREE="$OUTROOT/trees/chr1.50.raw_singer.tsz"
FAI="$OUTROOT/chr1.fai"
CALLS="$OUTROOT/calls"
RADIOCARBON="$PROJECT/batches/joint_recovery_matrix_20261029/oracle/baseline/recent_growth_r1/synthetic_radiocarbon.tsv"
LABELS=(raw_singer polegon)
LABEL="${LABELS[$SLURM_ARRAY_TASK_ID]}"
OUTPUT="$OUTROOT/$LABEL"

module load conda
conda activate dnaging
test -s "$RAW_TREE" && test -s "$FAI" && test -e "$CALLS/chr1.npz"
mkdir -p "$OUTROOT/checkpoints"
if [[ "$LABEL" == raw_singer ]]; then
  TREE="$RAW_TREE"
else
  TREE="$POLEGON_TREE"
fi

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
if [[ ! -s "$OUTPUT/draw/run.json" ]]; then
  python "$PROJECT/pipeline/insertion_real_data.py" run-draw \
    --tree "$TREE" --draw-id 50 --calls "$CALLS" --chrom-fai "$FAI" \
    --ne "$RUN_ROOT/demography.tsv" --ne-series recent_growth \
    --epsilon 0.01 --epsilon-grid 0:0.06:0.01 \
    --time-scale-grid 0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4 \
    --time-scale-chunk-size 3 --age-min 0 --age-max 10000 --age-step 100 \
    --block-size 5000000 --checkpoint "$OUTROOT/checkpoints/$LABEL.npz" \
    --checkpoint-every 5000 --output "$OUTPUT/draw"
fi
if [[ ! -s "$OUTPUT/merged/run.json" ]]; then
  python "$PROJECT/pipeline/insertion_real_data.py" merge \
    --draw-results "$OUTPUT/draw" --bootstrap-replicates 200 \
    --bootstrap-seed "$((2026110200+SLURM_ARRAY_TASK_ID))" --output "$OUTPUT/merged"
fi
if [[ ! -s "$OUTPUT/calibration/run.json" ]]; then
  python "$PROJECT/pipeline/joint_scale_epsilon_calibration.py" \
    --likelihood-dir "$OUTPUT/merged" --radiocarbon "$RADIOCARBON" --n-folds 5 \
    --bootstrap-replicates 2000 \
    --bootstrap-seed "$((2026110200+SLURM_ARRAY_TASK_ID))" \
    --output "$OUTPUT/calibration"
fi
