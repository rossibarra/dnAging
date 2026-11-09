#!/bin/bash -l
# One POLEGON-dated SINGER draw across all chromosomes for each demographic run.
#SBATCH --job-name=jointInfer
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=1-12:00:00
#SBATCH --array=0-8%9
#SBATCH --output=jointInfer-%A_%a.out
#SBATCH --error=jointInfer-%A_%a.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTROOT="${OUTROOT:-$PROJECT/batches/joint_recovery_matrix_20261029}"
ROW="$((SLURM_ARRAY_TASK_ID+2))"
RUN_NAME="$(awk -F '\t' -v n="$ROW" 'NR==n{print $2}' "$PROJECT/batches/nine_demographies_20260916/runs.tsv")"
MODEL="$(awk -F '\t' -v n="$ROW" 'NR==n{print $3}' "$PROJECT/batches/nine_demographies_20260916/runs.tsv")"
RUN_ROOT="$(awk -F '\t' -v n="$ROW" 'NR==n{print $5}' "$PROJECT/batches/nine_demographies_20260916/runs.tsv")"
OUTPUT="$OUTROOT/inferred/$RUN_NAME"
PREP="$RUN_ROOT/insertion_three_chrom/prepared"
DRAW_ID="$(awk -F '\t' 'NR==2{print $1}' "$PREP/draw_manifest.tsv")"
TREE="$(awk -F '\t' 'NR==2{print $2}' "$PREP/draw_manifest.tsv")"
DRAW="$OUTPUT/draw_$DRAW_ID"
CHECKPOINT="$OUTPUT/checkpoints/draw_$DRAW_ID.npz"
MERGED="$OUTPUT/merged"
CALIBRATION="$OUTPUT/calibration"
RADIOCARBON="$OUTROOT/oracle/baseline/$RUN_NAME/synthetic_radiocarbon.tsv"
test -n "$DRAW_ID" && test -f "$TREE" && test -f "$RADIOCARBON"
mkdir -p "$OUTPUT/checkpoints"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
module load conda
conda activate dnaging
if [[ ! -s "$DRAW/run.json" || ! -s "$DRAW/ll_by_chrom_scale_epsilon.npy" ]]; then
  python "$PROJECT/pipeline/insertion_real_data.py" run-draw \
    --tree "$TREE" --draw-id "$DRAW_ID" --calls "$PREP/calls" \
    --ne "$RUN_ROOT/demography.tsv" --ne-series "$MODEL" \
    --epsilon 0.01 --epsilon-grid "${EPSILON_GRID:-0:0.06:0.01}" \
    --time-scale-grid "${TIME_SCALE_GRID:-0.7,0.8,0.9,1.0,1.1,1.2,1.25,1.3,1.4}" \
    --time-scale-chunk-size "${TIME_SCALE_CHUNK_SIZE:-3}" \
    --age-min 0 --age-max 10000 --age-step 100 --block-size 5000000 \
    --checkpoint "$CHECKPOINT" --checkpoint-every 5000 --output "$DRAW"
fi
if [[ ! -s "$MERGED/run.json" || ! -s "$MERGED/ll_marginal_scale_epsilon.npy" ]]; then
  python "$PROJECT/pipeline/insertion_real_data.py" merge \
    --draw-results "$DRAW" --bootstrap-replicates 200 \
    --bootstrap-seed "$((2026102900+SLURM_ARRAY_TASK_ID))" --output "$MERGED"
fi
if [[ ! -s "$CALIBRATION/run.json" ]]; then
  python "$PROJECT/pipeline/joint_scale_epsilon_calibration.py" \
    --likelihood-dir "$MERGED" --radiocarbon "$RADIOCARBON" --n-folds 5 \
    --bootstrap-replicates 2000 --bootstrap-seed "$((2026102900+SLURM_ARRAY_TASK_ID))" \
    --output "$CALIBRATION"
fi
