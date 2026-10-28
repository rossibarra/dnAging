#!/bin/bash -l
# Exact Phase-B pilot: one inferred ARG draw and one chromosome on a coarse grid.
#SBATCH --job-name=jointPilot
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=jointPilot-%j.out
#SBATCH --error=jointPilot-%j.err
set -euo pipefail

: "${RUN_ROOT:?set RUN_ROOT to one nine-demography run directory}"
: "${OUTROOT:?set OUTROOT to a new joint-pilot directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
PREP="$RUN_ROOT/insertion_three_chrom/prepared"
DRAW_INDEX="${DRAW_INDEX:-0}"
ROW="$((DRAW_INDEX+2))"
DRAW_ID="$(awk -F '\t' -v n="$ROW" 'NR==n{print $1}' "$PREP/draw_manifest.tsv")"
TREE="$(awk -F '\t' -v n="$ROW" 'NR==n{print $2}' "$PREP/draw_manifest.tsv")"
test -n "$DRAW_ID" && test -f "$TREE"

FINAL="$OUTROOT/draw_$DRAW_ID"
CHECKPOINT="$OUTROOT/checkpoints/draw_$DRAW_ID.npz"
if [[ -f "$FINAL/run.json" && -f "$FINAL/ll_by_chrom_scale_epsilon.npy" ]]; then
  echo "already complete: $FINAL"
  exit 0
fi
mkdir -p "$OUTROOT/checkpoints" "$TMPDIR/calls"
# The pilot deliberately exposes only chr1 to the concatenated-ARG driver.
cp "$PREP/calls/chr1.npz" "$TMPDIR/calls/chr1.npz"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
module load conda
conda activate dnaging
python "$PROJECT/pipeline/insertion_real_data.py" run-draw \
  --tree "$TREE" --draw-id "$DRAW_ID" --calls "$TMPDIR/calls" \
  --ne "$RUN_ROOT/demography.tsv" --ne-series constant \
  --epsilon 0.01 --epsilon-grid "${EPSILON_GRID:-0:0.10:0.02}" \
  --time-scale-grid "${TIME_SCALE_GRID:-0.8:1.6:0.1}" \
  --time-scale-chunk-size "${TIME_SCALE_CHUNK_SIZE:-3}" \
  --age-min 0 --age-max 10000 --age-step "${AGE_STEP:-100}" \
  --block-size 5000000 --checkpoint "$CHECKPOINT" \
  --checkpoint-every "${CHECKPOINT_EVERY:-5000}" --output "$FINAL"
