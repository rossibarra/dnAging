#!/bin/bash -l
# Five controlled perturbations across all nine demographic simulation runs.
#SBATCH --job-name=jointRecover
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-12:00:00
#SBATCH --array=0-44%12
#SBATCH --output=jointRecover-%A_%a.out
#SBATCH --error=jointRecover-%A_%a.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTROOT="${OUTROOT:-$PROJECT/batches/joint_recovery_matrix_20261029}"
RUN_INDEX="$((SLURM_ARRAY_TASK_ID % 9))"
SCENARIO_INDEX="$((SLURM_ARRAY_TASK_ID / 9))"
ROW="$((RUN_INDEX+2))"
RUN_NAME="$(awk -F '\t' -v n="$ROW" 'NR==n{print $2}' "$PROJECT/batches/nine_demographies_20260916/runs.tsv")"
RUN_ROOT="$(awk -F '\t' -v n="$ROW" 'NR==n{print $5}' "$PROJECT/batches/nine_demographies_20260916/runs.tsv")"
SCENARIOS=(baseline compressed expanded discordance joint)
INPUT_SCALES=(1.0 0.8 1.25 1.0 0.8)
INJECTED_EPSILONS=(0.0 0.0 0.0 0.02 0.02)
SCENARIO="${SCENARIOS[$SCENARIO_INDEX]}"
OUTPUT="$OUTROOT/oracle/$SCENARIO/$RUN_NAME"
if [[ -s "$OUTPUT/recovery_run.json" && -s "$OUTPUT/calibration/run.json" ]]; then
  echo "already complete: $OUTPUT"
  exit 0
fi

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
module load conda
conda activate dnaging
python "$PROJECT/validation/joint_scale_epsilon_recovery.py" run \
  --run-root "$RUN_ROOT" --scenario "$SCENARIO" \
  --input-time-scale "${INPUT_SCALES[$SCENARIO_INDEX]}" \
  --injected-epsilon "${INJECTED_EPSILONS[$SCENARIO_INDEX]}" \
  --injection-seed "$((2026102900+RUN_INDEX))" \
  --time-scale-grid "${TIME_SCALE_GRID:-0.7,0.8,0.9,1.0,1.1,1.2,1.25,1.3,1.4}" \
  --epsilon-grid "${EPSILON_GRID:-0:0.06:0.01}" \
  --output "$OUTPUT"
