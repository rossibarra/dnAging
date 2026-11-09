#!/bin/bash -l
# Split true ARGs at production boundaries, date each chunk, and merge both controls.
#SBATCH --job-name=chunkTruePolegon
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-2%3
#SBATCH --output=chunkTruePolegon-%A_%a.out
#SBATCH --error=chunkTruePolegon-%A_%a.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUN_ROOT="$PROJECT/batches/nine_demographies_20260916/runs/recent_growth_r1"
OUTROOT="${OUTROOT:-$PROJECT/batches/per_chromosome_calibration_20261103}"
INDEX="$((SLURM_ARRAY_TASK_ID+1))"
REP="$(printf '%03d' "$INDEX")"
OUTPUT="$OUTROOT/true_topology_polegon_chunk_matched"

module load conda
conda activate dnaging
python "$PROJECT/validation/chunked_true_topology_polegon.py" \
  --tree "$RUN_ROOT/simulations/replicate_$REP/known_modern_arg.trees" \
  --parameter-dir "$RUN_ROOT/singer/output/chr$INDEX/chunks" \
  --chunk-output "$OUTPUT/chunks/chr$INDEX" \
  --true-output "$OUTPUT/trees/chr$INDEX.chunked_true_times.trees" \
  --polegon-output "$OUTPUT/trees/chr$INDEX.chunked_polegon_true_topology.trees" \
  --metadata "$OUTPUT/trees/chr$INDEX.chunked_comparison.json" \
  --draw 50 --expected-chunks 50 --mutation-rate 1e-8
