#!/bin/bash -l
#SBATCH --job-name=slim-edge-sim
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=24:00:00
#SBATCH --array=1-10000%1000

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTDIR="${OUTDIR:-$PROJECT/slim_single_site_edge_validation}"
REPLICATE="${SLURM_ARRAY_TASK_ID:?array task required}"
BASE_SEED="${BASE_SEED:-202609090000}"
SEED=$((BASE_SEED + REPLICATE))
stem="replicate_$(printf '%06d' "$REPLICATE").seed_${SEED}"
mkdir -p "$OUTDIR/trees" "$OUTDIR/samples" "$OUTDIR/logs" "$OUTDIR/raw"
raw="$OUTDIR/raw/${stem}.trees"
tree="$OUTDIR/trees/${stem}.trees"
sample="$OUTDIR/samples/${stem}.samples.tsv"
[[ ! -e "$tree" && ! -e "$sample" ]] || { echo "output exists: $stem" >&2; exit 3; }
module load conda
conda activate dnaging
slim -s "$SEED" -d NE=10000 -d MUT_RATE=1e-6 -d END_GEN=40000 \
  -d N_PANEL=26 -d "TREE_FILE=\"$raw\"" -d "SAMPLE_FILE=\"$sample\"" \
  "$PROJECT/slim/single_site_edge_validation.slim" >"$OUTDIR/logs/${stem}.log" 2>&1
python "$PROJECT/slim_edge_interval_validation.py" compact \
  --input-tree "$raw" --sample-file "$sample" --output-tree "$tree"
rm -f "$raw"
