#!/bin/bash -l
# Summarize fixed-topology raw-SINGER versus POLEGON dating pilot.
#SBATCH --job-name=sumSingerPolegon
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=sumSingerPolegon-%j.out
#SBATCH --error=sumSingerPolegon-%j.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTROOT="${OUTROOT:-$PROJECT/batches/singer_polegon_pilot_20261102}"
RUN_ROOT="$PROJECT/batches/nine_demographies_20260916/runs/recent_growth_r1"

module load conda
conda activate dnaging
python "$PROJECT/validation/singer_polegon_diagnostic.py" summarize \
  --singer-output "$RUN_ROOT/singer/output" --chromosome chr1 --draw 50 \
  --raw-tree "$OUTROOT/trees/chr1.50.raw_singer.tsz" \
  --polegon-tree "$RUN_ROOT/singer/output/chr1/trees/chr1.50.tsz" \
  --raw-calibration "$OUTROOT/raw_singer/calibration" \
  --polegon-calibration "$OUTROOT/polegon/calibration" \
  --output "$OUTROOT/summary.json"
