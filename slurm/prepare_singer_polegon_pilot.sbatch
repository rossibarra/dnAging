#!/bin/bash -l
# Reconstruct one raw-SINGER tree and stage immutable chr1 pilot inputs.
#SBATCH --job-name=prepSingerPolegon
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=prepSingerPolegon-%j.out
#SBATCH --error=prepSingerPolegon-%j.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUN_ROOT="$PROJECT/batches/nine_demographies_20260916/runs/recent_growth_r1"
OUTROOT="${OUTROOT:-$PROJECT/batches/singer_polegon_pilot_20261102}"
SINGER_OUTPUT="$RUN_ROOT/singer/output"
RAW_TREE="$OUTROOT/trees/chr1.50.raw_singer.tsz"
FAI="$OUTROOT/chr1.fai"
CALLS="$OUTROOT/calls"

module load conda
conda activate dnaging
mkdir -p "$OUTROOT/trees" "$CALLS" "$OUTROOT/checkpoints"
if [[ ! -e "$CALLS/chr1.npz" ]]; then
  ln -s "$RUN_ROOT/insertion_three_chrom/prepared/calls/chr1.npz" "$CALLS/chr1.npz"
fi
if [[ ! -s "$RAW_TREE" ]]; then
  python "$PROJECT/validation/singer_polegon_diagnostic.py" reconstruct \
    --singer-output "$SINGER_OUTPUT" --chromosome chr1 --draw 50 --output "$RAW_TREE"
fi
if [[ ! -s "$FAI" ]]; then
  python "$PROJECT/validation/singer_polegon_diagnostic.py" write-fai \
    --tree "$RAW_TREE" --chromosome chr1 --output "$FAI"
fi
