#!/bin/bash -l
#SBATCH --job-name=ewInfer
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=0-99%25

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/edge_width_ne_sweep}"
NE="${NE:?set NE}"
SIMS="${SIMS:?set SIMS}"
module load conda
conda activate dnaging
# One bin set across all three Ne: branch lengths scale with Ne, so bins tuned
# to one value pile the others into a single stratum and destroy the contrast.
python "$PROJECT/scripts/edge_width_scaling.py" infer \
    --input-root "$SIMS" \
    --frequency-table "$ROOT/tables/constant_ne_${NE}.freq_table.npz" \
    --output-dir "$ROOT/inference_ne${NE}" \
    --replicate "$((SLURM_ARRAY_TASK_ID + 1))" \
    --n-modern 26 --epsilon 0 --mutation-age-max-tau 3 \
    --width-edges 3000 10000 30000 100000 300000
