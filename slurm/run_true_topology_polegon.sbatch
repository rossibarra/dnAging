#!/bin/bash -l
# Date each existing recent-growth true ARG with POLEGON; topology is unchanged.
#SBATCH --job-name=truePolegon
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-2%3
#SBATCH --output=truePolegon-%A_%a.out
#SBATCH --error=truePolegon-%A_%a.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUN_ROOT="$PROJECT/batches/nine_demographies_20260916/runs/recent_growth_r1"
OUTROOT="${OUTROOT:-$PROJECT/batches/per_chromosome_calibration_20261103}"
INDEX="$((SLURM_ARRAY_TASK_ID+1))"
REP="$(printf '%03d' "$INDEX")"
NES=(7228.830449165209 7285.013149503907 7248.776124289828)
SEEDS=(2026110311 2026110312 2026110313)
TREE="$RUN_ROOT/simulations/replicate_$REP/known_modern_arg.trees"
OUTPUT="$OUTROOT/true_topology_polegon/trees/chr$INDEX.polegon_true_topology.trees"

module load conda
conda activate dnaging
mkdir -p "$(dirname "$OUTPUT")"
python "$PROJECT/validation/true_topology_polegon.py" date \
  --tree "$TREE" --output "$OUTPUT" --mutation-rate 1e-8 \
  --ne "${NES[$SLURM_ARRAY_TASK_ID]}" --seed "${SEEDS[$SLURM_ARRAY_TASK_ID]}"
