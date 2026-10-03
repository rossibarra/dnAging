#!/bin/bash -l
#SBATCH --job-name=vneInsMerge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-10%10
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
SCENARIO="${SCENARIO:?set SCENARIO to true_ne or estimated_ne}"
OUTPUT_TAG="${OUTPUT_TAG:-$SCENARIO}"
MERGED_TAG="${MERGED_TAG:-$OUTPUT_TAG}"
R=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
BASE="$PROJECT/msprime_variable_ne_error/insertion_posterior"
PREP="$BASE/prepared/replicate_$R"
module load conda
conda activate dnaging
EXTRA=()
if [[ "${SKIP_D0:-0}" == 1 ]]; then EXTRA+=(--skip-d0); fi
python "$PROJECT/insertion_real_data.py" merge --manifest "$PREP/draw_manifest.tsv" \
  --draw-root "$BASE/$OUTPUT_TAG/draws/replicate_$R" "${EXTRA[@]}" \
  --output "$BASE/$MERGED_TAG/merged/replicate_$R"
