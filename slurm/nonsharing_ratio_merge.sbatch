#!/bin/bash -l
# Step 2: merge all draw NPZs in $OUTROOT/draws into $OUTROOT/ages_merged.tsv (single writer).
# Submit after the array finishes:  sbatch --dependency=afterok:<ARRAY_JOBID> slurm/nonsharing_ratio_merge.sbatch
#SBATCH --job-name=nsrMerge
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=nsrMerge-%j.out
#SBATCH --error=nsrMerge-%j.err
set -euo pipefail
: "${OUTROOT:?set OUTROOT to the run directory used by the draw array}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
test ! -e "$OUTROOT/ages_merged.tsv"
shopt -s nullglob
DRAWS=("$OUTROOT"/draws/draw_*.npz)
test "${#DRAWS[@]}" -gt 0 || { echo "no draw files in $OUTROOT/draws" >&2; exit 2; }
source /etc/profile.d/modules.sh
module load conda
conda activate dnaging
python "$PROJECT/pipeline/nonsharing_ratio_real_data.py" merge \
  --draw-results "${DRAWS[@]}" --bootstrap-replicates "${BOOTSTRAPS:-1000}" \
  --output "$OUTROOT/ages_merged.tsv"
