#!/bin/bash -l
# Step 1: one array task per ARG draw -> $OUTROOT/draws/draw_<id>.npz (one writer per file).
# Draws are the *.tsz / *.trees files directly inside DRAW_DIR, sorted by name; the draw id is
# the trailing number of the file name (run.combined.100.tsz -> 100).
# Size the array to the number of draws, e.g.  sbatch --array=0-9 ...  (default 0-0 = first draw).
# See docs/nonsharing_ratio/RUNNING.md for the required variables.
#SBATCH --job-name=nsrDraw
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=06:00:00
#SBATCH --array=0-0
#SBATCH --output=nsrDraw-%A_%a.out
#SBATCH --error=nsrDraw-%A_%a.err
set -euo pipefail
: "${DRAW_DIR:?set DRAW_DIR to the directory of ARG draws (*.tsz or *.trees)}"
: "${CONTROLS_VCF:?set CONTROLS_VCF to the held-out modern controls VCF}"
: "${ANCIENT_VCF:?set ANCIENT_VCF to the ancient pseudo-haploid VCF}"
: "${CHROM_FAI:?set CHROM_FAI to the FAI giving the concatenation offsets}"
: "${A0:?set A0 (generations): keep segregating sites dated older than this}"
: "${OUTROOT:?set OUTROOT to a new output directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
mapfile -t TREES < <(find "$DRAW_DIR" -maxdepth 1 \( -name '*.tsz' -o -name '*.trees' \) | sort)
TREE="${TREES[$SLURM_ARRAY_TASK_ID]:-}"
test -n "$TREE" || { echo "no draw for array index $SLURM_ARRAY_TASK_ID in $DRAW_DIR" >&2; exit 2; }
DRAW_ID="$(basename "$TREE" | sed -E 's/.*[^0-9]([0-9]+)\.(tsz|trees)$/\1/')"
mkdir -p "$OUTROOT/draws"
EXTRA=()
[[ -n "${NE_TSV:-}" ]] && EXTRA+=(--ne "$NE_TSV" --ne-series "${NE_SERIES:-posterior_mean}")
[[ -n "${DEMES_YAML:-}" ]] && EXTRA+=(--demography "$DEMES_YAML")
[[ -n "${PANEL_FIXED:-}" ]] && EXTRA+=(--panel-fixed "$PANEL_FIXED" --fixed-polarity "${FIXED_POLARITY:-table}")
[[ -n "${POLARITY_TABLE:-}" ]] && EXTRA+=(--polarity-table "$POLARITY_TABLE")
[[ -n "${CHROMS:-}" ]] && EXTRA+=(--chroms $CHROMS)
[[ -n "${DEIDENTIFY:-}" ]] && EXTRA+=(--deidentify)
[[ -n "${PANEL_NAMES:-}" ]] && EXTRA+=(--panel-names "$PANEL_NAMES")
source /etc/profile.d/modules.sh
module load conda
conda activate dnaging
python "$PROJECT/pipeline/nonsharing_ratio_real_data.py" run-draw \
  --tree "$TREE" --draw-id "$DRAW_ID" --controls-vcf "$CONTROLS_VCF" \
  --controls-mode "${CONTROLS_MODE:-phased}" --ancient-vcf "$ANCIENT_VCF" \
  --chrom-fai "$CHROM_FAI" --a0 "$A0" --position-shift "${POSITION_SHIFT:-0}" \
  ${EXTRA[@]+"${EXTRA[@]}"} --output "$OUTROOT/draws/draw_$DRAW_ID.npz"
