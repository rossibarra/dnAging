#!/bin/bash -l
# Step 0 (once): stream the panel VCF and list sites fixed among panel haplotypes.
# One writer, one output file.  Submit from the repo root:
#   PANEL_VCF=/path/panel.vcf OUTROOT=/path/new_run_dir sbatch slurm/nonsharing_ratio_prepare_fixed.sbatch
#SBATCH --job-name=nsrFixed
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=nsrFixed-%j.out
#SBATCH --error=nsrFixed-%j.err
set -euo pipefail
: "${PANEL_VCF:?set PANEL_VCF to the (haploid) VCF of the ARG panel}"
: "${OUTROOT:?set OUTROOT to a new output directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
mkdir -p "$OUTROOT"
test ! -e "$OUTROOT/panel_fixed.npz"
source /etc/profile.d/modules.sh
module load conda
conda activate dnaging
python "$PROJECT/pipeline/nonsharing_ratio_real_data.py" prepare-fixed \
  --panel-vcf "$PANEL_VCF" --min-called "${MIN_CALLED:-21}" \
  --output "$OUTROOT/panel_fixed.npz"
