#!/bin/bash
#SBATCH --job-name=realIns1
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --time=32:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=realIns1-%j.out
#SBATCH --error=realIns1-%j.err

set -euo pipefail

source /etc/profile.d/modules.sh
module load conda
conda activate dnaging

cd "$SLURM_SUBMIT_DIR"

tree=/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/results/combined/tsz_files/run.combined.75.tsz
ne=radiocarbon_goodvcfs/insertion_pilot_20261012/demography/posterior_mean.from_zero.tsv
calls=radiocarbon_goodvcfs/insertion_pilot_20261012/calls
output=radiocarbon_goodvcfs/insertion_pilot_20261012/draws/draw_50

python pipeline/insertion_real_data.py run-draw \
    --tree "$tree" \
    --draw-id 50 \
    --calls "$calls" \
    --chrom-fai /quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/ref.fa.fai \
    --ne "$ne" \
    --ne-series posterior_mean \
    --ne-extrapolate-last \
    --epsilon 0.01 \
    --epsilon-grid 0:0.12:0.01 \
    --age-min 0 \
    --age-max 7000 \
    --age-step 25 \
    --block-size 5000000 \
    --checkpoint-every 25000 \
    --output "$output"
