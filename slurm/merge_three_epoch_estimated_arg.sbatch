#!/bin/bash -l
#SBATCH --job-name=merge3chr
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=02:00:00
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"; CONDITION="${CONDITION:?}"
PREP="$PROJECT/msprime_three_epoch/insertion_three_chrom/prepared"; ROOT="$PROJECT/msprime_three_epoch/insertion_three_chrom/$CONDITION"
module load conda; conda activate dnaging
python "$PROJECT/pipeline/insertion_real_data.py" merge --manifest "$PREP/draw_manifest.tsv" --draw-root "$ROOT" --skip-d0 --output "$ROOT/merged"
