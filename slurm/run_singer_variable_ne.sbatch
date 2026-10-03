#!/bin/bash -l
#SBATCH --job-name=vneSinger
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --array=1-10%3

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
SINGER_REPO="${SINGER_REPO:-/home/jri/src/singer-snakemake}"
REP=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
WORK="$PROJECT/msprime_variable_ne_error/singer/replicate_${REP}"

module load conda
conda activate singer-snakemake-env
# --directory keeps .snakemake/ inside the replicate so the ten runs cannot
# collide on a shared lock; the workflow finds its bundled SINGER and POLEGON
# binaries through workflow.basedir, which is independent of the working dir.
snakemake \
    --snakefile "$SINGER_REPO/workflow/Snakefile" \
    --directory "$WORK" \
    --configfile "$WORK/config.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-32}" \
    --rerun-incomplete
