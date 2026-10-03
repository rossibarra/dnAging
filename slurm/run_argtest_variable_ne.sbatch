#!/bin/bash -l
#SBATCH --job-name=vneArgtest
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-10%3

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ARGTEST_REPO="${ARGTEST_REPO:-/home/jri/src/argtest}"
REP=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
WORK="$PROJECT/msprime_variable_ne_error/argtest/replicate_${REP}"

module load conda
conda activate argtest
# argtest's Snakefile invokes `python scripts/...` relatively, so the working
# directory must expose a scripts/ dir; each replicate symlinks the repo's.
# That lets --directory give every replicate its own .snakemake lock, which a
# shared working directory would not.
snakemake \
    --snakefile "$ARGTEST_REPO/Snakefile" \
    --directory "$WORK" \
    --configfile "$WORK/config.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-16}" \
    --rerun-incomplete
