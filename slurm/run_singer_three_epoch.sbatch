#!/bin/bash -l
#SBATCH --job-name=3epSinger
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=10-00:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
W="$PROJECT/msprime_three_epoch/singer"
module load conda
conda activate singer-snakemake-env
# One run over all three chromosomes: 3 x 50 Mb, 100 haploid samples, constant
# recombination rate (no hapmap supplied), 100 MCMC samples with 50% burn-in.
snakemake \
    --snakefile /home/jri/src/singer-snakemake/workflow/Snakefile \
    --directory "$W" \
    --configfile "$W/config.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-64}" \
    --rerun-incomplete
