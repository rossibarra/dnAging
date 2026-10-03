#!/bin/bash -l
#SBATCH --job-name=3epSim
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=1-3

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_three_epoch}"
module load conda
conda activate dnaging
# Three independent 50 Mb chromosomes of ONE dataset: they share the demography
# (--ne-file) AND the ancient individuals (--age-seed), while ancestry and
# mutation seeds stay replicate-specific so the genealogies are independent.
# Without the shared age seed each chromosome would carry a different set of
# ancient people, and concatenating their VCFs would be meaningless.
python "$PROJECT/simulate_variable_ne_error.py" \
    --output-root "$ROOT/simulations" \
    --replicate "$SLURM_ARRAY_TASK_ID" \
    --base-seed "${BASE_SEED:-1750000000}" \
    --age-seed "${AGE_SEED:-1750200000}" \
    --ne-file "$ROOT/demography/three_epoch.tsv" \
    --n-modern 100 --n-ancient 10 --age-min 0 --age-max 10000 \
    --length 50000000 --mutation-rate 1e-8 --recombination-rate 1e-8 \
    --epsilon 0
