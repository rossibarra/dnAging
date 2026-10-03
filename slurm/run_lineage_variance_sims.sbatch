#!/bin/bash -l
#SBATCH --job-name=linvar
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --array=1-3

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="$PROJECT/msprime_lineage_variance"
# One age per replicate; ALL 20 ancient lineages sit at that same age in the SAME
# tree sequence, so the spread of T-hat across them isolates the lineage draw
# with tree, panel, demography and true age all held fixed.
case "$SLURM_ARRAY_TASK_ID" in
  1) AGE=2000 ;;
  2) AGE=5000 ;;
  3) AGE=8000 ;;
esac
module load conda
conda activate dnaging
python "$PROJECT/validation/simulate_variable_ne_error.py" \
    --output-root "$ROOT/simulations_age${AGE}" \
    --replicate 1 --base-seed "$((1760000000 + SLURM_ARRAY_TASK_ID))" \
    --ne-file "$ROOT/demography/constant_ne_50000.tsv" \
    --n-modern 33 --n-ancient 20 \
    --age-min "$AGE" --age-max "$AGE" \
    --length 10000000 --mutation-rate 1e-8 --recombination-rate 1e-8 \
    --epsilon 0
