#!/bin/bash -l
#SBATCH --job-name=vneInsSum
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
SCENARIO="${SCENARIO:?set SCENARIO to true_ne or estimated_ne}"
module load conda
conda activate dnaging
python "$PROJECT/summarize_variable_ne_posterior_insertion.py" --root "$PROJECT/msprime_variable_ne_error" \
  --scenario "$SCENARIO" --n-replicates "${N_REPLICATES:-10}" \
  --output "$PROJECT/msprime_variable_ne_error/insertion_posterior/$SCENARIO/summary"
