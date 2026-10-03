#!/bin/bash -l
#SBATCH --job-name=freqTrueARG
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=1-10%10
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"; R=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
module load conda; conda activate dnaging
python "$PROJECT/insertion_frequency_stratified_true_arg.py" --root "$PROJECT/msprime_variable_ne_error" --replicate "$SLURM_ARRAY_TASK_ID" --epsilon "${EPSILON:-0.01}" --output "$PROJECT/msprime_variable_ne_error/insertion_frequency/true_arg/replicate_$R"
