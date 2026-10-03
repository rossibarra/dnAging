#!/bin/bash -l
#SBATCH --job-name=freqCalls
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-10%10
set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"; module load conda; conda activate dnaging
python "$PROJECT/prepare_frequency_stratified_calls.py" --root "$PROJECT/msprime_variable_ne_error" --replicate "$SLURM_ARRAY_TASK_ID"
