#!/bin/bash -l
#SBATCH --job-name=b300god
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=0-299%40

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="${ROOT:-$PROJECT/msprime_batch300_insertion}"
REP=$((SLURM_ARRAY_TASK_ID + 1))
module load conda
conda activate dnaging
# n_modern varies across batch300 (panel 10-40 haploids), so read it per replicate.
NMODERN=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["n_modern_haploid"])' \
    "$ROOT/simulations/replicate_$(printf '%03d' "$REP")/metadata.json")
python "$PROJECT/godambe_block_insertion.py" \
    --input-root "$ROOT/simulations" \
    --output-dir "$ROOT/godambe" \
    --replicate "$REP" --n-modern "$NMODERN" \
    --n-blocks "${N_BLOCKS:-50}" \
    --age-min 0 --age-max 15000 --age-step 50
