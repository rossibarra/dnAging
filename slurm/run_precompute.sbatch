#!/bin/bash
#SBATCH --job-name=Tprecmp
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/Tprecmp_%j.out
#SBATCH --error=logs/Tprecmp_%j.err
#
# STEP 1 (once): build the age-conditioned frequency table E[p_T | d0, t_i] from
# the coalescence-Ne curve, via the exact neutral WF moment recursion. Cheap and
# demography-specific; independent of samples/sites, so it is built a single time
# and reused for every chromosome and every ancient sample.
#
#   NE=coalescence-ne-estimates.tsv OUT=freq_table.npz \
#     sbatch slurm/run_precompute.sbatch
#
# Env:  NE (required)  OUT (required)
#       NSAMPLE (default 26)  TMAX (default 30000)  NT (default 300)
#       AGEMAX (default 4e7)  NAGE (default 100)

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
mkdir -p logs

# Activate the project conda env BEFORE submitting; sbatch defaults to
# --export=ALL, so the job inherits it. Checked here so a forgotten activation
# fails with this message instead of an opaque ImportError further down.
python -c 'import numpy, scipy' 2>/dev/null || {
    echo "ERROR: 'python' cannot import numpy, scipy." >&2
    echo "  Activate the env before submitting: conda activate dnaging; sbatch $0" >&2
    exit 90
}

python "$PROJECT/precompute_freq_trajectory_moments.py" \
    --ne "$NE" \
    --n-sample "${NSAMPLE:-26}" \
    --t-min 0 --t-max "${TMAX:-30000}" --n-t "${NT:-300}" \
    --age-min 10 --age-max "${AGEMAX:-4e7}" --n-age "${NAGE:-100}" \
    --output "$OUT"
