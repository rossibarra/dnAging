#!/bin/bash
#SBATCH --job-name=Tprecmp
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/Tprecmp_%j.out
#SBATCH --error=logs/Tprecmp_%j.err
#
# STEP 1 (once): build the age-conditioned frequency table E[p_T | d0, t_i, n] from
# the coalescence-Ne curve, via the exact neutral WF moment recursion. Single core,
# demography-specific and independent of samples/sites, so it is built once and
# reused for every chromosome and every ancient sample.
#
# COST. No longer the old float64 cost model. Precision is chosen per age row and
# both extremes of the mutation-age axis are expensive, so cost is dominated by the
# OLDEST rows: at n=26 one (panel, 300 sample-age) row takes ~0.02 s at tau_i=5e-4
# and ~2.9 s at tau_i=2.3, but ~42 s at tau_i=2000 where ~930 digits are needed.
# The 2-day reservation above is sized for that tail across MIN_N..NSAMPLE panels.
#
# The cheapest way to cut it is AGEMAX. Inference discards mutation-age mass beyond
# --mutation-age-max (default tau_i = 3), so rows far past that are computed and
# never read. Set AGEMAX just beyond the generation age corresponding to tau_i = 3
# under this demography rather than leaving the historical 4e7 default, which at
# Ne ~ 1e4 reaches tau_i ~ 2000 and spends most of the build on unused rows.
#
#   NE=coalescence-ne-estimates.tsv OUT=freq_table.npz \
#     sbatch slurm/run_precompute.sbatch
#
# Env:  NE (required)  OUT (required)
#       NSAMPLE (default 26)  TMAX (default 30000)  NT (default 300)
#       AGEMAX (default 4e7)  NAGE (default 100)
#       MIN_N, PRECISION  (optional, pass-if-set; the Python CLI owns the defaults)

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

# pass-if-set, so no default is duplicated out of the Python CLI
# ${EXTRA[@]+...} not "${EXTRA[@]}": expanding an EMPTY array under `set -u` is an
# unbound-variable error on bash before 4.4, which kills the job when no optional
# setting is passed.
EXTRA=()
[[ -n "${MIN_N:-}" ]]     && EXTRA+=(--min-n "$MIN_N")
[[ -n "${PRECISION:-}" ]] && EXTRA+=(--precision "$PRECISION")

python "$PROJECT/precompute_freq_trajectory_moments.py" \
    --ne "$NE" \
    --n-sample "${NSAMPLE:-26}" \
    --t-min 0 --t-max "${TMAX:-30000}" --n-t "${NT:-300}" \
    --age-min 10 --age-max "${AGEMAX:-4e7}" --n-age "${NAGE:-100}" \
    ${EXTRA[@]+"${EXTRA[@]}"} \
    --output "$OUT"
