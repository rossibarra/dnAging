#!/bin/bash -l
#SBATCH --job-name=vneNe
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-10%5

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ARGTEST_REPO="${ARGTEST_REPO:-/home/jri/src/argtest}"
REP=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
WORK="$PROJECT/msprime_variable_ne_error/argtest/replicate_${REP}"
SIM="$PROJECT/msprime_variable_ne_error/simulations/replicate_${REP}"
KEEP="${KEEP_LAST:-50}"

module load conda
conda activate argtest

# Use only the LAST $KEEP MCMC reps. burnin-frac is a fraction, so derive it
# from the actual file count rather than assuming 100 combined ARGs: with
# exactly 100 the default 0.5 coincides with this, but a partial merge would
# silently keep the wrong number. find_tree_files sorts numerically, so
# run.combined.<n>.tsz orders 0..99 and "last" means highest MCMC index.
N=$(ls "$WORK/out/combined"/*.tsz 2>/dev/null | wc -l)
if [[ "$N" -lt "$KEEP" ]]; then
    echo "ERROR: only $N combined ARGs, need at least $KEEP" >&2
    exit 4
fi
FRAC=$(python3 -c "print(($N - $KEEP) / $N)")
echo "replicate $REP: $N combined ARGs, burnin-frac $FRAC keeps the last $KEEP"

# --num-bins/--num-quantiles/--time-bins-file has no default: the script
# raises unless exactly one is given. 50 log-spaced windows is the convention
# MATH.md records for the real analysis.
python "$ARGTEST_REPO/scripts/coalescence_ne_plots_from_ts.py" \
    --ts-dir "$WORK/out/combined" \
    --pattern "*.tsz" \
    --burnin-frac "$FRAC" \
    --num-bins "${NUM_BINS:-50}" \
    --out-dir "$WORK/ne"

# Convert to the schema insertion's PiecewiseConstantNe.from_tsv accepts. Uses
# the project env because it imports insertion_likelihood to validate the result.
conda activate dnaging
python "$PROJECT/validation/argtest_ne_to_insertion.py" \
    --estimates "$WORK/ne/coalescence-ne-estimates.tsv" \
    --output "$WORK/ne/estimated_ne.tsv" \
    --truth "$SIM/constant_ne_epochs.tsv"
