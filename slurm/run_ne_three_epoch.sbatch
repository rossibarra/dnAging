#!/bin/bash -l
#SBATCH --job-name=3epNe
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
W="$PROJECT/msprime_three_epoch/argtest"
KEEP="${KEEP_LAST:-50}"

module load conda
conda activate argtest
# Only the LAST $KEEP MCMC reps. burnin-frac is a fraction, so derive it from the
# actual file count rather than assuming 100; find_tree_files sorts numerically,
# so run.combined.<n>.tsz orders 0..99 and "last" means highest MCMC index.
N=$(ls "$W/out/combined"/*.tsz | wc -l)
if [[ "$N" -lt "$KEEP" ]]; then echo "ERROR: only $N ARGs, need $KEEP" >&2; exit 4; fi
FRAC=$(python3 -c "print(($N - $KEEP) / $N)")
echo "using burnin-frac $FRAC of $N ARGs -> last $KEEP"
# --num-bins has no default; 50 log-spaced windows is the MATH.md convention.
python /home/jri/src/argtest/scripts/coalescence_ne_plots_from_ts.py \
    --ts-dir "$W/out/combined" --pattern "*.tsz" \
    --burnin-frac "$FRAC" --num-bins "${NUM_BINS:-50}" \
    --out-dir "$W/ne"

conda activate dnaging
python "$PROJECT/argtest_ne_to_insertion.py" \
    --estimates "$W/ne/coalescence-ne-estimates.tsv" \
    --output "$W/ne/estimated_ne.tsv" \
    --truth "$PROJECT/msprime_three_epoch/simulations/replicate_001/constant_ne_epochs.tsv"
