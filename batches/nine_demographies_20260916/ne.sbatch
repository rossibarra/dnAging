#!/bin/bash -l
#SBATCH --job-name=nine_ne
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail
export PROJECT=/quobyte/jrigrp/beil/logan_collab/dnAging
export BATCH_ROOT=/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nine_demographies_20260916
export PYTHONPATH="$BATCH_ROOT/code:$PROJECT:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load conda
conda activate dnaging
RUN_INDEX="${RUN_INDEX:?}"
RUN_LINE=$(awk -F '\t' -v n="$((RUN_INDEX+2))" 'NR==n' "$BATCH_ROOT/runs.tsv")
IFS=$'\t' read -r INDEX RUN_NAME MODEL REPLICATE RUN_ROOT BASE_SEED SINGER_SEED <<< "$RUN_LINE"
export RUN_ROOT
test "$INDEX" = "$RUN_INDEX"
cd "$RUN_ROOT"
W="$RUN_ROOT/argtest"
KEEP="${KEEP_LAST:-50}"

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
python "$BATCH_ROOT/code/argtest_ne_to_insertion.py" \
    --estimates "$W/ne/coalescence-ne-estimates.tsv" \
    --output "$W/ne/estimated_ne.tsv" \
    --truth "$RUN_ROOT/simulations/replicate_001/constant_ne_epochs.tsv"
