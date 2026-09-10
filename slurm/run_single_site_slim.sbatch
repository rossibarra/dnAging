#!/bin/bash -l
#SBATCH --job-name=slim1site
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --array=1-10000%1000

set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
MODEL="${MODEL:-$PROJECT/slim/single_site_neutral.slim}"
OUTDIR="${OUTDIR:-$PROJECT/slim_single_site_output}"
NE="${NE:-10000}"
MUT_RATE="${MUT_RATE:-1e-6}"
END_GEN="${END_GEN:-40000}"
BASE_SEED="${BASE_SEED:-202609090000}"
REPLICATE="${SLURM_ARRAY_TASK_ID:?This script must be submitted as an array job}"
SEED=$((BASE_SEED + REPLICATE))

module load conda
conda activate dnaging
command -v slim >/dev/null 2>&1 || {
    echo "ERROR: slim is unavailable; update dnaging from environment.yml" >&2
    exit 90
}

mkdir -p "$OUTDIR/frequencies" "$OUTDIR/samples" "$OUTDIR/logs"

stem="replicate_$(printf '%06d' "$REPLICATE").seed_${SEED}"
freq_final="$OUTDIR/frequencies/${stem}.frequencies.tsv"
sample_final="$OUTDIR/samples/${stem}.samples.tsv"
freq_tmp="${freq_final}.tmp.${SLURM_JOB_ID}"
sample_tmp="${sample_final}.tmp.${SLURM_JOB_ID}"
slim_log="$OUTDIR/logs/${stem}.slim.log"

if [[ -e "$freq_final" || -e "$sample_final" ]]; then
    echo "ERROR: completed output already exists for replicate $REPLICATE and seed $SEED" >&2
    exit 3
fi

cleanup()
{
    rm -f "$freq_tmp" "$sample_tmp"
}
trap cleanup EXIT

slim \
    -s "$SEED" \
    -d "NE=${NE}" \
    -d "MUT_RATE=${MUT_RATE}" \
    -d "END_GEN=${END_GEN}" \
    -d "REPLICATE=${REPLICATE}" \
    -d "FREQ_FILE=\"${freq_tmp}\"" \
    -d "SAMPLE_FILE=\"${sample_tmp}\"" \
    "$MODEL" >"$slim_log" 2>&1

# Only completed replicates acquire the mergeable final suffixes.
mv "$freq_tmp" "$freq_final"
mv "$sample_tmp" "$sample_final"
trap - EXIT

printf 'completed replicate=%s seed=%s frequency_file=%s sample_file=%s\n' \
    "$REPLICATE" "$SEED" "$freq_final" "$sample_final"
