#!/bin/bash
#SBATCH --job-name=Tinfer
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/Tinfer_%A_%a.out
#SBATCH --error=logs/Tinfer_%A_%a.err
#
# STEP 2: per-chromosome inference for ALL ancient samples at once (one multi-
# sample VCF), by table lookup on (t_i from the store, d0 from the panel VCF).
# STEP 3: merge chromosomes into the genome-wide per-sample posteriors.
#
#   # per chromosome (array over chroms.txt):
#   TABLE=freq_table.npz STORE=... POLARITY=... PANELVCF=panel26.vcf.gz \
#   VCF=ancient_300.vcf.gz OUTROOT=results/Tage \
#     sbatch --array=0-9 slurm/run_infer.sbatch
#
#   # combine chromosomes -> results/Tage/genome/ages_table.tsv:
#   MERGE=1 TABLE=freq_table.npz OUTROOT=results/Tage \
#     sbatch slurm/run_infer.sbatch
#
# Env: TABLE STORE POLARITY PANELVCF VCF OUTROOT   (INCLUDE optional site list)
#      CHROMS_FILE (default chroms.txt)  EPSILON (default 1e-3)  SAMPLES_FILE (optional)

set -euo pipefail
PROJECT="${PROJECT:-$SLURM_SUBMIT_DIR}"
source "$PROJECT/slurm/slurm_conda_bootstrap.sh" 2>/dev/null || true
mkdir -p logs

EXTRA=()
[[ -n "${INCLUDE:-}" ]]       && EXTRA+=(--include-positions "$INCLUDE")
[[ -n "${SAMPLES_FILE:-}" ]]  && EXTRA+=(--samples-file "$SAMPLES_FILE")

if [[ "${MERGE:-0}" == "1" ]]; then
    python "$PROJECT/posterior_sample_age_infer.py" \
        --freq-table "$TABLE" --merge "$OUTROOT"/chr* \
        --output "$OUTROOT/genome"
    exit 0
fi

CHROMS_FILE="${CHROMS_FILE:-chroms.txt}"
mapfile -t CHROMS < "$CHROMS_FILE"
CHROM="${CHROMS[${SLURM_ARRAY_TASK_ID:-0}]}"

python "$PROJECT/posterior_sample_age_infer.py" \
    --freq-table "$TABLE" \
    --store "$STORE" --draw-polarity "$POLARITY" \
    --panel-vcf "$PANELVCF" --vcf "$VCF" --chrom "$CHROM" \
    --epsilon "${EPSILON:-1e-3}" \
    "${EXTRA[@]}" \
    --output "$OUTROOT/$CHROM"
