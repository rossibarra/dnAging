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
#      CHROMS_FILE (default chroms.txt)  SAMPLES_FILE (optional)
#      EPSILON (optional; defaults to the value in posterior_sample_age_infer.py)

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
mkdir -p logs

# Activate the project conda env BEFORE submitting; sbatch defaults to
# --export=ALL, so the job inherits it. Checked here so a forgotten activation
# fails with this message instead of an opaque ImportError further down.
python -c 'import numpy, normalize_tes' 2>/dev/null || {
    echo "ERROR: 'python' cannot import numpy, normalize_tes." >&2
    echo "  normalize_tes is NOT in the conda env; it comes from the normalizeTEs" >&2
    echo "  checkout on PYTHONPATH. Before submitting:" >&2
    echo "    conda activate dnaging" >&2
    echo "    export PYTHONPATH=/path/to/normalizeTEs:\${PYTHONPATH:-}" >&2
    echo "    sbatch $0" >&2
    exit 90
}

EXTRA=()
[[ -n "${INCLUDE:-}" ]]       && EXTRA+=(--include-positions "$INCLUDE")
[[ -n "${SAMPLES_FILE:-}" ]]  && EXTRA+=(--samples-file "$SAMPLES_FILE")
# Only pass --epsilon when set: the Python CLI owns the default, so it cannot
# drift out of sync with a value duplicated here.
[[ -n "${EPSILON:-}" ]]       && EXTRA+=(--epsilon "$EPSILON")

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
    "${EXTRA[@]}" \
    --output "$OUTROOT/$CHROM"
