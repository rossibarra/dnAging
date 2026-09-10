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
#      All of these are optional and pass-if-set; the Python CLI owns the
#      defaults, so nothing is duplicated here:
#        EPSILON  PLOIDY  MIN_N  MUT_AGE_MAX  MARGINALISE  PER_SAMPLE_TSV=1
#        SAVE_EPSILON_DATA=1  DRAW_IDS=0,8,16

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

# ${EXTRA[@]+...} not "${EXTRA[@]}": expanding an EMPTY array under `set -u` is an
# unbound-variable error on bash before 4.4, which kills the job when no optional
# setting is passed.
EXTRA=()
[[ -n "${INCLUDE:-}" ]]       && EXTRA+=(--include-positions "$INCLUDE")
[[ -n "${SAMPLES_FILE:-}" ]]  && EXTRA+=(--samples-file "$SAMPLES_FILE")
# Everything below is pass-if-set, never defaulted here: the Python CLI owns the
# defaults so they cannot drift out of sync with values duplicated in this file.
# Set them explicitly for a run whose settings must be pinned (see
# HPC_REAL_DATA_CODEX.md), and they are then recorded in each part's run.json.
[[ -n "${EPSILON:-}" ]]       && EXTRA+=(--epsilon "$EPSILON")
[[ -n "${PLOIDY:-}" ]]        && EXTRA+=(--ploidy "$PLOIDY")
[[ -n "${MIN_N:-}" ]]         && EXTRA+=(--min-n "$MIN_N")
[[ -n "${MUT_AGE_MAX:-}" ]]   && EXTRA+=(--mutation-age-max "$MUT_AGE_MAX")
[[ -n "${MARGINALISE:-}" ]]   && EXTRA+=(--marginalise "$MARGINALISE")
[[ "${PER_SAMPLE_TSV:-0}" == "1" ]] && EXTRA+=(--per-sample-tsv)
[[ "${SAVE_EPSILON_DATA:-0}" == "1" ]] && EXTRA+=(--save-epsilon-data)
[[ -n "${DRAW_IDS:-}" ]]       && EXTRA+=(--draw-ids "$DRAW_IDS")

CHROMS_FILE="${CHROMS_FILE:-chroms.txt}"
[[ -f "$CHROMS_FILE" ]] || { echo "ERROR: no such file: $CHROMS_FILE" >&2; exit 91; }
# read with a portable loop rather than `mapfile`, which needs bash 4 and is
# absent on bash 3.2 (macOS, and some older cluster images)
CHROMS=()
while IFS= read -r line; do CHROMS+=("$line"); done < <(awk 'NF && $1 !~ /^#/' "$CHROMS_FILE")
(( ${#CHROMS[@]} )) || { echo "ERROR: no chromosomes in $CHROMS_FILE" >&2; exit 92; }

if [[ "${MERGE:-0}" == "1" ]]; then
    # Derive the parts from CHROMS_FILE rather than globbing "$OUTROOT"/chr*.
    # The glob silently matched nothing whenever chromosomes are named 1..22
    # instead of chr1..chr22 -- and with `set -u` the unmatched pattern was then
    # passed through literally. Requiring every part to exist also enforces
    # "do not merge partial output", which a glob cannot.
    PARTS=(); MISSING=()
    for c in "${CHROMS[@]}"; do
        if [[ -d "$OUTROOT/$c" ]]; then PARTS+=("$OUTROOT/$c"); else MISSING+=("$c"); fi
    done
    if (( ${#MISSING[@]} )); then
        echo "ERROR: ${#MISSING[@]} of ${#CHROMS[@]} chromosome parts missing under $OUTROOT:" >&2
        printf '  %s\n' "${MISSING[@]}" >&2
        echo "Rerun those array indices; refusing to merge partial output." >&2
        exit 93
    fi
    echo "merging ${#PARTS[@]} chromosome parts"
    python "$PROJECT/posterior_sample_age_infer.py" \
        --freq-table "$TABLE" --merge ${PARTS[@]+"${PARTS[@]}"} \
        --output "$OUTROOT/genome"
    exit 0
fi

CHROM="${CHROMS[${SLURM_ARRAY_TASK_ID:-0}]}"
[[ -n "$CHROM" ]] || { echo "ERROR: no chromosome at array index ${SLURM_ARRAY_TASK_ID:-0}" >&2; exit 94; }

python "$PROJECT/posterior_sample_age_infer.py" \
    --freq-table "$TABLE" \
    --store "$STORE" --draw-polarity "$POLARITY" \
    --panel-vcf "$PANELVCF" --vcf "$VCF" --chrom "$CHROM" \
    ${EXTRA[@]+"${EXTRA[@]}"} \
    --output "$OUTROOT/$CHROM"
