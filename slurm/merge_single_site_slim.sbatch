#!/bin/bash -l
#SBATCH --job-name=merge1site
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
OUTDIR="${OUTDIR:-$PROJECT/slim_single_site_output}"
EXPECTED_REPLICATES="${EXPECTED_REPLICATES:-100}"

shopt -s nullglob
frequency_files=("$OUTDIR"/frequencies/*.frequencies.tsv)
sample_files=("$OUTDIR"/samples/*.samples.tsv)

if (( ${#frequency_files[@]} != EXPECTED_REPLICATES )); then
    echo "ERROR: expected $EXPECTED_REPLICATES frequency files, found ${#frequency_files[@]}" >&2
    exit 2
fi
if (( ${#sample_files[@]} != EXPECTED_REPLICATES )); then
    echo "ERROR: expected $EXPECTED_REPLICATES sample files, found ${#sample_files[@]}" >&2
    exit 2
fi

if [[ -e "$OUTDIR/all_replicates.frequencies.tsv" || \
      -e "$OUTDIR/all_replicates.samples.tsv" ]]; then
    echo "ERROR: merged output already exists in $OUTDIR" >&2
    exit 3
fi

merge_tsvs()
{
    local output="$1"
    shift
    local files=("$@")
    local tmp="${output}.tmp.${SLURM_JOB_ID}"

    awk 'FNR == 1 && NR != 1 {next} {print}' "${files[@]}" >"$tmp"
    mv "$tmp" "$output"
}

merge_tsvs "$OUTDIR/all_replicates.frequencies.tsv" "${frequency_files[@]}"
merge_tsvs "$OUTDIR/all_replicates.samples.tsv" "${sample_files[@]}"

printf 'merged %d replicates into %s\n' \
    "$EXPECTED_REPLICATES" "$OUTDIR"
