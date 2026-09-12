#!/bin/bash -l
#SBATCH --job-name=slim-edge-consolidate
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
MAIN="$PROJECT/slim_single_site_edge_validation"
RETRY="$PROJECT/slim_single_site_edge_validation_retry"
BASE_SEED=202609090000

for i in {323..384} {1148..1152} {1283..1290} {1296..1303}; do
    stem="replicate_$(printf '%06d' "$i").seed_$((BASE_SEED + i))"
    retry_tree="$RETRY/trees/${stem}.trees"
    retry_sample="$RETRY/samples/${stem}.samples.tsv"
    main_tree="$MAIN/trees/${stem}.trees"
    main_sample="$MAIN/samples/${stem}.samples.tsv"
    [[ -s "$retry_tree" && -s "$retry_sample" ]] || {
        echo "ERROR: incomplete retry output for $stem" >&2
        exit 4
    }
    [[ ! -e "$main_tree" ]] || {
        echo "ERROR: refusing to replace existing tree $main_tree" >&2
        exit 5
    }
    if [[ -e "$main_sample" ]]; then
        [[ ! -e "${main_sample}.bak" ]] || {
            echo "ERROR: backup already exists for $main_sample" >&2
            exit 6
        }
        mv "$main_sample" "${main_sample}.bak"
    fi
    cp -p "$retry_tree" "$main_tree"
    cp -p "$retry_sample" "$main_sample"
done

tree_count=$(find "$MAIN/trees" -maxdepth 1 -type f -name '*.trees' | wc -l)
sample_count=$(find "$MAIN/samples" -maxdepth 1 -type f -name '*.samples.tsv' | wc -l)
[[ "$tree_count" -eq 10000 && "$sample_count" -eq 10000 ]] || {
    echo "ERROR: expected 10000 pairs, found trees=$tree_count samples=$sample_count" >&2
    exit 7
}
echo "validated trees=$tree_count samples=$sample_count"
