#!/bin/bash -l
set -euo pipefail
cd "$(dirname "$0")"
BATCH_ROOT="$PWD"
test ! -e jobs.tsv || { echo "jobs.tsv exists; refusing duplicate submissions" >&2; exit 1; }
printf 'run\tstage\tjob_id\tdependency\n' > jobs.tsv
submit() {
 local stage="$1" dependency="$2" id
 local extra=(--partition=low --qos=jrigrp-low-qos)
 if [[ -n "$dependency" ]]; then extra+=(--dependency="afterok:$dependency"); fi
 id=$(sbatch --parsable --export="ALL,RUN_INDEX=$idx" --chdir="$BATCH_ROOT/runs/$name" \
     --output="$BATCH_ROOT/runs/$name/logs/$stage-%A_%a.out" "${extra[@]}" "$BATCH_ROOT/$stage.sbatch") || return 1
 id="${id%%;*}"
 [[ "$id" =~ ^[0-9]+$ ]] || { echo "Invalid sbatch job ID: $id" >&2; return 1; }
 printf '%s\t%s\t%s\t%s\n' "$name" "$stage" "$id" "$dependency" >> jobs.tsv
 printf '%s' "$id"
}
while IFS=$'\t' read -r idx name model rep root seed singer_seed; do
 [[ "$idx" == index ]] && continue
 sim=$(submit simulate "")
 singer=$(submit singer "$sim")
 argtest=$(submit argtest "$singer")
 ne=$(submit ne "$argtest")
 mask=$(submit mask "$argtest")
 prep=$(submit prepare "$mask")
 insertion=$(submit insertion "$ne:$prep")
 merge=$(submit merge "$insertion")
 echo "$name: sim=$sim singer=$singer argtest=$argtest ne=$ne mask=$mask prepare=$prep insertion=$insertion merge=$merge"
done < runs.tsv
