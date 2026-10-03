#!/bin/bash
set -euo pipefail
PROJECT="${PROJECT:-$PWD}"
N="${N_REPLICATES:-10}"
KEEP="${KEEP_LAST:-50}"
if [[ "$KEEP" -ne 50 ]]; then echo "ERROR: this suite is configured for KEEP_LAST=50" >&2; exit 2; fi
for ((rep=1; rep<=N; rep++)); do
  r=$(printf '%03d' "$rep")
  base="$PROJECT/msprime_variable_ne_error/argtest/replicate_$r"
  count=$(find "$base/out/combined" -maxdepth 1 -type f -name '*.tsz' | wc -l)
  [[ "$count" -eq 100 ]] || { echo "ERROR replicate_$r has $count/100 draws" >&2; exit 3; }
  [[ -s "$base/ne/estimated_ne.tsv" ]] || { echo "ERROR missing replicate_$r/ne/estimated_ne.tsv" >&2; exit 4; }
  [[ -s "$PROJECT/msprime_variable_ne_error/simulations/replicate_$r/constant_ne_epochs.tsv" ]] || \
    { echo "ERROR missing replicate_$r true Ne trajectory" >&2; exit 5; }
done
prep=$(sbatch --parsable --array="1-$N%5" "$PROJECT/slurm/prepare_variable_ne_insertion.sbatch")
draws=$((N*KEEP-1))
printf 'prepare=%s\n' "$prep"
for scenario in true_ne estimated_ne; do
  infer=$(sbatch --parsable --dependency="afterok:$prep" --array="0-$draws%50" \
    --export="ALL,SCENARIO=$scenario" "$PROJECT/slurm/run_variable_ne_insertion_draws.sbatch")
  merge=$(sbatch --parsable --dependency="afterok:$infer" --array="1-$N%10" \
    --export="ALL,SCENARIO=$scenario" "$PROJECT/slurm/merge_variable_ne_insertion_draws.sbatch")
  summary=$(sbatch --parsable --dependency="afterok:$merge" \
    --export="ALL,N_REPLICATES=$N,SCENARIO=$scenario" "$PROJECT/slurm/summarize_variable_ne_insertion.sbatch")
  printf '%s_inference=%s\n%s_merge=%s\n%s_summary=%s\n' \
    "$scenario" "$infer" "$scenario" "$merge" "$scenario" "$summary"
done
