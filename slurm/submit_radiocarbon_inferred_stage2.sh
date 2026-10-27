#!/bin/bash
set -euo pipefail
PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
RUN_ROOT="${RUN_ROOT:-$PROJECT/batches/nine_demographies_20260916/runs/constant_r1}"
OUTROOT="${OUTROOT:-$PROJECT/batches/radiocarbon_inferred_stage2_20261012/constant_r1}"
array_job=$(sbatch --parsable --export="ALL,RUN_ROOT=$RUN_ROOT,OUTROOT=$OUTROOT" \
  "$PROJECT/slurm/run_radiocarbon_inferred_stage2.sbatch")
merge_job=$(sbatch --parsable --dependency="afterok:$array_job" \
  --export="ALL,RUN_ROOT=$RUN_ROOT,OUTROOT=$OUTROOT" \
  "$PROJECT/slurm/merge_radiocarbon_inferred_stage2.sbatch")
printf 'array_job\t%s\nmerge_job\t%s\noutput\t%s\n' "$array_job" "$merge_job" "$OUTROOT/result"
