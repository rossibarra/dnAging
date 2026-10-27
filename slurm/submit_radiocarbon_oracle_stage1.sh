#!/bin/bash
set -euo pipefail

PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
RUN_ROOT="${RUN_ROOT:-$PROJECT/batches/nine_demographies_20260916/runs/constant_r1}"
OUTROOT="${OUTROOT:-$PROJECT/batches/radiocarbon_oracle_stage1_20261011/constant_r1/injected_0.01}"
INJECTED_EPSILON="${INJECTED_EPSILON:-0.01}"
INJECTION_SEED="${INJECTION_SEED:-2026101101}"

array_job=$(sbatch --parsable \
  --export="ALL,RUN_ROOT=$RUN_ROOT,OUTROOT=$OUTROOT,INJECTED_EPSILON=$INJECTED_EPSILON,INJECTION_SEED=$INJECTION_SEED" \
  "$PROJECT/slurm/run_radiocarbon_oracle_stage1.sbatch")
merge_job=$(sbatch --parsable --dependency="afterok:$array_job" \
  --export="ALL,OUTROOT=$OUTROOT" \
  "$PROJECT/slurm/merge_radiocarbon_oracle_stage1.sbatch")
printf 'array_job\t%s\nmerge_job\t%s\noutput\t%s\n' "$array_job" "$merge_job" "$OUTROOT/result"
