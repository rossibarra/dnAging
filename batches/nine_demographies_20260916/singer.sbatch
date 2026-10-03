#!/bin/bash -l
#SBATCH --job-name=nine_singer
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=72G
#SBATCH --time=3-00:00:00
#SBATCH --output=nine-singer-%j.out

set -euo pipefail
export PROJECT=/quobyte/jrigrp/beil/logan_collab/dnAging
export BATCH_ROOT=/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nine_demographies_20260916
export PYTHONPATH="$BATCH_ROOT/code:$PROJECT:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load conda
conda activate dnaging
RUN_INDEX="${RUN_INDEX:?}"
RUN_LINE=$(awk -F '\t' -v n="$((RUN_INDEX+2))" 'NR==n' "$BATCH_ROOT/runs.tsv")
IFS=$'\t' read -r INDEX RUN_NAME MODEL REPLICATE RUN_ROOT BASE_SEED SINGER_SEED <<< "$RUN_LINE"
export RUN_ROOT
test "$INDEX" = "$RUN_INDEX"
cd "$RUN_ROOT"

python "$BATCH_ROOT/code/validate_samples.py" "$RUN_ROOT"
conda activate singer-snakemake-env
# A single repair/resume writer per workflow. The caller must first verify that
# no older, unguarded workflow job is active in this directory.
exec 9>"$RUN_ROOT/singer/.resume.guard"
flock -n 9 || { echo "Another resume process owns this workflow" >&2; exit 1; }
if [[ -d "$RUN_ROOT/singer/.snakemake/locks" ]]; then
 cp -an "$RUN_ROOT/singer/.snakemake/locks" "$RUN_ROOT/singer/.snakemake/locks.pre_unlock_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}.bak"
fi
snakemake --snakefile /home/jri/src/singer-snakemake/workflow/Snakefile \
 --directory "$RUN_ROOT/singer" --configfile "$RUN_ROOT/singer/config.yaml" --unlock
snakemake --snakefile /home/jri/src/singer-snakemake/workflow/Snakefile \
 --directory "$RUN_ROOT/singer" --configfile "$RUN_ROOT/singer/config.yaml" \
 --cores "$SLURM_CPUS_PER_TASK" --rerun-incomplete
