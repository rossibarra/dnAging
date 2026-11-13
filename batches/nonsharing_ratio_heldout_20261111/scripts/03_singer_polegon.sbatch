#!/bin/bash -l
# SINGER (150 MCMC samples, thin 10) + POLEGON dating, chr10, panel minus held-out lines.
# Recipe = the real NAM+TE run (singer-snakemake from beil/te_evo, singer-0.1.8-beta custom, polegon-0.1.3-alpha),
# config = chr.9.yaml values (see singer_config.yaml).  Only the 10 post-burn-in draws in DRAWS are merged;
# SINGER itself must run the full chain.  Snakemake is resumable: re-submit this script after a preemption.
#SBATCH --job-name=nsrSinger
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --requeue
#SBATCH --output=logs/nsrSinger-%j.out
#SBATCH --error=logs/nsrSinger-%j.err
set -euo pipefail
B=/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nonsharing_ratio_heldout_20261111
SM=/quobyte/jrigrp/beil/te_evo/singer_analysis/singer-snakemake
CHR=chr.10.combined.snp.te.sorted
DRAWS="75 82 89 96 103 110 117 124 131 138"
source /etc/profile.d/modules.sh
module load conda
conda activate singer-snakemake-env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$B/singer_out"
cd "$B"
exec 9>"$B/singer_out/.resume.guard"
flock -n 9 || { echo "another writer owns singer_out" >&2; exit 1; }
TARGETS=()
for r in $DRAWS; do TARGETS+=("$B/singer_out/$CHR/trees/$CHR.$r.tsz"); done
if [[ -d "$B/singer_out/.snakemake/locks" ]]; then
  snakemake --snakefile "$SM/workflow/Snakefile" --directory "$B/singer_out" --configfile "$B/singer_config.yaml" --unlock || true
fi
snakemake --snakefile "$SM/workflow/Snakefile" --directory "$B/singer_out" \
  --configfile "$B/singer_config.yaml" --cores "$SLURM_CPUS_PER_TASK" --rerun-incomplete \
  --keep-going "${TARGETS[@]}"
