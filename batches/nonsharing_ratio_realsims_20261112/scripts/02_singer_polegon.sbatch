#!/bin/bash -l
# SINGER (150 MCMC samples, thin 10) + POLEGON on the 26 panel haplotypes of ONE replicate (sbatch --export=REP=n).
# Same snakemake recipe/settings as batches/nonsharing_ratio_heldout_20261111 (SINGER 0.1.8-beta custom, POLEGON 0.1.3-alpha).
# Only the 5 draws in DRAWS are merged; SINGER runs the whole chain.  Resumable: resubmit after preemption.
#SBATCH --job-name=nsrSim
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --requeue
#SBATCH --output=logs/nsrSim-%j.out
#SBATCH --error=logs/nsrSim-%j.err
set -euo pipefail
B=/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nonsharing_ratio_realsims_20261112
SM=/quobyte/jrigrp/beil/te_evo/singer_analysis/singer-snakemake
CHR=sim$REP
DRAWS="75 93 112 131 149"
source /etc/profile.d/modules.sh
module load conda
conda activate singer-snakemake-env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$B"
exec 9>"$B/singer_out/.resume.guard.$REP"
flock -n 9 || { echo "another writer owns singer_out/sim$REP" >&2; exit 1; }
TARGETS=(); for r in $DRAWS; do TARGETS+=("$B/singer_out/$CHR/$CHR/trees/$CHR.$r.tsz"); done
mkdir -p "$B/singer_out/$CHR"
if [[ -d "$B/singer_out/$CHR/.snakemake/locks" ]]; then
  snakemake --snakefile "$SM/workflow/Snakefile" --directory "$B/singer_out/$CHR" --configfile "$B/singer_config_rep$REP.yaml" --unlock || true
fi
snakemake --snakefile "$SM/workflow/Snakefile" --directory "$B/singer_out/$CHR" \
  --configfile "$B/singer_config_rep$REP.yaml" --cores "$SLURM_CPUS_PER_TASK" --rerun-incomplete \
  --keep-going "${TARGETS[@]}"
