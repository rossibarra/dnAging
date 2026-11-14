#!/bin/bash -l
# run-draw for every draw x A0 of ONE replicate (sbatch --export=ALL,REP=n), then merge per (rep, A0).
# One writer per output file; each rep is its own job.  Needs the 5 POLEGON trees of that rep.
#SBATCH --job-name=nsrRun
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=logs/nsrRun-%j.out
#SBATCH --error=logs/nsrRun-%j.err
set -euo pipefail
source /etc/profile.d/modules.sh
module load conda
conda activate dnaging
cd /quobyte/jrigrp/beil/logan_collab/dnAging
B=batches/nonsharing_ratio_realsims_bottleneck_20261112
CLI="python pipeline/nonsharing_ratio_real_data.py"
r=$REP
for d in 75 93 112 131 149; do
  for a in 10000 20000 50000; do
    out=$B/draws/rep$r/draw${d}_a0$a.npz
    [[ -e $out ]] && continue
    $CLI run-draw --tree $B/singer_out/sim$r/sim$r/trees/sim$r.$d.tsz --draw-id $d \
      --controls-vcf $B/sim/rep$r/controls.vcf.gz --controls-mode sample --ancient-vcf $B/sim/rep$r/ancient.vcf.gz \
      --panel-names $B/sim/rep$r/panel_names.txt --chrom-fai $B/chr1.fai --a0 $a --age-mode child \
      --panel-fixed $B/draws/rep$r/panel_fixed.npz --fixed-polarity ref --output $out > $B/logs/draw_rep${r}_d${d}_a0$a.log 2>&1
  done
done
# same draws but with c(t) from the TRUE demography (demes file, haploid counts -> --ploidy 1)
for d in 75 93 112 131 149; do
  for a in 10000 20000 50000; do
    out=$B/draws/rep$r/draw${d}_demes_a0$a.npz
    [[ -e $out ]] && continue
    $CLI run-draw --tree $B/singer_out/sim$r/sim$r/trees/sim$r.$d.tsz --draw-id $d \
      --controls-vcf $B/sim/rep$r/controls.vcf.gz --controls-mode sample --ancient-vcf $B/sim/rep$r/ancient.vcf.gz \
      --panel-names $B/sim/rep$r/panel_names.txt --chrom-fai $B/chr1.fai --a0 $a --age-mode child \
      --demography $B/demography_bottleneck.yaml --ploidy 1 \
      --panel-fixed $B/draws/rep$r/panel_fixed.npz --fixed-polarity ref --output $out > $B/logs/drawdemes_rep${r}_d${d}_a0$a.log 2>&1
  done
done
for a in 10000 20000 50000; do
  m=$B/merged/rep${r}_a0$a.tsv
  [[ -e $m ]] || $CLI merge --draw-results $B/draws/rep$r/draw{75,93,112,131,149}_a0$a.npz --output $m > $B/logs/merge_rep${r}_a0$a.log 2>&1
done
