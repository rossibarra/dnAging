#!/bin/bash -l
# run-draw for every draw x A0 x CONTROL ARM of ONE replicate (sbatch --export=ALL,REP=n).
# Arm A: controls.vcf.gz (controls from A, same pop as panel) -> draw{d}_a0{A0}.npz ; arm B: controlsB.vcf.gz (controls from B,
# same pop as ancient) -> drawB{d}_a0{A0}.npz.  Demes-c(t) variants add _demes.  The two arms differ ONLY in --controls-vcf.
# One writer per output file; each rep is its own job.  Needs the 5 POLEGON trees of that rep.
#SBATCH --job-name=nsrRun
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --output=logs/nsrRun-%j.out
#SBATCH --error=logs/nsrRun-%j.err
set -euo pipefail
source /etc/profile.d/modules.sh
module load conda
conda activate dnaging
cd /quobyte/jrigrp/beil/logan_collab/dnAging
B=batches/nonsharing_ratio_realsims_split6k_Bcontrols_20261113
CLI="python pipeline/nonsharing_ratio_real_data.py"
r=$REP
for arm in A B; do
  if [[ $arm == A ]]; then CV=controls.vcf.gz; X=""; else CV=controlsB.vcf.gz; X="B"; fi
  for d in 75 93 112 131 149; do
    for a in 10000 20000 50000; do
      out=$B/draws/rep$r/draw${X}${d}_a0$a.npz
      [[ -e $out ]] || $CLI run-draw --tree $B/singer_out/sim$r/sim$r/trees/sim$r.$d.tsz --draw-id $d \
        --controls-vcf $B/sim/rep$r/$CV --controls-mode sample --ancient-vcf $B/sim/rep$r/ancient.vcf.gz \
        --panel-names $B/sim/rep$r/panel_names.txt --chrom-fai $B/chr1.fai --a0 $a --age-mode child \
        --panel-fixed $B/draws/rep$r/panel_fixed.npz --fixed-polarity ref --output $out > $B/logs/draw${X}_rep${r}_d${d}_a0$a.log 2>&1
      out=$B/draws/rep$r/draw${X}${d}_demes_a0$a.npz
      [[ -e $out ]] || $CLI run-draw --tree $B/singer_out/sim$r/sim$r/trees/sim$r.$d.tsz --draw-id $d \
        --controls-vcf $B/sim/rep$r/$CV --controls-mode sample --ancient-vcf $B/sim/rep$r/ancient.vcf.gz \
        --panel-names $B/sim/rep$r/panel_names.txt --chrom-fai $B/chr1.fai --a0 $a --age-mode child \
        --demography $B/demography_split6k_bottleneck.yaml --deme A --ploidy 1 \
        --panel-fixed $B/draws/rep$r/panel_fixed.npz --fixed-polarity ref --output $out > $B/logs/draw${X}demes_rep${r}_d${d}_a0$a.log 2>&1
    done
  done
  for a in 10000 20000 50000; do
    m=$B/merged/rep${r}${X}_a0$a.tsv
    [[ -e $m ]] || $CLI merge --draw-results $B/draws/rep$r/draw${X}{75,93,112,131,149}_a0$a.npz --output $m > $B/logs/merge${X}_rep${r}_a0$a.log 2>&1
  done
done
