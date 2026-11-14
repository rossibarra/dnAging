#!/bin/bash
# From the repo root, via hpc_run. (1) prepare-fixed on the PRE-dropping full panel VCF; (2) sanity check (a):
# the estimator on the TRUE panel tree sequence (true node times, age-mode child) for each replicate, A0 and CONTROL ARM
# (armA = controls.vcf.gz from A, armB = controlsB.vcf.gz from B).  File prefix: true_ / trueB_ (c(t) from tree), truedemes_ / truedemesB_.
set -euo pipefail
B=batches/nonsharing_ratio_realsims_split6k_Bcontrols_20261113
CLI="python pipeline/nonsharing_ratio_real_data.py"
for r in ${REPS:-1 2 3 4 5}; do
  mkdir -p $B/draws/rep$r
  [[ -e $B/draws/rep$r/panel_fixed.npz ]] || $CLI prepare-fixed --panel-vcf $B/sim/rep$r/panel_full.vcf.gz --min-called 21 --output $B/draws/rep$r/panel_fixed.npz
  for arm in A B; do
    if [[ $arm == A ]]; then CV=controls.vcf.gz; X=""; else CV=controlsB.vcf.gz; X="B"; fi
    for a in 10000 20000 50000; do
      out=$B/draws/rep$r/true${X}_a0$a.npz
      [[ -e $out ]] || $CLI run-draw --tree $B/sim/rep$r/true_panel.trees --draw-id 0 \
        --controls-vcf $B/sim/rep$r/$CV --controls-mode sample --ancient-vcf $B/sim/rep$r/ancient.vcf.gz \
        --panel-names $B/sim/rep$r/panel_names.txt --chrom-fai $B/chr1.fai --a0 $a --age-mode child \
        --panel-fixed $B/draws/rep$r/panel_fixed.npz --fixed-polarity ref --output $out > $B/logs/true${X}_rep${r}_a0$a.log 2>&1
      out=$B/draws/rep$r/truedemes${X}_a0$a.npz
      [[ -e $out ]] || $CLI run-draw --tree $B/sim/rep$r/true_panel.trees --draw-id 0 \
        --controls-vcf $B/sim/rep$r/$CV --controls-mode sample --ancient-vcf $B/sim/rep$r/ancient.vcf.gz \
        --panel-names $B/sim/rep$r/panel_names.txt --chrom-fai $B/chr1.fai --a0 $a --age-mode child \
        --demography $B/demography_split6k_bottleneck.yaml --deme A --ploidy 1 \
        --panel-fixed $B/draws/rep$r/panel_fixed.npz --fixed-polarity ref --output $out > $B/logs/truedemes${X}_rep${r}_a0$a.log 2>&1
    done
  done
done
