#!/bin/bash
# From the repo root, via hpc_run. (1) prepare-fixed on the PRE-dropping full panel VCF; (2) sanity check (a):
# the estimator on the TRUE panel tree sequence (true node times, age-mode child) for each replicate and A0.
set -euo pipefail
B=batches/nonsharing_ratio_realsims_20261112
CLI="python pipeline/nonsharing_ratio_real_data.py"
for r in ${REPS:-1 2 3 4 5}; do
  mkdir -p $B/draws/rep$r
  [[ -e $B/draws/rep$r/panel_fixed.npz ]] || $CLI prepare-fixed --panel-vcf $B/sim/rep$r/panel_full.vcf.gz --min-called 21 --output $B/draws/rep$r/panel_fixed.npz
  for a in 10000 20000 50000; do
    out=$B/draws/rep$r/true_a0$a.npz
    [[ -e $out ]] && continue
    $CLI run-draw --tree $B/sim/rep$r/true_panel.trees --draw-id 0 \
      --controls-vcf $B/sim/rep$r/controls.vcf.gz --controls-mode sample --ancient-vcf $B/sim/rep$r/ancient.vcf.gz \
      --panel-names $B/sim/rep$r/panel_names.txt --chrom-fai $B/chr1.fai --a0 $a --age-mode child \
      --panel-fixed $B/draws/rep$r/panel_fixed.npz --fixed-polarity ref --output $out > $B/logs/true_rep${r}_a0$a.log 2>&1
  done
done
