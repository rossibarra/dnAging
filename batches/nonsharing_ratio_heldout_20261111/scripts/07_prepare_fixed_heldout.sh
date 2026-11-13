#!/bin/bash
# Panel-fixed sites for the held-out chr10 stand-in: all sites of the 20 panel columns (including the sites
# SINGER dropped as non-segregating), scanned for sites where every called panel haplotype carries one allele.
# Run from the repo root through hpc_run. Intermediate VCF is node-local and removed on exit.
set -euo pipefail
B=batches/nonsharing_ratio_heldout_20261111
SRC=/quobyte/jrigrp/beil/te_evo/singer_analysis/singer-snakemake/for_te_chr10/chr.10.combined.snp.te.sorted.vcf.gz
# bcftools by absolute path: 'module load' inside the hpc_run conda env swaps python
BCFTOOLS=/cvmfs/hpc.ucdavis.edu/sw/spack/environments/main/view/generic/bcftools-1.19/bin/bcftools
W=$(mktemp -d "${TMPDIR:-/tmp}/nsr-fixed.XXXXXX"); trap 'rm -rf "$W"' EXIT
$BCFTOOLS view -S $B/inputs/panel_samples.txt --force-samples "$SRC" -Oz -o "$W/panel_all.vcf.gz"
python pipeline/nonsharing_ratio_real_data.py prepare-fixed --panel-vcf "$W/panel_all.vcf.gz" \
  --min-called 17 --output $B/panel_fixed_all_sites.npz
