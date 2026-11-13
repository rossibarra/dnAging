#!/bin/bash
# Verify the stand-in outputs with the non-sharing-ratio CLI (run via hpc_run from the repo root, e.g.
#   HPC_MEM=16G HPC_TIME=01:00:00 ~/.claude/bin/hpc_run 'bash batches/nonsharing_ratio_heldout_20261111/scripts/04_verify_run_draw.sh 75')
# The real ancient genotypes are read from their existing location; the chr10 export is written ONLY to
# node-local scratch ($TMPDIR) and removed on exit, so nothing ancient enters the repo.
set -euo pipefail
DRAW=${1:?draw id}
A0=${A0:-12000}
B=/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nonsharing_ratio_heldout_20261111
FT=/quobyte/jrigrp/beil/logan_collab/dnAging/radiocarbon_goodvcfs/filtered_trusted_20261012
PREFIX=$FT/ARG_site_ancient_haps_mind95.trusted50.mappable_q37
CHR=chr.10.combined.snp.te.sorted
# plink 1.9 by absolute path: `module load plink` inside an activated conda env fails (module runs `conda deactivate`)
PLINK=/cvmfs/hpc.ucdavis.edu/sw/spack/environments/main/view/generic/plink-1.9-beta7.1/bin/plink
W=$(mktemp -d "${TMPDIR:-/tmp}/nsr-verify.XXXXXX")
trap 'rm -rf "$W"' EXIT
# same export recipe as export_chr1_vcf.sh, chromosome 10 = NC_050105.1, relabelled to "10" like the panel
"$PLINK" --bfile "$PREFIX" --allow-extra-chr --chr NC_050105.1 --keep-allele-order --recode vcf --out "$W/anc" >/dev/null
awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$1="10"; print}' "$W/anc.vcf" | gzip -c > "$W/ancient_chr10.vcf.gz"
cd /quobyte/jrigrp/beil/logan_collab/dnAging
mkdir -p "$B/verify"
OUT="$B/verify/run_draw_${DRAW}${TAG:-}.npz"
python pipeline/nonsharing_ratio_real_data.py run-draw \
  --tree "$B/singer_out/$CHR/trees/$CHR.$DRAW.tsz" --draw-id "$DRAW" \
  --controls-vcf "$B/controls/heldout_controls_chr10.vcf.gz" \
  --ancient-vcf "$W/ancient_chr10.vcf.gz" \
  --panel-names "$B/inputs/panel_samples.txt" \
  --chrom-fai "$B/chr10.fai" --a0 "$A0" ${EXTRA:-} --output "$OUT" 2>&1 | sed -E 's/^[^\t]+\tT_hat/<ancient sample>\tT_hat/' 
