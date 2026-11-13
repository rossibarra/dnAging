#!/bin/bash
# run via: A0=20000 ~/.claude/bin/hpc_run 'bash batches/.../scripts/09_ref_alt_split.sh 75'  (ancient chr10 export stays node-local)
set -euo pipefail
DRAW=${1:?draw}; A0=${A0:-20000}
B=/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nonsharing_ratio_heldout_20261111
PREFIX=/quobyte/jrigrp/beil/logan_collab/dnAging/radiocarbon_goodvcfs/filtered_trusted_20261012/ARG_site_ancient_haps_mind95.trusted50.mappable_q37
PLINK=/cvmfs/hpc.ucdavis.edu/sw/spack/environments/main/view/generic/plink-1.9-beta7.1/bin/plink
W=$(mktemp -d "${TMPDIR:-/tmp}/nsr-split.XXXXXX"); trap 'rm -rf "$W"' EXIT
"$PLINK" --bfile "$PREFIX" --allow-extra-chr --chr NC_050105.1 --keep-allele-order --recode vcf --out "$W/anc" >/dev/null
awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$1="10"; print}' "$W/anc.vcf" | gzip -c > "$W/ancient_chr10.vcf.gz"
cd /quobyte/jrigrp/beil/logan_collab/dnAging
python -I $B/scripts/09_ref_alt_split.py "$B/singer_out/chr.10.combined.snp.te.sorted/trees/chr.10.combined.snp.te.sorted.$DRAW.tsz" \
  "$B/controls/heldout_controls_chr10.vcf.gz" "$W/ancient_chr10.vcf.gz" "$B/inputs/panel_samples.txt" "$B/chr10.fai" "$A0"
