#!/bin/bash -l
# Build the SINGER+POLEGON input dir for "chr10 panel minus held-out lines" and the controls VCF.
#SBATCH --job-name=nsrPrep
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/nsrPrep-%j.out
#SBATCH --error=logs/nsrPrep-%j.err
set -euo pipefail
source /etc/profile.d/modules.sh
module load bcftools/1.19
B=/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nonsharing_ratio_heldout_20261111
SRC=/quobyte/jrigrp/beil/te_evo/singer_analysis/singer-snakemake/for_te_chr10
STEM=chr.10.combined.snp.te.sorted
IN=$B/inputs/for_te_chr10_panel          # SINGER input dir (same file layout as for_te_chr10)
CTRL=$B/controls
mkdir -p "$IN" "$CTRL"
[[ ! -e "$IN/$STEM.vcf.gz" ]] || { echo "refusing to overwrite $IN/$STEM.vcf.gz"; exit 3; }
cd "$B"
# Unchanged sidecars (genetic map, accessibility mask, omit list): plain copies.
cp "$SRC/$STEM.hapmap" "$SRC/$STEM.mask.bed" "$SRC/$STEM.omit.txt" "$IN/"
# Per-sample masks: drop held-out lines (column 4 = sample name).
awk 'NR==FNR{h[$1]=1;next} !($4 in h)' inputs/heldout_samples.txt "$SRC/$STEM.sample.mask.bed" > "$IN/$STEM.sample.mask.bed"
# Panel VCF: drop held-out columns, then drop sites that are no longer segregating
# (bcftools recomputes AC/AN after -s). REF column stays first.
bcftools view -S inputs/panel_samples.txt --force-samples "$SRC/$STEM.vcf.gz" -Ou \
  | bcftools view -c 1:minor -Ou | bcftools view -Oz -o "$IN/$STEM.vcf.gz"
# Controls VCF: held-out lines only, same record set/REF/ALT/coordinates/GT encoding as the panel VCF
# (haploid GT, chromosome label "10"), no site filtering.
bcftools view -S inputs/heldout_samples.txt --force-samples "$SRC/$STEM.vcf.gz" -Oz -o "$CTRL/heldout_controls_chr10.vcf.gz"
bcftools index -t "$CTRL/heldout_controls_chr10.vcf.gz"
{
 echo "panel_records $(bcftools view -H $IN/$STEM.vcf.gz | wc -l)"
 echo "source_records $(bcftools view -H $SRC/$STEM.vcf.gz | wc -l)"
 echo "controls_records $(bcftools view -H $CTRL/heldout_controls_chr10.vcf.gz | wc -l)"
 echo "panel_cols $(bcftools query -l $IN/$STEM.vcf.gz | wc -l)"
 echo "controls_cols $(bcftools query -l $CTRL/heldout_controls_chr10.vcf.gz | wc -l)"
 echo "sample_mask_names_remaining $(cut -f4 $IN/$STEM.sample.mask.bed | sort -u | wc -l)"
} | tee "$B/inputs/prep_summary.txt"
