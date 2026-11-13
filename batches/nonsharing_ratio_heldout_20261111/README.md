# Held-out-controls stand-in for the non-sharing-ratio estimator (chr10, maize NAM+TE panel)

Purpose: test `pipeline/nonsharing_ratio_real_data.py` end to end before a genuinely new modern VCF exists.
Design rule 1 (controls must be outside the ARG) is met by re-inferring ARGs from the modern panel MINUS k
held-out lines and using those lines as controls. When the real controls VCF arrives, keep using the existing
(full-panel) ARG draws; this directory is only a test bed.

## What was built (all paths relative to this directory)

| Item | Path |
|---|---|
| Held-out controls VCF (6 lines, chr10, bgzip+tbi) | `controls/heldout_controls_chr10.vcf.gz` |
| Panel sample names used for the ARG (for `--panel-names`) | `panel_names.txt` (= `inputs/panel_samples.txt`, 20 names incl. `REF`) |
| Held-out names / selection record | `inputs/heldout_samples.txt`, `inputs/selection.json` |
| SINGER input dir (panel minus held-out) | `inputs/for_te_chr10_panel/` |
| 10 POLEGON-dated SINGER draws (tszip, chr10 only) | `singer_out/chr.10.combined.snp.te.sorted/trees/chr.10.combined.snp.te.sorted.{75,82,89,96,103,110,117,124,131,138}.tsz` (about 80 MB each) |
| FAI for the single-chromosome draws (`--chrom-fai`) | `chr10.fai` |
| SINGER+POLEGON config | `singer_config.yaml` |
| Scripts / sbatch / exact commands / manifest | `scripts/`, `COMMANDS.md`, `manifest.json` |
| Verification output (de-identified) | `verify/run_draw_75.deidentified.npz` |
| First attempt that failed (k=5), logs only of interest | `failed_k5_attempt/` (safe to delete) |

## Choices

* **Panel.** The real SINGER+POLEGON inputs are the maize NAM+TE per-chromosome files
  (`/quobyte/jrigrp/beil/te_evo/singer_analysis/singer-snakemake/for_te_chr10/`): 26 haploid columns (`REF`, the
  B73 reference column, plus 25 inbred lines), SNPs and TE presence/absence sites together, per-sample masks,
  accessibility mask, genetic map, omit list. The existing recipe is `singer-snakemake` (SINGER 0.1.8-beta custom,
  POLEGON 0.1.3-alpha) with the `chr.9.yaml` settings: 150 MCMC samples, thin 10, POLEGON 150/10/5 scaling reps.
  (The repo's own `batches/singer_polegon_pilot_20261102` is a simulation pilot, so it was not the real-data recipe.)
* **k = 6, seed 20261111.** `random.Random(20261111).sample(sorted(25 non-REF columns), 6)` ->
  CML322, CML69, Il14H, Ki11, M162W, Tzi8. `REF` stays in the panel (SINGER's input needs it). 6/25 = 24%, inside the
  suggested 4-8. k=5 (20%) was tried first and failed: the workflow converts haploid VCFs to diploid and
  asserts an even number of samples (21 columns). k=4 or 6 give an even panel; 6 was chosen to give more controls.
  Caveat (not tested): the draw has 3 tropical/CML-type lines among the 6 held out, so the controls are not a random sample of ancestry groups in any
  balanced sense.
* **Chromosome 10 only.** It is the smallest nuclear chromosome (152.4 Mb vs 308 Mb for chr1), 174 chunks of
  about 1 Mb. chr1 would have roughly doubled the cost; chr10 already cost about 60 CPU-hours.
* **Draws.** SINGER must run its whole chain (150 samples x thin 10) whatever the number of draws kept. Of the 75
  post-burn-in samples, 10 evenly spaced ones (75, 82, ..., 138) were POLEGON-dated and merged.
* **Ne.** The workflow estimates a constant "ballpark" Ne per chromosome from diversity (96,672 here vs 101,546 in the
  26-line real run); this is not overridable in the workflow and there is no time-varying trajectory in this recipe
  (the `Ne` passed to SINGER and POLEGON is that constant).
* **Not done:** the argtest post-processing (low-recombination / low-access / mutation-load masks, per-sample trimming,
  chromosome concatenation) that produced the real `run.combined.N.tsz`. These draws are the direct
  singer-snakemake merged trees. See "Differences from the real draws".

## Controls VCF format

Produced by `bcftools view -S heldout_samples.txt` on the same file that the panel was cut from, so it has exactly the
panel's conventions: chromosome label `10`, 1-based POS, REF/ALT strings from the reference genome (not polarised
ancestral/derived: the polarity is set per ARG site by SINGER and the CLI polarises by allele strings), haploid GT
(`0`, `1`, `2` for multiallelic, `.`), TE presence/absence records retained, no site filtering (3,628,984 records).
The panel VCF has 3,070,359 records: sites no longer segregating after removing the 6 lines were dropped
(`bcftools view -c 1:minor`), the REF column stays first.

## Verification (observed, not assumed)

Jobs: prep 40590168 (2 CPU, 8G, 1 min); SINGER+POLEGON 40590169 (low/jrigrp-low-qos, 48 CPU, 64G, 3 h 13 min wall,
TotalCPU 2-11:43:11, about 59.7 CPU-hours, state COMPLETED; all 10 trees written). Budget: far below 500 CPU-hours, 48 CPUs
at a time. SINGER finished all 174 chunks; the 10 merges took the last 40 minutes or so.

Tree properties (draw 75): 20 haploid samples, sequence length 152,435,372 (FAI length + 1), 2,892,186 sites, 3,014,579
mutations, 1,831,370 trees, root time up to 1.7e7 generations (the real combined draw also reaches 1.9e7), time unit generations,
no `chrom_offsets` metadata.

`run-draw` on draw 75 with the real ancient VCF (chr10 export written to node-local scratch only, deleted on exit; the
PLINK source is read from its existing location, nothing ancient was copied into the repo), controls = this VCF, `--a0 12000`,
`--panel-names panel_names.txt`: exit 0, 196.6 s. The overlap check passed, so no control/ancient sample is in the panel.
* Controls: 2,856,338 of 2,892,186 ARG sites matched, 0 allele mismatches (position shift 0 is right for the controls and thus for the panel coordinates).
* Ancient: 347,199 sites matched, 16,650 allele mismatches (4.6% of 363,849 positions found); not investigated, mismatch in REF/ALT strings
  against the ARG's ancestral/derived strings.
* 2,856,338 usable sites, 2,619,604 selected at a0 = 12,000, 0 panel-fixed sites added (TODO item: needs the polarity source).
* 50 ancient targets + 6 controls; ancient T_hat: min 0, median 594, max 1,399 generations; the CLI flagged 2 targets with H<=0 (T_hat = 0).
  I do not know the true ages (radiocarbon data was deliberately not read), so I make no accuracy claim.
* The saved NPZ was rewritten with ancient names replaced by `anc001...` in VCF column order (original removed).

## Problems / observations for the CLI owner (pipeline files untouched)

1. A single-chromosome draw has no `chrom_offsets` metadata, so `--chrom-fai` is required; `chr10.fai` with the chromosome named `10`
   worked (the padding=1 branch of `chromosome_layout` reproduced the sequence length).
2. At `--a0 12000` 92% of usable sites (2.62 M of 2.86 M) were selected. I did not check what fraction of edge-midpoint ages exceeds 12,000;
   if POLEGON ages are this old for most sites, the selection may not be doing much. Unverified.
3. `module load plink` fails inside an activated conda env (module runs `conda deactivate`); the verify script calls plink by absolute path.
4. The CLI was edited while I worked (the output gained a timing field and a flagged-target warning), so the verified behaviour is that of the
   version on disk at about 02:00 on 2026-11-12.

## Differences from the real draws (what a user should know before comparing)

* 20-column panel (19 lines + REF) vs 26; chr10 only; single draws not concatenated across chromosomes; ARG sample nodes are haploid
  (20) as in the real combined draw (26).
* No argtest trimming/masking step; constant Ne estimated from this panel.
* The ARGs depend on the 6 held-out lines being absent, so ages will differ from the real-data ARGs by an unknown amount.

## Housekeeping

`singer_out/.../chunks/` holds the per-chunk SINGER/POLEGON intermediates for all 150 reps (82 GB, about 100k files).
Only `trees/*.tsz` are needed downstream; the chunks directory can be deleted by the user (I did not, since a recursive delete was
blocked by the safety check and the files are cheap to keep until the tests are done).

## Later additions (2026-11-12)

Everything after the first verification run (draw 75, midpoint ages) was added on 2026-11-12. Outputs in `verify/` are
de-identified (S01.. ancient, C01.. controls); no dates or sample IDs are stored here.

| Item | Where | What it shows |
|---|---|---|
| A0 as child-node time (`--age-mode child`), A0 30,000 midpoint, A0 20,000 child | `verify/run_draw_75_{a30000,child_a20000}.npz`, `scripts/04_verify_run_draw.sh` (env vars A0, TAG, EXTRA) | midpoint kept 88-92% of sites; child-node A0 20,000 keeps 34% |
| Panel-fixed sites from the 20 panel columns only (held-out samples are not used) | `scripts/07_prepare_fixed_heldout.sh`, `panel_fixed_all_sites.npz`, `verify/run_draw_75_child_a20000_fixed.npz` | 389,112 fixed sites; with REF polarity 75 kept; no effect on r |
| Ancestral alleles for fixed sites from the full 26-haplotype real draw | `scripts/08_polarity_from_full_arg.py`, `polarity_chr10_from_full_arg.tsv`, `verify/run_draw_75_child_a20000_fixedtable.npz` | 27,874 fixed sites kept; r 0.06 |
| REF/ALT class split | `scripts/09_ref_alt_split.{sh,py}` | excess non-sharing is positive in both classes (+0.060 derived=ALT, +0.079 derived=REF) |
| More draws (10 draws, child-node A0 20,000) | `more_draws_check/` (own README) | offset and per-sample f stable; averaging does not raise r |
| chr10 in chunks | `block_check/` (own README) | r stays near 0 as sites are pooled; offset about 1.3 in every chunk |

Summary of conclusions and caveats: `docs/nonsharing_ratio/REAL_DATA_OFFSET_REPORT.md`.
Large and removable: `singer_out/.../chunks/` (about 82 GB); only the `trees/*.tsz` draws are needed.
