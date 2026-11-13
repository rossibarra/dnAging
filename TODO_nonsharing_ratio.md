# TODO: non-sharing-ratio estimator

Method: "Dating ancient genomes with inferred ARGs" (Nov 11 2026).
Core code, unit tests and an msprime driver exist; real-data inputs and inferred ARGs do not yet.

## Done

- [x] `pipeline/nonsharing_ratio.py`: site selection, f per target, H = log(f_anc / mean f_mod),
      hazard inversion, c(t) from lineage counts.
- [x] `tests/test_nonsharing_ratio.py`: unit tests (24 across the three ratio test files), passing on a compute node.
- [x] `validation/nonsharing_ratio_validation.py`: msprime driver, true-age or edge-midpoint site selection.
- [x] 20 Mb constant-Ne validation (20 replicates, true ARG, `results/nonsharing_ratio/const20mb_seed*.tsv`).
      Mean bias per age is within about 1-2 SE (SE roughly 60-130 generations) of zero for 250-6,000;
      at 8,000 the c(t)-from-ARG arms are -203 +/- 120 (edge ages) and -334 +/- 126 (true ages).
      20 replicates; the 8,000 result is borderline and is worth a larger run.
- [x] `pipeline/nonsharing_ratio_demography.py`: c(t) from ARG lineage counts (default) or a demes file
      (`rates_from_demes`: constant/exponential/linear epochs, ancestor chains, years -> generations).
      Tests in `tests/test_nonsharing_ratio_demography.py`.
- [x] `pipeline/nonsharing_ratio_io.py`: stdlib VCF reader (sample or phased columns, biallelic SNPs),
      polarisation onto ARG sites by allele string, missing calls handled (per-target denominators).
- [x] `pipeline/nonsharing_ratio_real_data.py`: `run-draw` / `merge` CLI (controls VCF, ancient VCF, one SINGER draw per
      call; per-block pair counts; block bootstrap with the same blocks across draws).
      End-to-end test against a direct computation in `tests/test_nonsharing_ratio_real_data.py`.

## Target interface (real data)

Inputs:
1. **Controls VCF**: held-out modern genomes (outside the ARG panel).
2. **Unknown ancient VCF**: the samples to date.
3. **ARG set**: SINGER posterior draws (tszip `.trees`, built on the panel only).
   Open question: reuse the ARGtest-concatenated layout that `pipeline/insertion_real_data.py` reads
   (`chrom_offsets` metadata or `--chrom-fai`), or accept per-chromosome files.

Also needed: panel haplotype identity (which samples built each ARG), c(t) source
(ARG lineage counts or a separate demographic estimate), and cutoff A0.

Output: per-sample T_hat per draw and merged across draws, with block-bootstrap intervals.

## Work items

### Inputs and plumbing

- [~] **Panel-fixed sites missing from the ARG's site list.** (optional `prepare-fixed`/`--panel-fixed`/`--polarity-table` path built and tested; real data still lacks an invariant-site panel table and an ancestral-allele source, so real-data rule 2 is met only for ARG sites that are all-derived. Simulation: excluding fixed sites gives -559 at T=4000, including gives +3.)
      Background: without it the CLI only sees panel-fixed sites that are ARG sites, but design rule 2 needs all of them.
- [x] Degenerate targets (f=0, H<=0, no_data) are flagged in stdout, NPZ and the merged TSV `flag` column.
- [x] Position convention verified on real data: shift 0 (ancient VCF 71.5% of records matched at 0 vs ~2% at +/-1; panel VCF 0 allele mismatches). Open: 4.5% allele mismatches in the ancient VCF are unexplained.
- [ ] Optionally accept the existing `chr*.npz` call format as well as VCFs.
- [~] Overlap check exists (controls vs ancient, and vs `--panel-names`); the panel names are not read from the ARG itself yet.
- [~] Polarisation uses the ARG site's ancestral/derived strings; check it against the existing polarity handling on real data.
- [~] Missing calls handled; ancient diploid hets become missing. Pseudo-haploid ancient data still needs a test.

### Estimator with inferred ARGs
- [~] Dated site age = edge midpoint of the ARG's mutation edge (`edge_midpoint_ages`); confirm POLEGON re-dated ARGs flow through unchanged.
- [~] Current merge: per-draw T_hat, averaged over draws. Compare with pooling counts or averaging site ages on simulation.
- [x] c(t) from the draw's κ(t) or a demes file. [ ] Parametric fit option for sharp demographic change (report limitation).
- [ ] Pick A0 from the data (above any plausible sample age); report sensitivity to A0.
- [x] Block-bootstrap intervals (blocks resampled identically across draws).

### Validation
- [ ] Larger constant-Ne run (more replicates) to settle the 8,000-generation result.
- [x] Run on SINGER + POLEGON ARGs from msprime (report's main claim). `batches/nonsharing_ratio_realsims_20261112/`: 5 reps, 20 Mb, 26 panel + 10 controls + 10 ancient (50 total), constant Ne, child-node A0. ARG: r 0.948 / 0.904 / 0.917 and mean bias -714 / -611 / +4,025 at A0 10k / 20k / 50k. A0 50k is biased high on the true tree too (cause untested).
- [x] Bottleneck demography: `batches/nonsharing_ratio_realsims_bottleneck_20261112/` (Ne x0.3 at 500 gens ago, exponential recovery by 300). ARG: r 0.937 / 0.927 / 0.936, bias -438 / -352 / +5,370. True tree is more biased high than under constant size (+655 / +880 at A0 10k / 20k), cause untested. 5 replicates, so differences from constant size are not resolved. Variable-Ne result is for this one demography only.
- [ ] Panel-as-controls vs held-out controls (report's design rule 1) and with/without panel-fixed sites (rule 2).
- [ ] Ascertainment: simulate a joint call set that only has sites variable across its samples.
- [ ] Genotype error / damage / low coverage at the selected sites; calibrate on old alleles carried by every modern.
- [ ] Population structure check: ratio should be flat across mutation ages. Direct test launched 2026-11-12: bottleneck simulation with the ancient samples from a population that split 6,000 generations ago (`batches/nonsharing_ratio_realsims_split6k_bottleneck_20261112/`; read its README for the result). Only one split time; other split times and controls drawn from the ancient samples' population are not run.
- [~] Radiocarbon-dated samples (`radiocarbon_goodvcfs/`, `batches/radiocarbon_*`): first real-data run done, see "Real-data results" below.

### Real-data results (2026-11-12; see `docs/nonsharing_ratio/REAL_DATA_OFFSET_REPORT.md`)
- [x] `--age-mode child` added (A0 on the child-node time of the mutation edge). The default edge midpoint kept 88-92% of sites because tip-ward edges with deep parents have huge midpoints.
- [x] Held-out stand-in (`batches/nonsharing_ratio_heldout_20261111/`, chr10, draw 75, 20-haplotype panel, 6 controls, 46 dated ancient samples): r between T_hat and calibrated radiocarbon age is -0.03 to 0.07 in every setting; with child-node A0 20,000 every sample is near 10,000 generations (ancient f 0.418 vs controls 0.318, ratio 1.31).
- [x] Read depth does not explain the offset (median f by depth tercile 0.4165 / 0.4189 / 0.4190; Spearman +0.26, wrong sign). Simple reference bias does not (excess positive for both derived=ALT, +0.060, and derived=REF, +0.079). Panel-fixed sites (REF polarity or full-ARG polarity table) change nothing.
- [ ] Offset cause is open. Leading hypothesis (untested): ancient samples are more diverged from the panel/controls than controls are from the panel. Predictions are in the report.
- [x] More chromosomes? `batches/nonsharing_ratio_heldout_20261111/block_check/` (chr10 draw 75 in 10 chunks of 15-20 Mb, n=46): r vs age per chunk ranges -0.31 to +0.45 (mean -0.03); pooling more chunks shrinks the noise in per-sample f about 1/sqrt(k) but r stays near 0 (k=1 -0.02, k=5 +0.03, all of chr10 +0.07) and the offset stays about 1.3 in every chunk. Inference (one draw, one chromosome): more chromosomes of the same data probably will not help, though per-chunk SD of r (0.24) cannot rule out a true r of 0.2-0.3. Per-sample f is poorly reproducible between disjoint halves of chr10 (mean correlation 0.13).
- [x] More draws? `batches/nonsharing_ratio_heldout_20261111/more_draws_check/` (10 draws, 75-138, same child-node A0 20,000 settings, n=46 dated): per-sample f and T_hat are stable across draws (ICC 0.95 / 0.91; between-draw Pearson about 0.96), the offset is stable (ancient median f 0.4189 +/- 0.0013 vs controls 0.3176 +/- 0.0014; every ancient sample above the highest control in every draw), and r with calibrated age does not improve by averaging (T_hat averaged over 10 draws 0.043; single draws 0.004 to 0.072; flat for 1, 2, 5, 10 draws). Inference (not tested): more draws probably will not help, since the variation between samples is reproducible but unrelated to age and the 10 draws come from one SINGER chain. Only A0 20,000 with child-node ages was run across draws.
- [ ] Reconcile the 4 of 50 trusted samples with no calibrated date (3 uncalibrated, 1 waiting for a date). Sample-to-date mapping verified by name on 2026-11-12 (0 mismatches over the 50 trusted50 VCF columns; `private/check_vcf_order.py`).
- [ ] Insertion comparison: per-sample insertion ages for the same samples were not located, so r cannot yet be compared directly.
- [ ] Optional: controls as pseudo-targets on real data; A0 sweep on real data.
- [ ] Cleanup of large intermediates once no longer needed: `singer_out/chunks/` (about 82 GB) in the held-out batch; the first sims batch's `singer_out/`.
- Test 4 (asymmetry tests) stays on hold.

### New open items from the real-data wiring (2026-11-11)
- [ ] Real draws found are raw SINGER (75 draws, 26 haploid samples, whole genome, no `chrom_offsets`; needs `--chrom-fai`). No real POLEGON-dated draws found in any run record. Decide raw SINGER vs POLEGON-dated and which time scale.
- [ ] Choose A0 for real data (smoke test with A0=20,000 selected ~90% of sites; placeholder).
- [ ] Export a genome-wide numeric-label ancient VCF (only chr1 exists) from the PLINK files.
- [ ] Supply the held-out controls VCF (the stand-in builder is running as `batches/nonsharing_ratio_heldout_20261111/`).
- [ ] Cache the VCF reads across draws (each draw re-reads both VCFs; genome-wide per-draw time unmeasured).
- [ ] Time one full-genome draw (chr1 smoke: ~6-11 min, 21 GB; full-genome extrapolation unverified).

### Housekeeping
- [x] Slurm wrappers `slurm/nonsharing_ratio_{prepare_fixed,run_draws,merge}.sbatch` (lint clean; only run_draws run for real). [x] Entry in `PROJECT_MAP.md` (2026-11-12).
- [ ] Link from `tests.md` (asymmetry tests) once Tests 1 and 3 report.
