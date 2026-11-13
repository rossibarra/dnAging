# Non-sharing-ratio estimator on SINGER+POLEGON ARGs of msprime simulations with a BOTTLENECK (2026-11-12)

Repeat of `batches/nonsharing_ratio_realsims_20261112` (constant Ne) with ONE change: the simulated demography.
Everything else (seeds, sample design, SINGER/POLEGON settings, estimator settings, analysis) is copied from that batch; see its README for the full design.
The seeds are identical to the template, so the ancient ages (`sim/rep*/truth.tsv`) are the same as in the constant-Ne run (not diffed; same seed+code path, unverified).

## The demography (the only change)
`demography_bottleneck.yaml` (demes), single panmictic population, no structure; ancient samples are sampled at their true ages in that population.
* Ancestral Ne = 20,000 diploid-equivalent (msprime ploidy=1, so haploid size 40,000 in the file; pair coalescence rate 1/(2Ne) = 2.5e-5 as in the template).
* **"70%" was interpreted as reduced BY 70%**: the size drops instantaneously, 500 generations ago, to **30% of ancestral Ne** (0.3 Ne = 6,000).
* **Exponential recovery (forward in time)**: from 0.3 Ne at 500 generations ago to Ne at 300 generations ago (200 generations), constant Ne afterwards (300 -> 0) and before (500 -> older).
* Check (observed): `rates_from_demes` gives c(t)*N = 1 on 0-300 and >500, and 1.37 / 2.5 per haploid-N unit in the 300-400 / 400-500 bins (bin-averaged), i.e. the file is read as intended. msprime applied the same file via `Demography.from_demes`; I did not run an independent check of the msprime-simulated coalescence rates (the mean panel pairwise diversity of 9.78e-4 to 9.94e-4 over the 5 reps is consistent with Ne ~ 20k, as the bottleneck is short).

## Design (as in template)
5 reps, 20 Mb, 36 modern + 10 ancient haploids per rep (50 ancient in total, stratified-uniform ages over 0-6000, one per 600-gen bin per rep); ARG built on 26 panel haplotypes only; 10 modern held-out controls; ancient/controls not in the ARG; sites with derived allele only in ancient samples dropped from the ancient VCF; SINGER+POLEGON (150 samples, thin 10, draws 75/93/112/131/149); `run-draw --age-mode child --a0 {10000,20000,50000}`, `--controls-mode sample`, `--panel-fixed`, `--fixed-polarity ref`; mean of T_hat over 5 draws.
* **SINGER Ne**: not a config setting in this workflow (same as the template). The workflow estimates a "ballpark Ne" from the panel diversity. Observed for reps 1-5: 19,565 / 19,703 / 19,887 / 19,670 / 19,804 (template: ~19,990), i.e. effectively the ancestral Ne 20,000, not the bottleneck size. mu 1.25e-8, rho 1e-8.
* c(t) for the estimator: (default) from the ARG lineage counts; ALSO (new here) from the true demography via `--demography demography_bottleneck.yaml --ploidy 1` (rates_from_demes). Everything else in those runs is identical, so the difference isolates ARG-c(t) error.
* True-tree check: estimator on the true panel tree sequence, with c(t) from that tree and with c(t) from the demes file.

## Results (50 ancient samples; bias = inferred - true, generations)
Files: `summary_by_a0.tsv`, `per_sample_results.tsv`, `site_selection.tsv`, `controls_pseudotargets.tsv`, figures `inferred_vs_true.png/.svg` (ARG c(t), the 3-panel figure) and `inferred_vs_true_truedemography_ct.png/.svg` (ARG, true c(t)).

| source | A0 | mean bias | median bias | RMSE | r | slope | 0-1500 | 1500-3000 | 3000-4500 | 4500-6000 | flagged |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ARG, c(t) from ARG | 10,000 | -438 | -329 | 740 | 0.937 | 0.847 | -212 | -282 | -463 | -882 | 4 |
| ARG, c(t) from ARG | 20,000 | -352 | -397 | 730 | 0.927 | 0.895 | -157 | -331 | -385 | -579 | 4 |
| ARG, c(t) from ARG | 50,000 | +5370 | +4830 | 6573 | 0.936 | 2.943 | +983 | +4373 | +6961 | +9880 | 5 |
| ARG, c(t) true demography | 10,000 | -251 | -253 | 649 | 0.936 | 0.903 | -151 | -134 | -219 | -567 | 4 |
| ARG, c(t) true demography | 20,000 | -164 | -179 | 681 | 0.927 | 0.955 | -101 | -189 | -141 | -250 | 4 |
| ARG, c(t) true demography | 50,000 | +5835 | +5291 | 7119 | 0.935 | 3.101 | +1120 | +4743 | +7514 | +10758 | 5 |
| true tree, c(t) from tree | 10,000 | +655 | +622 | 1049 | 0.957 | 1.286 | -77 | +562 | +946 | +1282 | - |
| true tree, c(t) from tree | 20,000 | +880 | +774 | 1401 | 0.944 | 1.413 | -48 | +587 | +1223 | +1924 | - |
| true tree, c(t) from tree | 50,000 | +7461 | +7157 | 8916 | 0.926 | 3.500 | +1716 | +6523 | +9311 | +13282 | - |
| true tree, c(t) true demography | 10,000 | +741 | +813 | 1154 | 0.955 | 1.320 | -61 | +614 | +1060 | +1458 | - |
| true tree, c(t) true demography | 20,000 | +980 | +857 | 1530 | 0.941 | 1.455 | -34 | +653 | +1349 | +2135 | - |
| true tree, c(t) true demography | 50,000 | +7884 | +7367 | 9472 | 0.925 | 3.697 | +1783 | +6761 | +9814 | +14266 | - |

Flagged = targets flagged in any draw (mostly H<=0 at young ages); their truncated T_hat is included in all statistics, as in the template.

### Comparison with the constant-Ne batch (ARG, c(t) from ARG)
| A0 | constant Ne: r, bias | bottleneck: r, bias |
|---|---|---|
| 10,000 | 0.948, -714 | 0.937, -438 |
| 20,000 | 0.904, -611 | 0.927, -352 |
| 50,000 | 0.917, +4025 | 0.936, +5370 |
Constant-Ne true-tree rows (template): +198 / +371 / +5876 -> bottleneck true-tree: +655 / +880 / +7461.

Sites selected (mean over draws; `site_selection.tsv`): about 27.4k (36%) / 16.4k (22%) / 3.1k (4%) of ~75.5-76.1k usable ARG sites at A0 10k / 20k / 50k, essentially as in the template.

### Controls as pseudo-targets (`controls_pseudotargets.tsv`; leave-one-out vs other 9 controls; 250 control x draw values per A0 for ARG, 50 for true tree)
| source | A0 10k median / mean | 20k | 50k |
|---|---|---|---|
| ARG, c(t) from ARG | 54 / 246 | 33 / 306 | 0 / 783 |
| ARG, c(t) true demography | 62 / 272 | 37 / 334 | 0 / 841 |
| true tree, c(t) from tree | 139 / 287 | 185 / 399 | 469 / 1012 |
| true tree, c(t) true demography | 158 / 292 | 177 / 402 | 465 / 1022 |
Positive means are expected from truncation at 0. Controls are near zero for the ARG at A0 10k/20k; ~800 mean at 50k (template: ~900-1000).

## Observations and hypotheses (hypotheses are UNTESTED)
* Measured: with the bottleneck the ARG-based estimator still tracks true age (r 0.93-0.94) with a modest negative bias at A0 10k/20k (-440 / -350, slope 0.85-0.90) and a large positive bias at A0 50k (+5400, slope ~2.9) -- qualitatively the same as constant Ne; bias magnitudes differ (smaller negative bias, larger positive at 50k). With only 5 replicates of 10 samples and large per-sample scatter, I did not test whether the differences from the constant-Ne run are beyond noise.
* Measured: using the true demography instead of ARG lineage counts shifts the ARG bias upward by ~+190 (A0 10k) and ~+190 (20k) and ~+470 (50k) and slightly improves RMSE at 10k/20k (740->649, 730->681), worsens it at 50k. So ARG-c(t) error accounts for only a small part of the error here (about 190 of ~440 / ~350 at 10k/20k).
* Measured, suspicious: the TRUE tree gives a clearly positive bias (+655 / +880) that is much larger than in the constant-Ne run (+198 / +371), and it is the same with c(t) from the true demography (+741 / +980), so it is not due to c(t) estimation. Therefore the bottleneck seems to induce a bias in the estimator itself (older samples overestimated: +1282 / +1924 in the 4500-6000 bin). Cause not investigated; a guess is that the time-varying c(t) interacts with the child-age definition of A0 or the age-dependent selection of sites, but this is not tested. The ARG negative bias then is plausibly ARG-inference error offsetting this positive estimator bias (also a guess).
* A0 = 50k is badly biased even with the true tree (+7461), as in the template (constant Ne); this is estimator behaviour at ~3k sites, not ARG error.
* Single-simulation caveat: one demographic scenario, 5 replicates, 50 targets; ancient ages (0-6000) span the bottleneck window (300-500), so samples are on both sides of it.

## Resources (sacct, observed; all COMPLETED)
Simulation: hpc_run srun job (5 reps, one job). SINGER+POLEGON per rep, 12 CPU/24G: jobs 40598362-66, elapsed 18:17, 23:44, 25:19, 26:36, 19:24 (rep1..5; MaxRSS of rep 1 5.8 GB). run-draw (15 ARG-c(t) + 15 demes-c(t) runs + 3 merges per rep, 1 CPU): jobs 40598778 (rep1), 40598781 (rep2), 40598782 (rep3), 40598798 (rep4), 40598779 (rep5), 3:32-4:30 each. prepare-fixed + true-tree runs: hpc_run job (logs/03_*.out). Analysis: hpc_run job 40598805. Heavy work always through hpc_run/sbatch (account jrigrp, low partition).

## Cleanup / layout
After the 25 POLEGON `.tsz` draws (draws 75/93/112/131/149 x 5 reps) were present, `singer_out/sim*/sim*/{chunks,plots,sim*.vcf}` were deleted (singer_out 19 GB -> 338 MB). Kept: `singer_out/sim*/sim*/trees/*.tsz` + their logs and sim*.log. `sim/` VCFs + truth, `draws/rep*/` (draw NPZs: `draw{d}_a0{A0}.npz` ARG c(t); `draw{d}_demes_a0{A0}.npz` true c(t); `true_*`, `truedemes_*` true tree), `merged/` (merge output, ARG c(t) only), `logs/`, `scripts/`.
Nothing outside this directory was modified; nothing committed.
