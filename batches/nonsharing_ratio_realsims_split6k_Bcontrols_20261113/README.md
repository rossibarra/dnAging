# Do controls from the ANCIENT population (B) remove the structure offset? (2026-11-13)

Question: when the ancient samples come from a population (B) diverged from the panel's (A), does using modern controls from B
(rather than A) cancel the constant ancient-vs-control offset, and does the estimator then recover ages?

Design = `batches/nonsharing_ratio_realsims_split6k_bottleneck_20261112/` (read its README) with one addition. Ancestral population (haploid 40,000)
splits 6,000 generations ago into A and B, no migration, both with a 0.3 Ne bottleneck 500 generations ago and exponential recovery by 300
(`demography_split6k_bottleneck.yaml`, copied unchanged). Per replicate (5 reps, 20 Mb, seeds 20261112 + 100*rep as in the other batches) msprime
samples, in this order: 26 panel (A, time 0; the ARG is built on these only), 10 controls-A (A, time 0), **10 controls-B (B, time 0)**, 10 ancient (B, ages
stratified-uniform over 0-6,000, same ages as the earlier batches). Because of the extra samples the trees differ from the earlier split batch, so these are new
simulations (not the same data). SINGER + POLEGON on the 26 panel haplotypes, draws 75/93/112/131/149, copied from the earlier batch.
Estimator: `run-draw --age-mode child --a0 {10000,20000,50000} --controls-mode sample --panel-fixed ... --fixed-polarity ref`, mean of 5 draws, c(t) from ARG
lineage counts (and, separately, from the true deme-A demography). **The two arms differ only in `--controls-vcf`**: arm A = `controls.vcf.gz` (A),
arm B = `controlsB.vcf.gz` (B). Same ARG draws, same ancient VCF, same panel-fixed file. The true-tree check was run with both control sets.
The ancient VCF keeps sites with a derived allele in the panel or in either control set (shared by both arms; the earlier batch used panel + A controls).

Limits: 5 replicates per condition; one split time (6,000, older than every ancient age); no migration; no masking/damage/read noise. B-controls here are drawn from
exactly the ancient population (best case); real moderns from the ancient samples' region would not be exact descendants. Nothing here is a statement about the real data.

## Results (50 ancient samples; bias = inferred - true, generations; ARG, c(t) from ARG)
Files: `summary_by_a0.tsv` (all sources incl. demes c(t) and true tree), `per_sample_results.tsv`, `offset_by_rep.tsv`, `offset_summary.tsv`,
`site_selection.tsv`, `controls_pseudotargets.tsv`, figures `inferred_vs_true.png/.svg` (top row arm A, bottom row arm B) and `inferred_vs_true_truedemography_ct.png/.svg`.
`*.tsv` is gitignored, so these tables carry the numbers.

| arm (controls) | A0 | mean bias | median bias | RMSE | r | slope | targets flagged (H<=0, any draw) | median ancient f | mean control f | offset ratio, ARG (range over reps) | ratio, true tree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 10,000 | +1905 | +1830 | 2586 | 0.065 | 0.021 | 0 | 0.395 | 0.346 | 1.144 (1.138-1.149) | 1.216 |
| A | 20,000 | +2281 | +2411 | 2884 | 0.108 | 0.045 | 0 | 0.329 | 0.285 | 1.155 (1.145-1.162) | 1.239 |
| A | 50,000 | +14094 | +13899 | 14290 | 0.083 | 0.088 | 0 | 0.260 | 0.163 | 1.596 (1.519-1.660) | 1.744 |
| B | 10,000 | -2716 | -2700 | 3198 | 0.111 | 0.018 | 29 | 0.392 | 0.391 | 1.004 (1.000-1.009) | 1.005 |
| B | 20,000 | -2576 | -2600 | 3085 | 0.122 | 0.031 | 29 | 0.325 | 0.323 | 1.007 (1.000-1.013) | 1.006 |
| B | 50,000 | -1632 | -1912 | 2603 | 0.071 | 0.052 | 23 | 0.249 | 0.241 | 1.034 (1.022-1.048) | 1.036 |

(offset ratio = per replicate median-ancient-f / mean-control-f, averaged over draws then replicates; f columns are replicate averages.) Bias by true-age bin (0-1500 / 1500-3000 / 3000-4500 / 4500-6000), ARG c(t):
A10k +3989/+2736/+925/-331; B10k -582/-1977/-3639/-4995; A20k +4359/+3000/+1369/+85; B20k -462/-1889/-3449/-4839. True tree: A r 0.16/0.18/0.09 and B r 0.13/0.09/0.13 at 10k/20k/50k, B bias -2706/-2538/-1353.
Per-replicate r (ARG, 10k/20k/50k): arm A rep1 0.65/0.56/0.25, rep2 0.22/0.26/0.31, rep3 0.03/0.16/0.10, rep4 -0.24/-0.26/-0.05, rep5 -0.24/-0.04/-0.04; arm B rep1 0.71/0.61/0.29, rep2 0.07/0.12/0.31, rep3 0.04/0.12/0.00, rep4 -0.19/-0.24/-0.16, rep5 -0.11/0.04/-0.07.
Reference without structure (bottleneck batch, ARG c(t)): r 0.937/0.927/0.936, bias -438/-352/+5370.

### Controls as pseudo-targets (`controls_pseudotargets.tsv`; 250 control x draw values per A0 for ARG, 50 for the true tree)
A-controls vs the other 9 A-controls: mean H -0.0001/-0.0002/-0.0014 (A0 10k/20k/50k), median T_hat 0. B-controls vs the other 9 B-controls: the same, ~0
(mean H -0.0001/-0.0001/-0.0008, median T_hat 0-2). This within-set test cannot show whether B-controls are "modern relative to the ancient samples", so I also computed
**cross pseudo-targets** (control f values from the arm-B and arm-A runs of the same draw and A0, separate run-draw calls over the same ARG):

| pseudo-target vs reference | A0 | mean H (ARG) | median T_hat (ARG) | mean H / median T_hat (true tree) |
|---|---|---|---|---|
| B-controls vs mean of A-controls | 10,000 | +0.122 | 4,439 | +0.182 / 6,869 |
| B-controls vs mean of A-controls | 20,000 | +0.126 | 4,533 | +0.196 / 7,318 |
| B-controls vs mean of A-controls | 50,000 | +0.388 | 14,273 | +0.471 / 17,443 |
| A-controls vs mean of B-controls | 10,000 | -0.123 | 0 (truncated) | -0.183 |

So modern B samples read as "4,400 generations old" against A-controls in the ARG (about 6,900 on the true tree) although their true age is 0; this has the same size as the offset of the B ancient samples
(ratio 1.14, H ~0.13-0.14 in arm A).

## What was measured and what I infer (inference is not tested here)
* Measured: controls from A leave a constant excess non-sharing of the ancient samples (median ancient f / mean control f about 1.14-1.16 at A0 10k/20k; the true tree gives 1.22-1.24); ages are not recovered (r 0.07-0.11) and biased high (+1.9k to +2.3k).
* Measured: controls from B remove the offset (ratio 1.004-1.007 at A0 10k/20k, 1.03 at 50k).
* Measured: with B-controls ages are still not recovered (r 0.07-0.12, slope 0.02-0.05), and the estimate collapses toward zero (bias -2.7k at 10k/20k; 29 of 50 targets have H<=0 in at least one draw and are truncated to 0).
* Inference (untested): in this design every ancient age (<6,000) is younger than the split, and the old alleles the estimator uses (above A0 >= 10,000 generations) are older than the split, so the age of a B sample may change its sharing with the A panel by little compared with the structure term; with matched controls nothing else is left, so the age signal is not detectable at 5 draws x 5 replicates. I did not test this explanation (for instance with ancient ages older than the split, or with controls only partly from B).
* Inference: matched controls are necessary to cancel a constant structure offset but not sufficient to recover ages when the panel itself is from the diverged population; whether any intermediate design (migration, a younger split, B-controls plus a B panel) recovers ages was not tried.
* Arm A here (r 0.07-0.11) differs from the earlier split batch (r 0.18-0.21): different simulations (new sample order, extra samples, ancient VCF rule) and only 5 replicates with large replicate-to-replicate variation (per-rep r from -0.26 to +0.65). I did not test whether the two differ beyond noise.

## Commands and jobs
See `COMMANDS.md`. Simulation hpc_run (5 reps; first run died at the end on a script bug writing `sim_info.json`; fixed; VCFs were verified identical on a re-run in scratch, and `sim_info.json` copied from it); SINGER+POLEGON 40627402-06 (reps 1-5; 18:44-28:59 elapsed, COMPLETED);
prepare-fixed + true-tree runs (both arms) hpc_run; per-draw runs (both arms, 30 ARG-c(t) + 30 demes-c(t) estimator runs + 6 merges per rep) 40627541-45 (4:42-9:48 elapsed, COMPLETED); analysis hpc_run 40627663 (an earlier attempt 40627662 failed on a missing `ARMS` line, fixed).
Nothing outside this directory was changed except documentation; nothing committed. `merged/` holds merge output per arm (`rep{r}_a0*.tsv` arm A, `rep{r}B_a0*.tsv` arm B), ARG c(t) only (explicit draw lists).
