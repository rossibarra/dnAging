# Non-sharing-ratio estimator: ancient samples from a SPLIT population (B), panel and controls from A (2026-11-12/13)

Same design as `batches/nonsharing_ratio_realsims_bottleneck_20261112/` (see its README for the design) with one change: structure.
Ancestral population (Ne 20,000 diploid-equivalent = 40,000 haploid) splits **6,000 generations ago** into A and B, no migration.
Panel (26 haplotypes, the ARG) and the 10 held-out controls are sampled at time 0 from **A**; the 10 ancient samples per replicate
are sampled from **B** at their true ages (stratified-uniform, one per 600-generation bin over 0-6,000, same seeds as the other batches).
Both A and B are independently reduced by 70% (to 0.3 Ne) 500 generations ago, with exponential recovery to Ne by 300 generations ago
(`demography_split6k_bottleneck.yaml`). The follow-up batch with controls drawn from B is
`batches/nonsharing_ratio_realsims_split6k_Bcontrols_20261113/`.

Limits: 5 replicates; one split time; no migration; no read-level noise, masks or damage. All numbers below are measurements on those 5 replicates.

## Results (50 ancient samples; bias = inferred - true, generations; mean of 5 draws per replicate)
Files: `summary_by_a0.tsv`, `per_sample_results.tsv`, `offset_by_rep.tsv`, `offset_summary.tsv`, `site_selection.tsv`, `controls_pseudotargets.tsv`,
`inferred_vs_true.png/.svg` (ARG c(t)) and `inferred_vs_true_truedemography_ct.png/.svg` (c(t) from the true deme-A demography). `*.tsv` is gitignored; the tables here carry the numbers.

| source | A0 | mean bias | median bias | RMSE | r | slope | 0-1500 | 1500-3000 | 3000-4500 | 4500-6000 | flagged |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ARG, c(t) from ARG | 10,000 | +1569 | +1743 | 2291 | 0.181 | 0.042 | +3603 | +2331 | +782 | -811 | 0 |
| ARG, c(t) from ARG | 20,000 | +1914 | +2065 | 2543 | 0.214 | 0.076 | +3904 | +2637 | +1102 | -324 | 0 |
| ARG, c(t) from ARG | 50,000 | +13877 | +14249 | 14069 | 0.205 | 0.238 | +15570 | +14694 | +12935 | +12110 | 0 |
| ARG, c(t) true demography | 10,000 | +1927 | +2111 | 2550 | 0.184 | 0.044 | +3956 | +2694 | +1144 | -455 | 0 |
| ARG, c(t) true demography | 20,000 | +2298 | +2429 | 2849 | 0.210 | 0.080 | +4283 | +3018 | +1489 | +68 | 0 |
| ARG, c(t) true demography | 50,000 | +15243 | +15415 | 15471 | 0.196 | 0.276 | +16837 | +16080 | +14327 | +13538 | 0 |
| true tree, c(t) from tree | 10,000 | +3810 | +3931 | 4169 | 0.149 | 0.043 | +5940 | +4481 | +2998 | +1453 | - |
| true tree, c(t) from tree | 20,000 | +4451 | +4400 | 4770 | 0.202 | 0.092 | +6455 | +5074 | +3727 | +2185 | - |
| true tree, c(t) from tree | 50,000 | +17472 | +17533 | 17652 | 0.162 | 0.206 | +19070 | +18452 | +16720 | +15347 | - |
| true tree, c(t) true demography | 10,000 | +4077 | +4316 | 4411 | 0.169 | 0.047 | +6185 | +4755 | +3285 | +1712 | - |
| true tree, c(t) true demography | 20,000 | +4742 | +4629 | 5040 | 0.210 | 0.095 | +6712 | +5388 | +4044 | +2454 | - |
| true tree, c(t) true demography | 50,000 | +19135 | +19493 | 19357 | 0.143 | 0.223 | +20643 | +20182 | +18423 | +16984 | - |

Compare the same estimator without structure (bottleneck batch): ARG c(t) r 0.937 / 0.927 / 0.936, bias -438 / -352 / +5370 at A0 10k / 20k / 50k.
Per-replicate r (ARG c(t), A0 10k / 20k / 50k): rep1 0.64 / 0.57 / 0.77, rep2 -0.14 / -0.17 / -0.23, rep3 0.08 / 0.42 / 0.62, rep4 0.10 / -0.04 / 0.23, rep5 0.11 / 0.18 / -0.04 (pooled r above is therefore driven partly by between-replicate differences; replicates are few).

### Offset: median ancient f vs mean control f (f = fraction of selected (site, reference) pairs not shared; `offset_summary.tsv`)
| A0 | median ancient f | mean control f | ratio, ARG (mean over 5 reps; range) | ratio, true tree |
|---|---|---|---|---|
| 10,000 | 0.389 | 0.343 | 1.135 (1.132-1.139) | 1.199 |
| 20,000 | 0.323 | 0.282 | 1.146 (1.136-1.164) | 1.218 |
| 50,000 | 0.251 | 0.158 | 1.588 (1.490-1.663) | 1.752 |
(ratios are per replicate, f from the draw-summed counts, averaged over draws, then over replicates; the ancient/control f columns are averages of the per-replicate values.)
Measured: with controls from A and ancient from B, ancient samples show a roughly constant excess non-sharing of about 1.13-1.15 at A0 10k/20k, similar in all 5 replicates. The real-data value in the report is about 1.31 at A0 = 20,000 (different data; not directly comparable).

### Controls as pseudo-targets (A controls, leave-one-out vs other 9; `controls_pseudotargets.tsv`; 250 control x draw values per A0 for the ARG, 50 for the true tree)
| source | A0 | mean H | mean T_hat | median T_hat |
|---|---|---|---|---|
| ARG | 10,000 | -0.0001 | 217 | 0 |
| ARG | 20,000 | -0.0002 | 269 | 30 |
| ARG | 50,000 | -0.0011 | 689 | 88 |
| true tree | 10,000 | -0.0001 | 248 | 0 |
| true tree | 20,000 | -0.0002 | 297 | 129 |
| true tree | 50,000 | -0.0019 | 945 | 23 |
Controls from A date as modern (H ~ 0, positive T_hat means come from truncation at 0), so the A-controls are consistent with each other and with the panel while the B ancient samples are not.

## Interpretation (inference, not tested beyond this batch)
Measured: the split removes age recovery (r ~0.2 versus ~0.93 without structure) and adds a large positive bias, and the offset ratio is ~1.14 and flat across replicates.
My reading: population structure between the ancient samples and the panel/controls can produce a constant, age-independent offset of the kind seen on real data, and when the controls come from the panel population the estimator reads that offset as age. Whether controls from B remove it is tested in the Bcontrols batch.
Observed on the true tree as well as the ARG, so the offset is not an ARG-inference artefact in this simulation.

## Commands and jobs (repo root; hpc_run = ~/.claude/bin/hpc_run)
See `COMMANDS.md`. Jobs: simulation hpc_run 40619415; SINGER+POLEGON 40619425-29 (reps 1-5; 18:26-37:57 elapsed, all COMPLETED); prepare-fixed + true-tree hpc_run 40619430; per-draw runs `04_run_draws.sbatch` 40627391-95 (reps 1-5, 4:16-4:31 each, COMPLETED); analysis hpc_run 40627525 (`05_analyze.py` then `06_offset.py`).
Setup/SINGER stage (jobs up to 40619430) was run by an earlier agent; I did not re-verify its logs beyond reading the job states above and confirming 25 `.tsz` files. SINGER's "ballpark Ne" for rep 1 was 19,295 (log line seen); not checked for other reps.

## Notes
* `merged/` (from `04_run_draws.sbatch`) globs `draw*_a0*.npz`, which also matches the `_demes` NPZs, so those merged files mix ARG-c(t) and demes-c(t) draws. They are not used by `05_analyze.py`/`06_offset.py`. The Bcontrols batch uses explicit draw lists.
* Large intermediates (`singer_out/`, `sim/` VCFs, `draws/` NPZ) are gitignored (`*.tsz`, `*.vcf.gz`, `*.npz`, `*.tsv`); nothing was committed.
