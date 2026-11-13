# Does averaging over SINGER+POLEGON draws rescue the real-data non-sharing ratio? (chr10, held-out controls)

Answer (inference from the numbers below): **No, probably not.** Per-sample f and T_hat are very stable across
the 10 draws (ICC about 0.9-0.95), so averaging removes little variance, and the averaged T_hat is no better
correlated with calibrated age (Pearson 0.043, n = 46) than single draws (0.004-0.072). The ancient-vs-control
offset in f is the same in every draw.

## What was run (observed)
* Draws present: 10 files `singer_out/chr.10.combined.snp.te.sorted/trees/*.{75,82,...,138}.tsz`, all with a matching
  POLEGON log. Per the batch README they are the 10 POLEGON-dated draws (I did not re-inspect the trees).
* Settings identical to `verify/run_draw_75_child_a20000.npz` (meta: a0 = 20000, age_mode child, controls = held-out VCF,
  ancient = trusted50 chr10 export in node-local scratch): `A0=20000 EXTRA="--age-mode child"` through
  `scripts/run_draw_deid.sh` (a copy of the batch's `04_verify_run_draw.sh` with output sent to node-local scratch and then
  relabelled S01-S50 / C01-C06 in VCF column order, the same convention as the existing de-identified NPZs).
* Draws 82-138 (9 jobs): `scripts/submit_all.sh`, each `hpc_run` with 2 CPU / 32G / 1.5 h, low partition, jobs
  40611546-40611554 all started on low; all 9 produced `out/draw_<id>.deid.npz`; each log shows "selected ~988k of ~2.855M
  sites" at about 124 s. I did not run sacct for memory use. Draw 75 is a copy of the existing
  `verify/run_draw_75_child_a20000.npz` (not rerun).
* Aggregation: `scripts/aggregate.py out dates.txt`. The dates file (46 dated of 50, from the scratchpad cmp.py list) lives only in the
  scratchpad, not here. f per target = sum(num)/sum(den) over blocks.

## Results (n = 50 ancient, 6 controls, 10 draws)

| Quantity | f (ancient) | T_hat (ancient) | f (controls) |
|---|---|---|---|
| mean within-sample SD across draws | 0.0021 | 210 gen | 0.0022 |
| SD across sample means | 0.0096 | 698 gen | 0.0079 |
| within/between ratio | 0.22 | 0.30 | 0.28 |
| ICC(1) | 0.950 | 0.911 | 0.925 |
| between-draw Pearson, mean (range) | 0.962 (0.931-0.981) | 0.960 (0.924-0.981) | |
| between-draw Spearman, mean (range) | 0.956 (0.918-0.977) | 0.956 (0.918-0.977) | |

Correlation of T_hat with calibrated age (46 dated):

| | Pearson | Spearman |
|---|---|---|
| single draws (10) | 0.004-0.072 (mean 0.043) | 0.069-0.112 |
| T_hat averaged over 10 draws | 0.043 (p = 0.77) | 0.101 (p = 0.50) |
| f averaged over 10 draws | 0.051 | 0.100 |
| random subsets of 1 / 2 / 5 draws (200 resamples, mean) | 0.044 / 0.045 / 0.043 | |

Mean T_hat over draws: median 10,219, range 8,242-11,484 (all 50 ancient; no T_hat = 0).

Offset stability (per draw, ancient median f vs control mean f): ancient 0.4189 +- 0.0013 (SD over draws),
controls 0.3176 +- 0.0014, difference 0.1012 +- 0.0020 (range 0.0974-0.1041). In every draw 100% of ancient samples have f above the
highest control. Per-control f over draws 0.308-0.326.

## Interpretation (inference, not tested)
* Draw-to-draw noise is small relative to between-sample spread, so there is little for averaging to remove. The ancient
  samples' f differs consistently among themselves across draws (so it is not MCMC noise), yet that variation is unrelated
  to age; averaging cannot fix a bias that every draw shares.
* The 10 draws come from one SINGER chain on the same panel and data (thinned samples, likely correlated); they do not
  sample ARG uncertainty independent of panel choice, ancient data or Ne. More draws from the same chain would plausibly
  behave the same (extrapolated from 10 draws; the r flattening across subset sizes 1, 2, 5, 10 supports this).
* The mean r of 0.043 is not significantly different from 0 at n = 46; I did not test whether a larger n changes this.

## Caveats
* Results are for chr10, 20-haplotype panel, child-node A0 = 20,000 only.
* The age pairing (VCF column order to dates) comes from the existing scratchpad cmp.py; reproduced the draw-75 r = 0.07 of
  the offset report, which supports it.
* Not run: other A0 or edge-midpoint settings across draws.
