# Block / pooling check: could more sites rescue the real-data estimator?

Input: `../verify/run_draw_75_child_a20000.npz` (draw 75, chr10, child-node A0=20000). It DOES keep per-block counts:
`num`,`den` are 31 blocks x 56 targets (50 ancient + 6 held-out controls), blocks = 5 Mb (`10:0`..`10:30`, `--block-size` default),
so no rerun was needed (step 3 not used). T_hat per subset is recomputed exactly like the pipeline
(`nonsharing_ratio.ages_from_counts` on the summed subset blocks, controls pooled for the baseline).
Chunks = 10 disjoint groups of 3-4 consecutive blocks (15-20 Mb each).

Command (head node, about 5 s, small arrays only; ages are a file OUTSIDE the repo, see below):

    module load conda; conda activate dnaging
    python -I block_check.py ../verify/run_draw_75_child_a20000.npz . --ages <ages.txt outside repo>

Without `--ages` the age-dependent columns are NaN. Full console output: `block_check_output.txt`.
This directory holds only aggregates: no dates, no sample IDs. The ages file lives in the head-node scratchpad only.

## 1. Per chunk (n = 46 dated ancient samples for r; 50 for f)

| chunk | Mb | control mean f | ancient median f | offset (anc/ctrl) | median T_hat | r (Pearson) | rho (Spearman) |
|---|---|---|---|---|---|---|---|
| 0 | 0-20 | 0.361 | 0.466 | 1.29 | 9353 | -0.31 | -0.22 |
| 1 | 20-35 | 0.339 | 0.386 | 1.14 | 4775 | -0.03 | -0.02 |
| 2 | 35-50 | 0.232 | 0.373 | 1.61 | 14579 | +0.45 | +0.43 |
| 3 | 50-65 | 0.348 | 0.383 | 1.10 | 2809 | -0.23 | -0.23 |
| 4 | 65-80 | 0.355 | 0.395 | 1.11 | 3324 | -0.02 | 0.00 |
| 5 | 80-95 | 0.323 | 0.428 | 1.33 | 10228 | -0.14 | -0.09 |
| 6 | 95-110 | 0.223 | 0.360 | 1.62 | 14709 | -0.10 | -0.18 |
| 7 | 110-125 | 0.305 | 0.374 | 1.23 | 7647 | +0.16 | +0.11 |
| 8 | 125-140 | 0.332 | 0.452 | 1.36 | 10973 | -0.26 | -0.15 |
| 9 | 140-155 | 0.417 | 0.459 | 1.10 | 2673 | +0.21 | +0.22 |

* r across the 10 chunks: mean -0.03, SD 0.24, range -0.31 to +0.45. Across the 31 single 5 Mb blocks: mean +0.01, SD 0.24, 15/31 positive.
  (For reference, with n = 46 and no true association the SD of r is about 1/sqrt(45) = 0.15; this formula is textbook, not computed here.)
  The observed SD of 0.24 is larger than that, I have not tested why (chunk-to-chunk heterogeneity is one candidate, unverified).
* Offset (ancient median f / control mean f) is positive in every chunk (1.10 to 1.62; mean 1.29, SD 0.20) and shows no trend along the chromosome
  (Spearman vs position -0.05, p = 0.88). It is not specifically elevated in any contiguous region; the largest offsets (chunks 2 and 6, about 1.6) coincide
  with the lowest control f (0.23, 0.22), i.e. the offset varies mostly because the control baseline varies by chunk, and the ancient median varies less (0.36-0.47).
  Whether chunks 2 and 6 are pericentromeric/low-recombination I did not check.
* Per-sample f is very noisy per chunk: mean over samples of SD across chunks 0.057 on mean f 0.408 (CV 0.14), versus a between-sample SD of the all-chunk f of 0.0096
  (CV 0.024). Part of the across-chunk SD is a chunk effect shared by all samples (not removed), so it overstates sample-specific noise.
* Split-half reliability of per-sample f (random 5 vs 5 chunks, 200 splits, 50 ancient samples): mean r = 0.13, SD 0.13. Differences between ancient samples in f
  are mostly not reproducible between disjoint halves of the chromosome. (Inference: Spearman-Brown for full length would be about 0.23, still weak.)

## 2. Pooling curve (random k chunks out of 10; 500 draws for k<10)

| k chunks (about Mb) | mean r vs age | SD of r over draws | fraction r > 0.3 | mean rho | between-sample CV of f | mean offset |
|---|---|---|---|---|---|---|
| 1 (15) | -0.02 | 0.22 | 0.09 | -0.01 | 0.072 | 1.28 |
| 2 (30) | -0.03 | 0.24 | 0.18 | -0.01 | 0.049 | 1.29 |
| 3 (46) | -0.04 | 0.26 | 0.21 | -0.01 | 0.038 | 1.31 |
| 5 (76) | +0.03 | 0.25 | 0.13 | +0.06 | 0.031 | 1.31 |
| 10 (152, one set) | +0.07 | n/a | n/a | +0.09 | 0.023 | 1.31 |

Per-sample CV of f across random k-chunk subsets: k=1 0.125, k=2 0.081, k=5 0.039 (roughly 1/sqrt(k), so adding sites does shrink sampling noise in f).
The k=5 and k=10 subsets overlap heavily with each other, so those rows are not independent; the SD for k=5 is from subsets of the same 10 chunks.

## Answer

Observed: r vs age sits near zero at every pooling level (mean -0.04 to +0.03 for k=1..5, +0.07 for all of chr10), and does not trend upward with k.
It varies a lot between independent chunks (-0.31 to +0.45), with the sign flipping. The ancient-vs-control offset is present in all chunks (1.10-1.62) with no
positional trend. The noise in per-sample f falls roughly as 1/sqrt(sites), but the reproducible between-sample signal is small (split-half 0.13).

Would more chromosomes plausibly help? **Probably not, unclear at the margin (inference).** If the estimator carried a real age signal diluted by sampling noise,
r should rise with k; it did not between 15 Mb and 152 Mb (a tenfold increase in sites). That suggests the limiting factor is not site count but a bias/signal problem
(for example the shared ancient-vs-control offset, ARG error, damage, or panel/ancestry mismatch), which more chromosomes of the same data would not remove. This
is a 1-draw, 1-chromosome, n=46 result with wide spread in r (SD about 0.24 per chunk), so a modest true correlation (about 0.2-0.3) cannot be excluded; the data
rule out neither that nor a purely null result. Adding chromosomes would shrink f noise (observed) but would shrink the offset not at all (offset is constant at about 1.3 across the pooling levels).

## Notes
* Only draw 75 was used; draw-to-draw variation was not examined here.
* `sd_of_T...` column in the console header is a placeholder and not computed.
