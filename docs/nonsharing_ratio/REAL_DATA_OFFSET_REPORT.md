# Real-data check of the non-sharing ratio: a constant ancient-vs-modern offset (2026-11-12)

Status (updated 2026-11-12, end of day): real data = one chromosome (chr10), 10 draws. Simulations = constant size and bottleneck.
Aggregates only; no sample-level ages or dates are stored here.
The population-structure explanation below is a hypothesis that has **not** been tested on real data. Simulations with ancient
samples from a population that split 6,000 generations ago are in section 8 (controls from the panel population and from the
ancient population); they show that structure can produce a constant offset in simulation, not that it explains the real one.

If you are coming back to this cold: read this file, then sections 6 to 8 (simulations; whether more data helps; structure), then the open
GitHub issues (label `non-sharing-ratio`).

## Summary

- Run on real data (chr10, held-out modern controls, 50 trusted ancient samples), the estimator did
  **not** recover radiocarbon ages: r between T_hat and calibrated age is about 0 in every setting tried.
- With child-node dating, every ancient sample comes out near 10,000 generations whatever its true age
  (true ages are roughly 0.2-2.2k). The estimate looks like a constant offset, not a clock.
- The offset: ancient samples fail to share old derived alleles ~31% more often than held-out modern
  controls (f = 0.418 vs 0.318, H = log ratio = 0.27).
- Read depth does not explain it. Simple reference bias does not explain it. A population-structure
  (or drift) difference between the ancient samples and the panel/controls is the leading hypothesis.
- In msprime simulations with SINGER+POLEGON ARGs (constant size and a bottleneck, no structure) the same estimator
  recovers ancient ages at A0 = 10k and 20k (r 0.90-0.94, bias about -350 to -700 generations). The real-data null is
  therefore not a failure of the method on clean data (section 6).
- More posterior draws (10 tried) and more sites (chr10 split into chunks) do not move r (section 7).

## Setup

- ARG: stand-in SINGER+POLEGON draw 75 on chr10, built from 20 panel haplotypes (6 lines held out as
  controls). Files: `batches/nonsharing_ratio_heldout_20261111/` (see its README.md).
- Targets: the 50 "trusted" ancient samples (`trusted50`); 46 have a calibrated date in the metadata table
  (3 have only an uncalibrated 14C age, 1 is waiting for a date), so age comparisons use n = 46.
- Site selection: sites dated above A0, by edge midpoint (original) or child-node time (`--age-mode child`).

## 1. Age recovery (46 dated samples, draw 75)

| Setting | Sites kept | Pearson r vs true age | Median T_hat (true median 852) |
|---|---|---|---|
| A0 = 12,000, edge midpoint | 92% | -0.03 | 609 |
| A0 = 30,000, edge midpoint | 88% | -0.02 | 732 |
| A0 = 20,000, child node | 34% | 0.07 | 9,969 |
| + panel-fixed sites, REF polarity | 34% | 0.07 | 9,970 |
| + panel-fixed sites, ancestral allele from the full 26-haplotype ARG | 35% | 0.06 | 10,316 |

Edge-midpoint dating kept most sites because a tip-ward edge with a very deep parent has a huge midpoint
even when its child node is young (child node > 30,000 for only 29% of sites). Fixed sites changed nothing
(75 and 27,874 sites added, of ~1M selected).

## 2. The offset (A0 = 20,000, child node)

Controls f = 0.318 (range 0.307-0.329 over 6 controls); ancient median f = 0.418 (range 0.391-0.438).
Ratio 1.31, H = 0.27, for essentially every ancient sample.

## 3. Evidence against read depth

Median f by mean depth of the ancient sample (depths from the metadata table):

| Mean depth | n | Median f | Ratio to controls |
|---|---|---|---|
| under 0.3x | 14 | 0.4165 | 1.31 |
| 0.3-1x | | 0.4189 | |
| 1x and over | 21 | 0.4190 | 1.32 |

Spearman(depth, f) = +0.26 (p = 0.07), the wrong sign for a low-depth artefact. The 7.8x sample is no
cleaner than the 0.04x ones. Caveat: for pseudo-haploid calls (one read per site) error and mapping bias
barely depend on depth, so this test cannot exclude depth-independent artefacts.

## 4. Evidence against simple reference bias

Sites split by whether the derived allele is the VCF REF or ALT (REF/ALT from the controls' VCF):

| Derived allele is | Sites | Controls f | Ancient median f | Ratio | Excess (ancient - controls) |
|---|---|---|---|---|---|
| ALT (ancestral = REF) | 433,619 | 0.552 | 0.612 | 1.11 | +0.060 |
| REF | 548,925 | 0.223 | 0.302 | 1.36 | +0.079 |

Reference bias (calls pulled toward REF) predicts extra non-sharing where the derived allele is ALT and
**less** where it is REF. The excess is positive in both classes and larger where the derived allele is REF.
Symmetric genotype error fits poorly too (rough arithmetic: excess = e(1-2f), about 0 at f = 0.55, so it
cannot give +0.06 there while giving +0.079 at f = 0.22). Not testable here: error rate measured at sites
where all moderns are derived and the ancient is ancestral, because ancient sites were ascertained from
ARG SNPs and such sites are not in the ancient VCF.

## 5. Hypothesis: population structure / drift (untested)

The ancient samples (South American maize) are more diverged from the inbred-line panel and controls than
the controls are from the panel. Old derived alleles carried by the panel references are then missing in
the ancient genomes more often, independent of the ancient sample's age, depth or allele orientation.
That would give a positive, roughly constant, age-independent offset as observed. Predictions to test:

- The offset should be flat across mutation age (TODO item "population structure check").
- It should shrink for ancient samples genetically closer to the panel, and for controls chosen from a
  more distant ancestry group (3 of the 6 held-out lines are tropical types; only a narrow spread of
  f among controls was seen, 0.307-0.329).
- Using controls from the same population as the ancient samples, if available, should remove it (in the section 8 simulation it removes the offset but not the age error).
- Simulated data with no structure or artefacts should recover ages: it does (section 6).
- Simulated data with ancient samples from a diverged population should show a flat positive offset:
  it does in simulation (ratio 1.14-1.16 at A0 10k/20k; section 8). Controls from the ancient population cancel it in simulation but age
  recovery was not restored.

If structure is the cause, the ratio method as formulated assumes the ancient and control samples differ
only by the age of the ancient lineage, which real ancient data violate.

## 6. Simulations without structure: the method works (SINGER + POLEGON, 2026-11-12)

Design (both batches): 5 replicates, 20 Mb, 36 modern haploid samples (26 in the ARG, 10 held-out controls) and 10 ancient
haploid samples each (50 in total, true ages stratified-uniform over 0-6,000 generations, not in the ARG). Mutations
private to one ancient sample are dropped before inference. 5 POLEGON-dated draws per replicate, averaged. Child-node A0.
Polarity: singer-snakemake randomises REF/ALT before SINGER (stated by the user, not re-read in the recipe), so polarity
is treated as in real data; the config has `polarised: False`.

Mean over draws, 50 ancient samples, bias = inferred - true age in generations (c(t) from ARG lineage counts):

| Demography | A0 | Mean bias | r | Slope | Bias in oldest bin (4.5-6k) |
|---|---|---|---|---|---|
| Constant size | 10k | -714 | 0.948 | 0.81 | -1,077 |
| Constant size | 20k | -611 | 0.904 | 0.83 | -807 |
| Constant size | 50k | +4,025 | 0.917 | 2.58 | +8,253 |
| Bottleneck (0.3 Ne at 500 ago, back to Ne by 300 ago) | 10k | -438 | 0.937 | 0.85 | -882 |
| Bottleneck | 20k | -352 | 0.927 | 0.90 | -579 |
| Bottleneck | 50k | +5,370 | 0.936 | 2.94 | +9,880 |

On the true tree the bias is positive: +198 / +371 / +5,876 (constant) and +655 / +880 / +7,461 (bottleneck) at A0 10k / 20k / 50k.
The bottleneck raises the true-tree bias, mostly in the oldest bin (cause untested). A0 = 50k is biased high even on the true
tree (only about 4% of sites survive; cause untested). Using the true demography for c(t) changes the ARG bias by about 190
generations at A0 10k and 20k. Controls used as pseudo-targets date as modern (median T_hat 0-54). Each condition has only
5 replicates, so small differences between the two demographies are not resolved.

Where: `batches/nonsharing_ratio_realsims_20261112/` (constant size) and
`batches/nonsharing_ratio_realsims_bottleneck_20261112/` (bottleneck); each has README.md, COMMANDS.md, `summary_by_a0.tsv`,
`inferred_vs_true.png`. The simulations have no ancient-vs-modern structure, masks, missing data, damage or variable
Ne other than the one bottleneck.

## 7. Would more draws or more chromosomes help the real data? Probably not (checked 2026-11-12)

More draws: `batches/nonsharing_ratio_heldout_20261111/more_draws_check/` (10 POLEGON-dated draws, 75-138, same child-node
A0 = 20,000 settings, 46 dated samples). Per-sample f and T_hat are stable across draws (ICC 0.95 and 0.91; between-draw
Pearson about 0.96). The offset is stable: ancient median f 0.4189 +/- 0.0013 (SD over draws) against controls 0.3176 +/- 0.0014,
and every ancient sample is above the highest control in every draw. r with calibrated age is not improved by averaging:
0.043 for T_hat averaged over 10 draws, 0.004-0.072 for single draws, flat for 1, 2, 5 and 10 draws.

More chromosomes: `batches/nonsharing_ratio_heldout_20261111/block_check/` (chr10 draw 75 in 10 chunks of 15-20 Mb).
r with age ranges from -0.31 to +0.45 across chunks (mean -0.03). Pooling k random chunks shrinks the spread of per-sample f
about as 1/sqrt(k), but mean r stays near zero (k=1 -0.02, k=5 +0.03, all of chr10 +0.07) and the offset stays about 1.3 in
every chunk. Per-sample f is poorly reproducible between disjoint halves of chr10 (mean correlation 0.13). The
offset varies along the chromosome mostly because control f varies.

Inference (not tested): the offset is shared by all draws and chunks, and the sample-to-sample variation that is reproducible
across draws is unrelated to age, so more data of the same kind probably will not recover ages. Limits: one chromosome; the
10 draws come from one SINGER chain; the per-chunk spread of r (SD 0.24) does not rule out a true r of 0.2-0.3; only A0 = 20,000
with child-node ages was run across draws.

## 8. Structure simulations: split population, controls from A versus from B (2026-11-12/13)

Design: ancestral Ne 20,000 (diploid-equivalent) splits 6,000 generations ago into A and B, no migration; both have the 0.3 Ne
bottleneck 500 generations ago with recovery by 300. Panel (26, the ARG) and the held-out controls from A; the 10 ancient samples per
replicate from B (ages stratified-uniform 0-6,000). 5 replicates, one split time, SINGER+POLEGON, 5 draws, child-node A0, ARG c(t).
Batches: `batches/nonsharing_ratio_realsims_split6k_bottleneck_20261112/` (controls from A only) and
`batches/nonsharing_ratio_realsims_split6k_Bcontrols_20261113/` (new simulations with 10 controls from A and 10 from B; the two arms differ only in
`--controls-vcf`). Offset = median ancient f / mean control f.

| Batch / arm | A0 | Offset (ARG) | r | Mean bias | Slope |
|---|---|---|---|---|---|
| No structure (bottleneck batch) | 10k / 20k | not computed | 0.937 / 0.927 | -438 / -352 | 0.85 / 0.90 |
| Split, controls A (earlier batch) | 10k | 1.135 | 0.181 | +1569 | 0.04 |
| Split, controls A (earlier batch) | 20k | 1.146 | 0.214 | +1914 | 0.08 |
| Split, controls A (earlier batch) | 50k | 1.588 | 0.205 | +13877 | 0.24 |
| Split, controls A (new batch) | 10k | 1.144 | 0.065 | +1905 | 0.02 |
| Split, controls A (new batch) | 20k | 1.155 | 0.108 | +2281 | 0.05 |
| Split, controls A (new batch) | 50k | 1.596 | 0.083 | +14094 | 0.09 |
| Split, controls B (new batch) | 10k | 1.004 | 0.111 | -2716 | 0.02 |
| Split, controls B (new batch) | 20k | 1.007 | 0.122 | -2576 | 0.03 |
| Split, controls B (new batch) | 50k | 1.034 | 0.071 | -1632 | 0.05 |

On the true tree the A-control offsets are 1.20 / 1.22 / 1.75 (earlier batch) and 1.22 / 1.24 / 1.74 (new batch); B-control 1.005 / 1.006 / 1.036.
Controls as pseudo-targets: A-controls and B-controls each date as modern against their own set (median T_hat 0-2). B-controls dated against the A-control mean
read as about 4,400 generations old in the ARG (about 6,900-7,400 on the true tree; A0 10k and 20k), so a modern B sample carries the same offset as the B ancient samples.
With B-controls, 29 of 50 targets have H<=0 in at least one draw at A0 10k and 20k (0 with A-controls).

Measurements: (1) with controls from A the ancient samples from B show a nearly constant excess non-sharing of 1.14-1.16 at A0 10k/20k (similar in all
5 replicates in both batches) and ages are not recovered (r 0.07-0.21, versus 0.93 without structure); (2) controls from B cancel the offset (ratio 1.00-1.01) but
ages are still not recovered (r 0.07-0.12, estimates collapse toward 0).
Inferences, not tested: structure alone can produce a roughly constant offset like the one seen on the real data, and matched controls remove the offset
but did not restore age signal here; a possible reason is that every ancient age is younger than the split while the old alleles used (A0 >= 10k)
are older than it, so age changes sharing with the A panel by little (untested; would need ancient ages older than the split or a panel from B).
Caveats: 5 replicates per condition (per-replicate r runs from -0.26 to +0.71); one split time; no migration; B-controls are
exactly the ancient population, the best case; real moderns from the ancient samples' region would not be exact. The earlier and new controls-A rows differ
(r 0.18-0.21 vs 0.07-0.11) because the simulations differ; whether this is beyond noise was not tested. Nothing here is a statement about the real data:
real offsets were 1.31 at A0 20,000 and the real ancient and panel populations, split time and gene flow are unknown.

## Caveats

- Real-data runs: chr10 only, a 20-haplotype panel (not the full 26), 6 unbalanced controls. Radiocarbon-age comparisons
  were first made on draw 75 alone and then on all 10 draws (section 7).
- Ancient sample identities were first mapped to calibrated dates by VCF column order. On 2026-11-12 the mapping was
  checked by sample name against the metadata table (all 50 columns match; `private/check_vcf_order.py`).
  4 of 50 samples have no calibrated date.
- 4.6% of ancient allele matches were dropped for allele mismatch; not investigated.
- Depths and dates come from the metadata table (`meta_v10_9July2026.csv`), kept only in the gitignored `private/` directory.
- All comparisons use time in generations equal to years (1 generation = 1 year assumed).

Scripts: `batches/nonsharing_ratio_heldout_20261111/scripts/` (`04_verify_run_draw.sh`,
`07_prepare_fixed_heldout.sh`, `08_polarity_from_full_arg.py`, `09_ref_alt_split.{sh,py}`),
`validation/nonsharing_ratio_a0_childtime_probe.py`. Simulation scripts: `batches/nonsharing_ratio_realsims_*/scripts/`.
Checks in section 7: `more_draws_check/scripts/`, `block_check/block_check.py`. Dates and the sample-to-date mapping used in
the radiocarbon comparisons are in the gitignored `private/` directory (`cmp.py`, `depth.py`, `dates.txt`, `ages.txt`); they are
not in git. The metadata table is in `private/meta_v10_9July2026.csv` (converted from `/quobyte/jrigrp/beil/logan_collab/original_files/meta_v10_9July2026.xlsx`); the sample-to-date mapping was checked by sample name (the VCF sample names are the metadata `FAM_Alias` values; all 50 columns match, 0 mismatches, `private/check_vcf_order.py`).
