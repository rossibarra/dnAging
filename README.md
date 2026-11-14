# sample_age_dating

Date ancient genomes using ancestral recombination graphs (ARGs) of modern genomes.

The current approach is the **non-sharing ratio** (method: "Dating ancient genomes with inferred ARGs",
Yun Deng, 2026-11-11). An ancient lineage of age $T$ cannot coalesce with modern lineages during its first
$T$ generations, so at mutations the ARG says are much older than $T$, an ancient sample fails to share the
derived allele with modern references more often than a modern control does. That excess is the clock.

![Graphical abstract. a: an ancient lineage starts late and cannot coalesce in (0, T). b: at old mutations the
ancient target fails to share the derived allele more often than a modern control. c: the ratio of
non-sharing gives age through the integral of c(t). d: simulated ancient samples, inferred versus true age, for
constant size and for a bottleneck.](docs/nonsharing_ratio/graphical_abstract.png)

*Panels a-c are schematic. Panel d is simulated data (SINGER + POLEGON ARGs, c(t) from the inferred ARG).
Regenerate with `python scripts/build_nonsharing_ratio_abstract.py`.*

## Status (2026-11-13)

| | result |
|---|---|
| Simulated ancient samples, ARG inferred, c(t) inferred | **works**: r 0.90 to 0.95 with true age, mean bias a few hundred generations (table below) |
| Real ancient maize genomes (one chromosome, radiocarbon-dated samples) | **does not yet recover ages**: r about 0 with radiocarbon age; every ancient sample shows a constant ~31% excess non-sharing vs held-out modern controls |
| Cause of the real-data offset | open. Not read depth, not simple reference bias; more posterior draws and more chromosomes of the same data do not help (inferred from one chromosome). Population structure between the ancient samples and the modern panel is the leading hypothesis and has **not** been tested on real data. In simulation, ancient samples from a population that split from the panel's 6,000 generations ago give a constant offset of about 1.14-1.16 when controls come from the panel's population. Controls from the ancient samples' population remove the offset (ratio 1.00-1.01) but ages are still not recovered (r about 0.1). So controls from the ancient samples' region may fix the offset without fixing ages; this is untested on real data (5 replicates, one split time). |

So: the estimator is validated on simulated data, including with ARG inference error and with c(t) taken from the
inferred ARG. It has not been shown to work on real ancient genomes, and a clean-data simulation succeeding does not
tell us why the real data differ.

## Simulation results

Setup (both demographies): msprime, 20 Mb, 5 replicates; 26 modern haplotypes build the ARG and 10 held-out modern
haplotypes are controls; 10 ancient haploid samples per replicate (50 per demography, true ages stratified-uniform over
0-6,000 generations, not in the ARG). The ARG is inferred by SINGER and re-dated by POLEGON (5 draws per replicate,
estimates averaged). **c(t) comes from the inferred ARG's lineage counts**, not from the true demography. Only
mutations dated above the cutoff A0 (on the child-node time of the mutation edge) are used. REF/ALT are randomised before
SINGER, as for real data. Mean bias is inferred minus true age, in generations.

| Demography | A0 | Pearson r | Mean bias | Slope |
|---|---|---|---|---|
| Constant size | 10,000 | 0.948 | -714 | 0.81 |
| Constant size | 20,000 | 0.904 | -611 | 0.83 |
| Bottleneck (0.3 Ne at 500 gens ago, back to Ne by 300) | 10,000 | 0.937 | -438 | 0.85 |
| Bottleneck | 20,000 | 0.927 | -352 | 0.90 |

- Using the true demography for c(t) instead of the ARG's own changes the mean bias by about 190 generations at
  A0 10,000 and 20,000 (bottleneck run: -438 to -251 and -352 to -164), and r is unchanged. So estimating c(t) from the ARG
  costs a little accuracy, not the signal.
- Controls used as pseudo-targets date as modern (median estimate 0 to 54 generations).
- **A0 = 50,000 is biased high** (bias +4,000 to +5,400) even on the true tree sequence, where only about 4% of sites survive the
  cutoff; the cause is untested. Use A0 near 10,000-20,000 for ages up to a few thousand generations.
- Caveats: 5 replicates per demography, so the two demographies are not resolved from each other; one bottleneck only; no
  population structure, missing data, damage or genotype error in these simulations.

Full numbers, per-bin biases, true-tree comparisons and figures: [docs/nonsharing_ratio/REAL_DATA_OFFSET_REPORT.md](docs/nonsharing_ratio/REAL_DATA_OFFSET_REPORT.md)
(section 6), `batches/nonsharing_ratio_realsims_20261112/` (constant size) and
`batches/nonsharing_ratio_realsims_bottleneck_20261112/` (bottleneck). The per-sample data behind panel d is
[results/nonsharing_ratio/sim_inferred_vs_true.csv](results/nonsharing_ratio/sim_inferred_vs_true.csv).

## More documentation

- [docs/nonsharing_ratio/REAL_DATA_OFFSET_REPORT.md](docs/nonsharing_ratio/REAL_DATA_OFFSET_REPORT.md): the real-data result, what was ruled out, and the simulation results in full.
- [docs/nonsharing_ratio/RUNNING.md](docs/nonsharing_ratio/RUNNING.md): how to run the estimator on real data.
- [PROJECT_MAP.md](PROJECT_MAP.md): where things are in the repository.
- Open work is tracked in the [issues](../../issues?q=is%3Aissue+label%3Anon-sharing-ratio) (label `non-sharing-ratio`).

Code: `pipeline/nonsharing_ratio.py` (site selection, $f$, $H$, inversion), `pipeline/nonsharing_ratio_demography.py` (c(t) from
lineage counts or a demes file), `pipeline/nonsharing_ratio_io.py` (VCF reading, polarisation),
`pipeline/nonsharing_ratio_real_data.py` (`run-draw`, `merge`, `prepare-fixed` CLI). Tests: `tests/test_nonsharing_ratio*.py`
(run on a compute node). Environment: [environment.yml](environment.yml) (`conda env create -f environment.yml; conda activate dnaging`).

## The method in one paragraph

For each ancient target and each modern reference carrying the derived allele at a site the ARG dates above A0, record whether
the target also carries it. $f$ is the fraction of (site, reference) pairs where it does not. The same quantity for held-out modern
controls, $\bar f_{\rm mod}$, is the baseline. Under neutrality and a known coalescence rate $c(t)$,
$f_{\rm anc} / \bar f_{\rm mod} = \exp\!\int_0^T c(t)\,dt$, so $\hat H = \log(f_{\rm anc}/\bar f_{\rm mod})$ and $\hat T$ solves
$\int_0^{\hat T} c(t)\,dt = \hat H$. Posterior ARG draws are averaged; block bootstrap gives intervals. Time is in generations.

## Data and privacy

Sample-level real data (ancient VCFs, radiocarbon dates, sample names, metadata) are **not** in this repository. They are kept in the
gitignored `private/` and `radiocarbon_goodvcfs/` directories. Everything committed for the real-data work is aggregate. The
simulated data are fully shareable.
