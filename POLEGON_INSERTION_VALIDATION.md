# POLEGON and ancient-lineage insertion validation

This document explains the ancient-lineage insertion likelihood, the
constant-population-size simulations used to test it, the diagnostics used to
separate likelihood errors from ARG-dating errors, and the initial results.

Development and validation are on branch
`feature/polegon-insertion-validation`. Large simulation and posterior files
are not stored in Git. They are distributed separately in
`polegon_constant_size_bundle_20261109.tgz`, whose root README describes every
bundled artifact and whose `SHA256SUMS` file checks its contents.

## 1. Ancient-lineage insertion approach

### Question and conditioning

We observe pseudo-haploid genotypes from an ancient sample whose age $T$ is
unknown. At each polymorphic site, a modern-panel ARG supplies a local tree and
a mutation-bearing edge. Rather than reducing the modern genealogy to an
allele count, the insertion model conditions on the complete local modern tree
and integrates where an additional ancient lineage first coalesces into it.

For site $i$, let the mutation-bearing edge have child $v_i$, parent
$w_i$, and time interval

$$
E_i=[t_{c,i},t_{p,i}),\qquad t_{c,i}<t_{p,i}.
$$

At time $u$ before the present, define

- $\kappa(u)$: the number of modern-tree lineages extant at $u$; and
- $\delta_i(u)$: the number of those lineages descended from $v_i$.

With diploid effective size $N_e(u)$, the inserted lineage coalesces with
each extant modern lineage at rate $1/[2N_e(u)]$. Starting at ancient sample
age $T$, its probability of surviving without coalescing until $u$ is

$$
S_T(u)=\exp\left[-\int_T^u \frac{\kappa(s)}{2N_e(s)}\,ds\right].
$$

If the mutation occurred at exact time $z_i$, the ancient haplotype is
derived only if $T<z_i$ and its first coalescence before $z_i$ enters one
of the $\delta_i(u)$ descendant lineages:

$$
q_i(T\mid z_i,E_i,\mathcal T)=
\mathbf 1\{T<z_i\}
\int_T^{z_i}S_T(u)\frac{\delta_i(u)}{2N_e(u)}\,du.
$$

Under the infinite-sites model, mutation time is uniform on its assigned edge,
so the implemented probability marginalizes it:

$$
q_i(T\mid E_i,\mathcal T)=
\frac{1}{t_{p,i}-t_{c,i}}
\int_{t_{c,i}}^{t_{p,i}}
q_i(T\mid z,E_i,\mathcal T)\,dz.
$$

`pipeline/insertion_likelihood.py` evaluates this analytically between tree
events and demographic epoch boundaries. It supports constant or
piecewise-constant $N_e(t)$.

### Discordance and per-site likelihood

For effective symmetric discordance probability $\epsilon$, the probability
of observing a derived allele is

$$
r_i(T)=\epsilon+(1-2\epsilon)q_i(T).
$$

For pseudo-haploid observation $y_i\in\{0,1\}$, the site log likelihood is

$$
\ell_i(T)=y_i\log r_i(T)+(1-y_i)\log[1-r_i(T)].
$$

On inferred ARGs, $\epsilon$ is an effective discordance parameter. It can
absorb genotype error, polarization error, topology error, or mutation-edge
error and should not automatically be interpreted as sequencing error.

### ARG time-scale calibration

Candidate ARG times can be multiplied by a positive correction
$\varsigma$. The ARG node times and mutation-edge endpoints are scaled inside
the insertion likelihood, while candidate sample ages and $N_e(t)$ remain in
real generations. Post-hoc multiplication of final MAP ages is not equivalent.

For dated training samples $D$, each with radiocarbon density $R_j(T)$, the
joint calibration surface is

$$
L(\varsigma,\epsilon)=
\prod_{j\in D}\int
L_j(T\mid\varsigma,\epsilon)R_j(T)\,dT.
$$

For a held-out sample, calibration uncertainty is marginalized rather than
replaced by the joint MAP:

$$
p(T_j\mid G_j,D)\propto p(T_j)
\sum_{\varsigma,\epsilon}
L_j(T_j\mid\varsigma,\epsilon)
p(\varsigma,\epsilon\mid D).
$$

Cross-validation estimates the calibration surface using training dates only;
the held-out radiocarbon density is used only for scoring.

### Posterior ARG draws and chromosomes

Sites are multiplied within one ARG draw before draws are averaged:

$$
\log\mathcal L_c(T)=
\operatorname{logsumexp}_g
\left[\sum_{i\in\mathcal I_c}\ell_{ic}^{(g)}(T)\right]-\log G_c.
$$

Independently inferred chromosomes are then combined by adding their
log-likelihoods. Averaging each site's probability across ARG draws before the
site product would destroy chromosome-wide posterior dependence and is not the
implemented model.

## 2. Constant-size simulations

### Design

The constant-size suite contains three independent panels:

- `constant_r1`
- `constant_r2`
- `constant_r3`

Every panel uses the same model:

| Property | Value |
|---|---:|
| Diploid effective population size | 20,000 |
| Chromosomes per panel | 3 |
| Chromosome length | 50 Mb |
| Modern haploid samples | 100 |
| Serial samples | 21 |
| Serial-sample ages | 0, 500, ..., 10,000 generations |
| Mutation rate | $10^{-8}$ per base per generation |
| Recombination rate | $10^{-8}$ per base per generation |
| Ancient-call error in truth simulation | 0 |
| Ancestry model | Hudson coalescent |
| Mutation model | Infinite-sites nucleotides |

The age-zero serial sample is a control, leaving 20 genuinely ancient samples.
The three chromosomes within a panel are independent ancestry and mutation
realizations for the same design. The `r1`--`r3` panels use different random
seeds and provide three independent repeats.

### Files produced per chromosome

`validation/simulate_variable_ne_error.py` writes:

- `truth_with_ancient.trees`: full truth containing modern and ancient nodes;
- `known_modern_arg.trees`: the same truth simplified to modern samples, with
  true node and mutation times;
- `all_samples.vcf.gz`: modern and ancient calls;
- `mutation_truth.tsv.gz`: mutation times, nodes, frequencies, and calls;
- `metadata.json`: seeds, sample-node IDs, ages, rates, and versions; and
- `constant_ne_epochs.tsv`: the demographic history used by the likelihood.

POLEGON should be tested on `known_modern_arg.trees`. The full
`truth_with_ancient.trees` file supplies complete simulation truth but is not
the modern-panel ARG that the production dating workflow receives.

### Inferred ARG posterior

For production inference, every 50 Mb chromosome was divided into 50 roughly
1 Mb chunks. For each posterior draw:

1. SINGER inferred each chunk.
2. POLEGON dated each chunk separately.
3. Chunks were merged into a chromosome tree sequence.
4. The three independently sampled chromosome ARGs were concatenated for
   storage, while their likelihoods remained separate.

There are 100 saved posterior draws for every chromosome and 100 combined ARGs
per panel. The share bundle contains all 1,200 posterior tree files, plus the
raw SINGER and POLEGON text tables and logs for representative draw 50.

## 3. Validation workflow

### A. Unit and integration tests

The focused suite covers exact insertion probabilities, time scaling,
radiocarbon integration, leakage-safe cross-validation, discrete quantiles,
chromosome-specific marginalization, true-topology POLEGON conversion, and
production-matched chunk merging:

```bash
pytest -q \
  tests/test_radiocarbon_calibration.py \
  tests/test_joint_scale_epsilon_calibration.py \
  tests/test_joint_scale_epsilon_recovery.py \
  tests/test_age_dependent_scale_diagnostic.py \
  tests/test_singer_polegon_diagnostic.py \
  tests/test_chromosome_scale_profile.py \
  tests/test_per_chromosome_calibration.py \
  tests/test_true_topology_polegon.py \
  tests/test_chunked_true_topology_polegon.py
```

This suite currently passes 21 tests on Farm.

### B. Controlled true-ARG recovery

`validation/joint_scale_epsilon_recovery.py` imposes known scale compression or
expansion and known discordance on true ARGs. The recovery matrix includes
baseline, 0.8 compression, 1.25 expansion, $\epsilon=0.02$, and a joint
compression/discordance condition across constant, bottleneck/recovery, and
recent-growth simulations.

All 45 recovery gates passed. This shows that the scale/epsilon grid,
radiocarbon integration, and cross-validation arithmetic work when topology
and node times are known.

### C. Inferred-ARG demographic matrix

The same calibration was then applied to POLEGON-dated SINGER draw 50 across
three chromosomes in nine demographic runs. Constant-size results were
consistent across `r1`--`r3` but substantially biased, despite selecting the
same scale and $\epsilon$.

### D. Raw SINGER versus POLEGON on fixed inferred topology

`validation/singer_polegon_diagnostic.py` reconstructs the raw-SINGER tree with
the production merge script and compares it with the POLEGON-dated tree while
holding topology, mutations, masks, and chunks fixed. In the recent-growth
pilot, POLEGON improved age estimation relative to raw SINGER, even though
both preferred scale 0.9. This ruled out a simple claim that POLEGON uniformly
over-corrects all SINGER times.

### E. Chromosome-specific calibration

`validation/chromosome_scale_profile.py` showed that chromosomes individually
could prefer scale 0.9 while their joint profile preferred a much smaller
global scale. `validation/per_chromosome_calibration.py` therefore learns a
separate $p_c(\varsigma,\epsilon)$ from training dates for each chromosome,
marginalizes each chromosome over its own calibration surface, and only then
multiplies the age likelihoods. This improved inferred-ARG error but did not
eliminate it.

### F. Production-matched true-topology POLEGON test

`validation/chunked_true_topology_polegon.py` provides the strongest
discriminator. It takes one true modern ARG and:

1. splits it at the exact production SINGER chunk boundaries;
2. retains true times in one arm;
3. runs POLEGON independently per chunk in the other arm;
4. merges both arms using shared sample nodes and duplicated chunk-specific
   internal nodes, matching production; and
5. runs identical chromosome-specific insertion calibration on both arms.

This paired test holds topology, mutations, chunk boundaries, genotype calls,
and likelihood code fixed. Only the node times differ. It has been completed
for recent growth; applying it to the constant suite is the next direct test.

## 4. Initial results

### Constant-size insertion calibration

The table reports five-fold held-out results for 20 dated samples in each
panel. “Oracle” uses true modern ARG node times. “Inferred” uses SINGER topology
and production POLEGON dates for draw 50 across all three chromosomes.

| Panel | ARG | Scale MAP | $\epsilon$ MAP | Bias | MAE | RMSE | Slope |
|---|---|---:|---:|---:|---:|---:|---:|
| constant_r1 | Oracle | 1.0 | 0.00 | -12 | 78 | 88 | 1.009 |
| constant_r2 | Oracle | 1.0 | 0.00 | +18 | 38 | 72 | 1.008 |
| constant_r3 | Oracle | 1.0 | 0.00 | -22 | 68 | 93 | 1.000 |
| constant_r1 | Inferred | 0.8 | 0.01 | -2,402 | 2,402 | 2,822 | 0.495 |
| constant_r2 | Inferred | 0.8 | 0.01 | -2,417 | 2,417 | 2,841 | 0.490 |
| constant_r3 | Inferred | 0.8 | 0.01 | -2,502 | 2,502 | 2,930 | 0.478 |

The true-ARG controls are point-calibrated, showing that the insertion model
can recover age under the simulated constant-size process. The inferred ARGs
remain severely biased young. A fitted global scale and effective discordance
do not solve that bias.

### Strongest current POLEGON discriminator

The production-matched fixed-topology experiment has so far been run under
recent growth:

| Timing arm | Bias | RMSE | Slope | RC MAP within genetic 95% |
|---|---:|---:|---:|---:|
| Chunked true times | +13 | 69 | 1.010 | 30% |
| Same topology and chunks, POLEGON times | -1,247 | 1,459 | 0.762 | 0% |

Chunking and merging did not create the MAP-age bias: the chunked true-time arm
was accurate. Replacing only node times with POLEGON estimates created a large
young bias. Median POLEGON/true node time was close to one, but the 5th--95th
percentile range was approximately 0.59--1.56, demonstrating that global means
and one multiplicative scale can conceal structured local distortions.

This result localizes the MAP-age error to POLEGON dating or the statistical
handoff from POLEGON point dates to the insertion likelihood. It does not by
itself demonstrate a POLEGON software bug. POLEGON used a scalar $N_e$, while
the recent-growth truth was nonconstant, and our likelihood treats its output
node times as exact rather than integrating dating uncertainty.

The true-time arm also had overly narrow intervals despite accurate MAP ages.
That is a separate composite-likelihood uncertainty problem, not the cause of
the systematic MAP bias.

## Code map

- `pipeline/insertion_likelihood.py`: insertion probability equations.
- `pipeline/insertion_real_data.py`: per-draw, per-chromosome likelihoods and
  ARG-draw marginalization.
- `pipeline/joint_scale_epsilon_calibration.py`: joint calibration and
  leakage-safe cross-validation.
- `pipeline/radiocarbon_calibration.py`: radiocarbon integration, summaries,
  and scoring metrics.
- `validation/simulate_variable_ne_error.py`: truth simulations.
- `validation/joint_scale_epsilon_recovery.py`: controlled recovery matrix.
- `validation/singer_polegon_diagnostic.py`: raw SINGER versus POLEGON.
- `validation/chromosome_scale_profile.py`: chromosome-subset profiles.
- `validation/per_chromosome_calibration.py`: chromosome-specific calibration.
- `validation/true_topology_polegon.py`: true-topology POLEGON adapter and
  likelihood grids.
- `validation/chunked_true_topology_polegon.py`: production-matched paired test.

## Current interpretation and next test

The evidence supports three separate conclusions:

1. The core insertion likelihood and joint calibration arithmetic recover
   known perturbations on true ARGs.
2. Inferred ARG node times produce large, reproducible age bias that a global
   scale and $\epsilon$ cannot repair.
3. In the recent-growth paired control, POLEGON timing—not chunking or topology
   alone—was sufficient to introduce most of that bias.

The next decisive experiment is the same production-matched true-time versus
POLEGON-time test on `constant_r1`--`constant_r3`. If POLEGON remains biased
under correctly specified constant $N_e$, demographic mismatch is not a
sufficient explanation. If the bias largely disappears, the scalar-$N_e$
assumption is the leading mechanism.
