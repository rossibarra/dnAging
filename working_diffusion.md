# Working notes: the diffusion approach (`main`)

Status as of 2026-09-10. Live working document: **status lines are the point**, so
update them rather than appending. Companion: [working_betabinom.md](working_betabinom.md)
for the alternative likelihood. Spec is [MATH.md](MATH.md); code is
[precompute_freq_trajectory_moments.py](precompute_freq_trajectory_moments.py)
(table build) and [posterior_sample_age_infer.py](posterior_sample_age_infer.py)
(inference).

## Scope and priority

**Solve the bias in perfect simulated data first.** Real-data concerns —
genotyping error, contamination, reference bias, damage — are deferred until the
estimator is unbiased on error-free simulations with true ARGs. See "Deferred:
real data" at the end for what has already been worked out and parked.

This has a sharp consequence for the eps hypotheses. On perfect simulated data
there is **no** genotyping error and the ARG is true, so eps has no physical role
at all: it is purely a numerical regulariser for the hard "carriage is impossible"
wall. And that wall never fires — on the `betabinom` branch, mutations postdating
the sample produced **0 carriers in 41,395 opportunities**. So on this data eps
should simply be set as small as is numerically safe, and any residual bias is
**not** eps. H0 reduces from a research question to a one-line setting.

## Where this stands

The likelihood is sound and, at adequate precision, is the **better of the two
approaches** — see the head-to-head below. What was wrong with it was arithmetic
and plumbing, not model. Six commits today; the numerics are now believed correct
and the inference module has had one review pass.

The bias is **not** resolved, but two new controls localise it sharply. The shared
Bernoulli age likelihood is calibrated when given the true population frequency
at every candidate time, and the production diffusion frequency lookup is
essentially unbiased when given the exact simulated mutation time. The remaining
target is therefore the replacement of an exact mutation time by an ARG edge
interval, including what information from that edge/tree is retained and how the
unknown time is marginalised. The paired run confirms that this transition alone
creates a +364-generation MAP bias at Ne=50K.

### Known-frequency end-to-end control — PASSED

Ten thousand independent one-base SLiM 5.2 Wright-Fisher simulations were run at
constant diploid Ne = 10,000 for 40,000 generations, with neutral mutation rate
1e-6 and no selection. Each replicate was treated as one independent biallelic
SNP. For each candidate generation the analysis used the **true simulated
derived frequency**, collapsed across recurrent mutation IDs, and the exact
haploid likelihood

    log L(T) = sum_i [g_i log p_i(T) + (1-g_i) log(1-p_i(T))].

There was no diffusion, beta-binomial approximation, ARG uncertainty, or
genotyping-error term. All 10,000 seeds were retained and unique
(202609090001--202609100000), in both filenames and TSV columns. The simulation
actually produced seven sample ages, not six:

| true age | MAP | equal-tailed 95% interval |
|---:|---:|---:|
| 6000 | 5993 | 5948--6079 |
| 5000 | 5018 | 4979--5051 |
| 4000 | 3891 | 3840--4001 |
| 3000 | 3005 | 2970--3056 |
| 2000 | 1985 | 1930--2008 |
| 1000 | 1000 | 941--1060 |
| 500 | 479 | 376--541 |

Coverage was 7/7, MAP mean absolute error 25 generations, and MAP RMSE 43. The
4000-generation posterior is the closest to a miss (true age at CDF 0.9741), so
this is a strong end-to-end sanity check, not a claim of formal 95% coverage from
seven non-independent age points. Code and tests are
`direct_frequency_age_infer.py` and `tests/test_direct_frequency_age_infer.py`;
the 100-way SLURM calculation used
`slurm/run_direct_frequency_age_chunks.sbatch` and
`slurm/merge_direct_frequency_age.sbatch`. Compact committed artifacts are the
[summary](slim_single_site_direct_age/results/age_summary.tsv),
[posterior plot](slim_single_site_direct_age/results/age_posteriors.png), and
[MAP calibration plot](slim_single_site_direct_age/results/age_calibration.png).

### Exact-mutation-time diffusion controls — PASSED

Two independent forward/coalescent controls now test the frequency conditional
rather than bypassing it.

**SLiM focal-lineage control.** The recurrent-mutation one-base simulations were
reanalysed by drawing a joint 26-haplotype modern panel and selecting at most one
panel-polymorphic mutation lineage per replicate. Each focal lineage used its
exact `origin_generation`. Of 10,000 independent replicates, 1,236 contained a
usable modern SNP; 6,656 had no extant lineage and 2,108 had no lineage
polymorphic in the sampled panel. Mean estimated-minus-true frequency was
+0.00197 (MAE 0.0577, RMSE 0.0953). All seven true sample ages lay inside their
95% intervals. MAPs were 5730, 5020, 3800, 2950, 2120, 870 and 330 for true ages
6000, 5000, 4000, 3000, 2000, 1000 and 500. This is supportive but not decisive:
there are only seven shared age points, and recurrent mutation makes the focal
allele versus all other states a less exact match to the infinite-sites model.

**msprime exact-time control.** One hundred new independent 10 Mb simulations
used constant diploid Ne = 50,000, 26 modern haploid samples, one ancient haploid
at a uniform random age from 0--10,000 generations, mutation and recombination
rates 1e-8, Hudson ancestry, and `InfiniteSites(NUCLEOTIDES)` on a continuous
genome. Each replicate retains the full truth tree sequence, a modern-only known
ARG with the ancient haplotype removed, a combined 27-haplotype VCF, exact
`Mutation.time` values, and all three seeds. Inference used eps = 0 and the
production `load_table`, `phi_lookup(tab, d0, t_i, t_i)` and `summarize` functions.

Across 7,354,754 used sites the result is essentially unbiased:

| statistic | result |
|---|---:|
| MAP bias (estimated - true) | +10.1 generations |
| posterior-mean bias | -5.7 generations |
| MAP MAE | 197 generations |
| MAP RMSE | 243 generations |
| MAP regression slope | 0.9915 |
| MAP regression intercept | +51.7 generations |

The ordinary composite-likelihood 95% intervals covered 68/100. That is a
precision result, not a point-bias result: tens of thousands of linked sites per
10 Mb replicate are treated as independent. The near-zero mean bias and near-unit
slope reject the diffusion conditional itself as the source of the sustained age
offset when mutation time is known. Results and code are under
`msprime_exact_time_ne50k/results/` and `msprime_exact_time_validation.py`.

**Paired edge-interval test — COMPLETE.** Jobs 38223776/38223803 reused those
exact simulations, VCF calls, known ARGs, table and eps = 0. The sole change was
replacing `t_i` with the known ARG interval `[time(child), time(parent)]` and
applying production's uniform marginalisation. That raised MAP bias from +10.1
to **+363.6 generations** (posterior-mean bias +351.3, MAE 418.8, RMSE 523.0,
slope 1.0194, intercept +268.0; 7,631,501 sites). Naive 95% coverage fell from
0.68 to 0.44. Thus the true edge interval is sufficient to recreate most of the
old Ne=50K production bias (+423 generations), without ARG error, variable Ne,
allele mapping error or epsilon.

## Head-to-head against the beta-binomial

300 msprime simulations (`~/Projects/mutrates/simbatch300`): Ne uniform 10k–100k,
mu log-uniform 1e-9–1e-8, rec = mu x {0.5,1,2}, panel 10–40 haplotypes, true age
100–10,000 generations, 10 Mb, true ARGs, eps = 1e-6. Both estimators on identical
data. Raw:

|                | bias  | med err | SD   | RMSE | slope | r     | 95% cov |
|----------------|------:|--------:|-----:|-----:|------:|------:|--------:|
| beta-binomial  | +1573 |   +1404 | 1666 | 2289 | 1.048 | 0.874 |    0.60 |
| diffusion float64 | -511 | +186   | 2535 | 2582 | 0.841 | 0.692 |    0.60 |
| diffusion 50 digits | +538 | +422 | 1054 | 1182 | 1.062 | 0.946 |    0.78 |

After each estimator gets its own leave-one-out linear calibration in diffusion
time with 1/panel as a covariate (so neither is charged for a bias a user could
calibrate away):

|                | bias | med abs err | RMSE | 95% cov |
|----------------|-----:|------------:|-----:|--------:|
| beta-binomial  | +114 |         758 | 1407 |    0.87 |
| diffusion float64 | +359 |      1418 | 2279 |    0.42 |
| diffusion 50 digits | +80 |       460 |  974 |    0.83 |

The high-precision diffusion's **uncalibrated** RMSE (1182) beats the
beta-binomial's **calibrated** RMSE (1407). Paired, it wins 194/300 (mean abs-err
difference 345 +- 51 generations, t = 6.8). It also leaves no residual age
dependence (t = -0.11 against -1.64), and its panel-size dependence is gone
(t = -1.09).

**The float64 column is the trap.** Its accuracy got *worse* as the panel grew,
which no model can do — that was cancellation, not modelling. Comparing
mpmath-backed beta-binomial code against float64-backed diffusion code was not a
fair test of the math, and the first version of this comparison reached the wrong
conclusion because of it.

## State of the bias

**The sign disagreement is resolved, and the old baseline was measuring a
different estimator.** `bias_ideas.md`'s numbers were produced by a temporary
harness, `/tmp/dnAging-bigsims-estimate.py`, which independently implemented
**denominator-weighted** marginalisation:

```python
def denominator_rows(age, n=26):          # its own alternating-sum denominators
    moments = eng._moms(a / (2 * NE), x0)
    value = sum(c * moments[k] for k, c in eng.coeff[d].items())
...
cden = cumulative_trapezoid(den, fine_age)
cnum = cumulative_trapezoid(den[:, None] * phi, fine_age)
p = numer / denom[:, None]                # a ratio of integrals
```

It imported production `load_table` but never called production `phi_lookup`,
normalizeTEs, or the ARG-draw pipeline. So **its table never described the
shipping uniform estimator**, and `bias_ideas.md` is not a valid baseline for the
bias under diagnosis. (Provenance established from the codex session transcript
and confirmed by reading the script; the file was written 2026-09-08 09:04 and
swept into 54b9113 at 09:47.)

That also explains what looked like a coincidence. Reproducing the weighted arm
through a completely independent route — format5 `log_den` plus analytic per-knot
integration, against dense trapezoid with recomputed denominators — agrees to
within one generation:

| Ne | `bias_ideas.md` | T4 `D_weighted` | diff |
|---:|---:|---:|---:|
| 10,000 | -836 / 890 | -836 / 889 | -0 |
| 50,000 | -564 / 653 | -563 / 653 | +1 |
| 100,000 | -611 / 705 | -612 / 705 | -1 |

Six quantities across three conditions. That is **mutual validation of the
weighted implementation**, not a red flag.

**The actual production (uniform) bias is positive and grows with Ne.** Ten
replicates per Ne, true ARGs, eps = 0, 10 Mb, identical sites and table (T4):

| Ne | bias | RMSE | SD | bias/RMSE | bias in tau |
|---:|---:|---:|---:|---:|---:|
| 10,000 | +250 | 484 | 436 | 0.52 | +0.0125 |
| 50,000 | +434 | 628 | 479 | 0.69 | +0.0043 |
| 100,000 | +651 | 749 | 391 | 0.87 | +0.0033 |

This is consistent with the independent 300-simulation reimplementation (+538 at
50 digits, uniform), so both uniform measurements agree in sign and rough
magnitude. The earlier "(a) versus (b) disagree in sign" puzzle was
uniform-versus-weighted all along.

Three things to note about the real target:

- **It is no longer bias-dominated at small Ne.** bias/RMSE runs 0.52, 0.69,
  0.87. At Ne=10,000 variance is the larger term, so replicate counts now matter
  in a way they did not for the old -836 figure (which was 94% bias).
- **It is neither constant in generations nor in diffusion time.** It grows
  sub-linearly with Ne: +250, +434, +651 in generations while falling 0.0125 ->
  0.0033 in tau. So it is not the fixed-drift offset the beta-binomial shows, and
  not the fixed-generation offset the old table suggested. Whatever causes it
  scales with something else.
- **Weighting makes it worse, not better** — see H2b.

**Measurement caveat carried over from `bias_ideas.md`:** the candidate-age grid
stops at **15,000 generations**. With true ages up to ~9,300 and a positive bias
growing with Ne, estimates approach that ceiling as Ne rises, so an upward bias
measured at high Ne is **partly censored and may be understated**. At Ne=100,000
(+651) there is still headroom; at 200K and above there is not, which is a second
reason the 200K and 500K rows should not be trusted. Widen the grid before
measuring bias at large Ne.

## Hypotheses for the bias

Numbering follows `bias_ideas.md`; **H0 is new**, and the priority order below
differs from that file's.

| # | Hypothesis | Status |
|---|---|---|
| **H0** | **The eps floor, acting asymmetrically in T.** r = eps + (1-2 eps) p, and phi *decreases* with T, so the floor bites hardest where phi is smallest — at large T. It therefore compresses the likelihood's ability to discriminate among old ages rather than shifting it uniformly. Measured below: it removes 44% of the large-T signal at d0=1 while leaving 80% of the small-T signal, and it acts almost entirely on **carried** sites (a 42% change there against 1% on absent sites). | **RETIRED for the zero-error simulation bias once eps=0 is used.** A nonzero eps is misspecification for these simulations and cannot explain residual bias in the corrected eps=0 runs. Retain it as a real-data sensitivity and identifiability question; see T2 and D1. |
| 5b | **Table numerics.** float64 discarded 11% of entries at n=12 and 59% at n=40, wrote 18% of a test build as silent zeros, and was already 0.8% wrong where it passed its own guard. | **FIXED** (d75b9b9, 54b9113, 33b23d6), and **its contribution to this bias is now measured at -0.2 generations** (T4: a legacy float64 table on the identical grid gives +251 against the exact table's +250, with zero dropped groups — the NaN cells lie in table regions this data never queries). The fixes are correct and worth having; they are not the bias. |
| 1 | **Leverage from rare carried alleles.** log p moves fast when p is small, so a few carried rare sites can outweigh many singleton absences. | **TESTED. Leverage is real; rare-site miscalibration is not supported.** Removing all d0=1 sites left 10K essentially unchanged (+241 vs +248 bias) and worsened 50K (+649 vs +423). Carried-only removal moved estimates older and absent-only removal moved them younger, as expected from deleting opposite likelihood terms. T3 independently finds singleton carriage essentially calibrated at 50K--200K. |
| 2a | **Double conditioning on d0.** Once the ARG edge is observed, d0 is determined, so reweighting candidate mutation ages by P(d0 \| t) may condition on the modern count twice. | **DEMOTED** from "leading structural hypothesis". The `betabinom` branch is essentially this fix carried to its limit, and it *loses* at matched precision (calibrated RMSE 1407 vs 974). The exact test in `bias_ideas.md` is cheap and still worth running; the reasoning is sound, but the empirical direction is against it. |
| **2b** | **Marginalisation ORDER over the edge: integral of ratios vs ratio of integrals.** Distinct from 2a. `phi_lookup` averages the conditional uniformly along the branch (an integral of ratios); the alternative weights candidate ages by P(d0 \| t_i), giving a ratio of integrals. | **TESTED AND REJECTED as the bias explanation (T1).** In matched 10 Mb infinite-sites simulations, denominator weighting moved estimates strongly younger and generally increased RMSE. Across the complete 10K/50K/100K sets its bias was -836/-563/-612 generations, versus +248/+423/+632 for uniform; RMSE was 889/653/705 versus 479/613/726. Eight completed 200K replicates agreed (weighted bias -1270, RMSE 1326; uniform +480, 742). The `betabinom` calibration result does not transfer because its weight conditions on k observed at T, whereas this one conditions on d0 observed at the present. Retain uniform as the default. |
| **3** | **Insufficient conditioning on the ARG beyond d0 and the age interval.** The conditional keeps only the *count* of descendants (eq. 4). Two edges with the same d0 can sit in quite different local genealogies, and under the structured coalescent the mutant class coalesces at rate proportional to 1/x, so branch lengths *within* the mutant clade also carry frequency information. That term is dropped. | **CONFIRMED at the level of localisation; the mechanism within edge treatment remains open.** Exact `Mutation.time` gives +10 generation bias, while replacing only that point with its true ARG edge interval gives +364 at Ne=50K. This removes ARG error, allele mapping, variable Ne and eps. T6 proves the marginalisation composition; T1 tested the alternative mutation-time measure. The surviving candidate is information in the tree beyond `(d0, interval)`. Specified as **T5**. |
| 4 | Wrong diffusion conditioning / boundary behaviour. | **TESTED AND REJECTED.** T3 found no predicted-probability error with the magnitude or Ne pattern needed to explain the bias. More decisively, 100 constant-Ne msprime simulations using exact mutation times gave +10 generation MAP bias, -6 generation posterior-mean bias and slope 0.992 across 7.35 million sites. The conditional is calibrated when `t_i` is known; any remaining failure is introduced by representing `t_i` as an edge interval or conditioning on that representation. |
| 5 | Edge quadrature and interpolation. | **IMPLEMENTATION FIXED; contribution to the bias now measured and negligible (T4).** The former 16-node boundary case was 30x high in a constructed regression case, but on real simulated data the switch to knot-split analytic integration moves the estimate by only **-2.7 / -11.3 / -18.9 generations** at Ne = 10K / 50K / 100K. The pathological geometry is rare enough not to matter in aggregate. Worth keeping fixed; not a bias explanation. |
| 6 | Ne scaling / haploid-diploid convention mismatch. | **RETIRED.** A factor-of-two convention error would produce a clean factor-of-two displacement in diffusion time and an error in generations proportional to Ne. The observed offset is roughly generation-scale across Ne and has neither signature. |
| 7 | Modern-polymorphism ascertainment mismatch. | **RETIRED as an explanation for the simulation bias.** The simulated data are generated and analyzed under the same modern-polymorphism ascertainment, yet the bias remains. Ascertainment differences may still matter when transferring the method to real data, but they cannot cause the bias under diagnosis here. |
| 8 | Composite-likelihood dependence. | OPEN, keep last. Dependence inflates precision without biasing calibrated marginals. My block-bootstrap SDs came out *conservative* (947 estimated vs 717 actual scatter), so on that evidence coverage failures here are bias, not underestimated variance. |

Retired in `bias_ideas.md` already and not revisited: T-grid resolution, ARG
inference error (true ARGs used), genotype error (none simulated).

## Open bugs

| Severity | Bug | Notes |
|---|---|---|
| Medium | **tau_i = 3 cutoff recovered 6.3% early across a demographic boundary.** With Ne=10,000 to generation 59,000 then Ne=1,000, the exact cutoff is 59,100 generations; interpolating log-spaced age rows returns 55,352. Edges starting in between are wrongly excluded. | MATH.md admits the approximation. Fix: store the breakpoints or the exactly inverted cutoff with the table. Exact within a single constant-Ne window. |
| Medium | **Build cost is concentrated in rows inference discards.** One (n=26, 300 sample-age) row costs 0.02 s at tau_i=5e-4 and 2.9 s at tau_i=2.33, but **42 s at tau_i=2000**, where 930 digits are needed — and `--mutation-age-max` throws away everything past tau_i=3. | Cap the age grid nearer the cutoff. Deliberately not done: it is a semantics change, not a bug. |
| Low | **Store coordinate convention is assumed, not validated.** The adapter subtracts one because this repo's converter declares one-based ARG coordinates. A store preserving ordinary zero-based tree positions would silently select the neighbouring site. | Not an unconditional bug for current stores. Fix: validate coordinate-convention metadata at ingestion. This project has been bitten by POS conventions before. |

**Not a bug, but the standing engineering option:** reformulating the conditioning
in a numerically stable basis (orthogonal-polynomial or spectral moments rather
than raw power moments) would cut the digits required and could retire the whole
per-row precision apparatus, on both branches. It is **bias-neutral by
construction** — it computes the same expectations with less cancellation, so it
is deliberately *not* a bias hypothesis. Noted in MATH.md; not needed at ARG panel
sizes now that precision is handled.

Recently fixed, for the record: five inference bugs in bc5268d (square chunks
transposing sites against samples; resolver eligibility never applied; p_T = 0
violated below table coverage; exact age knots poisoned by a zero-weight NaN;
invalid moment pairs becoming confident likelihoods) and three numerics bugs in
d75b9b9 / 54b9113 / 33b23d6.

## Judgement calls

Decisions taken that a reader might reasonably want back.

- **`--precision` is a floor, not a cap.** Accuracy is not negotiable and the
  digits a row needs are a property of the row. A user passing `--precision 80` to
  bound cost will not get it. Metadata records the range actually used.
- **Precision is chosen per age row, then verified from the answer.** The formula
  only estimates the loss; the exact moment identities are checked on the computed
  values and a violation triggers recomputation higher. An underestimate costs
  time, never accuracy.
- **Ceiling of 1400 digits, then NaN.** Covers tau_i up to ~3000 at n=26, far past
  the tau_i=3 cutoff. Beyond it, entries are NaN rather than aborting the build.
- **NaN now means two things** — d0 above the panel size, or precision exhausted.
  Inference treats NaN in any draw as disqualifying the whole site, since eq. (11)
  is a mixture over all G draws and one draw cannot be dropped.
- **A constraint violation is treated as a measurement, not a bug.** It says the
  precision fell short. The corollary is that these identities can no longer serve
  as an independent check of correctness.
- **Genotype-code orientation is derived from counts and errors on ambiguity**,
  rather than being guessed from shape. Test fixtures were transposed relative to
  the real reader and now match it.
- **eps = 0.01 remains the real-data default**, but zero-error simulation tests
  use eps=0. A nonzero value has no physical role in those tests.
- **MATH.md eq. (10b) specifies an integral of ratios.** The code implements that
  faithfully and retains it as the default. T1 added the discarded denominator to
  the table as `log_den` and exposed the ratio-of-integrals alternative as
  `--marginalise weighted`; the matched simulation comparison strongly worsened
  the downward bias under weighting, so H2b is retired as an explanation.

  The same question appears in three places, which is worth seeing as one issue:
  here as eq. (10b); on the `betabinom` branch as the prior-versus-posterior
  straddling weight w (MATH2.md eq. 8 against 8a, where the honest reason for the
  prior is that the posterior is not well defined — below T the allele does not
  exist, so there is no binomial factor to weigh against); and in `bias_ideas.md`
  as hypothesis 2. On the `betabinom` branch the *order* question is settled in
  favour of the weighted form and verified; the *weight for the straddling
  fraction* remains approximate there. That result does not transfer to `main`:
  T1 directly tested it here and favours retaining the uniform conditional.

## Confirmed sound

Verified, so they need not be re-litigated: eq. (11)'s draw mixture is composed
correctly (sites accumulated first, then a stabilised log-mean-exp over draws, not
a per-site average); polarity complements sit on the correct branches; eps is
applied consistently to the carried and absent branches; multiply-mapped sites are
excluded; whole-site rejection precedes accumulation, so no partially-accepted
site contributes. The closed-form partial-fraction expansion of e^{B tau} matches
`scipy.linalg.expm` and the independent 80-digit reference to the last digit
tested, at every n and tau_i tried.

Two of those are now measurements rather than assertions. **The within-edge age
marginalisation commutes with the likelihood exactly** — averaging the frequency
along the edge and then forming the site likelihood is the same number as
averaging the likelihood, because the site term is affine in the tabulated
moments (T6). **The draw mixture is doing real work at chromosome scale**: its
effective sample size is 7.8 of 10 draws, and collapsing it to the per-site
average costs a median total variation of 0.025 with 90% intervals 1.9% *wider*,
so the wrong order overstates uncertainty rather than understating it (T7). That
settles a direction REVIEW.md asserted without derivation, in the opposite sense.

The shared terminal likelihood is also confirmed independently of all frequency
approximations: with 10,000 independent SLiM loci and the true p_i(T), all seven
true ages were inside the 95% posterior intervals (MAP MAE 25, RMSE 43
generations). Therefore a residual age bias cannot be attributed to multiplying
the per-site Bernoulli terms, normalising the age posterior, or converting forward
generation to age-before-present. This control does not validate estimated
p_i(T), polarity, coordinate matching, or ARG-to-frequency lookup.

## TODO tests

Specified in enough detail to hand off. Add to this list rather than keeping test
plans in chat.

### T1. Marginalisation order: uniform average vs den-weighted (bias H2b) — COMPLETE

**Question.** MATH.md eq. (10b) averages the conditional uniformly along the
mutation's edge (an integral of ratios). The alternative weights each candidate
age by the probability of the observed panel count, P(d0 | t_i), giving a ratio
of integrals. On the `betabinom` branch the uniform form miscalibrates by up to
11x; this test asked whether that result transfers to `main`.

**Implementation.** den was already computed and discarded, and it does *not*
depend on T -- in `ExactMomentEngine.grid` it is built outside the `for it, tT`
loop, so it is a function of (n, d0, age) alone. Since num = phi * den, the
weighted form is a **den-weighted average of the phi already in the table**:

    p(T) = integral[ den(a) * phi(a,T) da ] / integral[ den(a) da ]

**Prerequisite.** den underflowing to zero was exactly the
silent-zero bug fixed in 33b23d6 (576 of 3120 entries). Weighting by den before
that fix would have been unreliable, which is part of why this is newly
approachable rather than long-neglected.

**Completed changes.**

1. `ExactMomentEngine.grid` returns log den and writes it as an extra table plane
   of shape (n_panel, n_sample, n_age) -- **no T axis**, so about 1/300 the size
   of `table` at the default n_t=300. Store **log den**: it is essentially
   P(d0 | t_i) and will underflow float32 for rare configurations.
2. `phi_lookup` now integrates both alternatives analytically within each
   log-age table segment, split at every age knot. This also fixes the 16-node
   edge-quadrature boundary error (the regression case was previously 30x high).
3. `--marginalise {uniform,weighted}` exposes both from one table; `uniform`
   remains the default.

**No separate w factor is needed.** The existence dilution falls out: for
placements with a <= T, phi = 0 contributes nothing to the numerator while den
stays in the denominator. That is the diffusion analogue of betabinom's explicit
w -- see MATH2.md eq. (8).

**Comparison run.** Both settings used identical SNPs, true ARGs, no genotype
error, and the same likelihood settings in the 10 Mb infinite-sites Ne sweep.
The 10K, 50K, and 100K sets completed all ten replicates. Eight of ten 200K
replicates were sufficient to establish the same direction, so the last two were
stopped; 500K was excluded by design.

| Ne | Uniform bias / RMSE (generations) | Weighted bias / RMSE (generations) |
|---:|---:|---:|
| 10K | +248 / 479 | -836 / 889 |
| 50K | +423 / 613 | -563 / 653 |
| 100K | +632 / 726 | -612 / 705 |
| 200K (8/10) | +480 / 742 | -1270 / 1326 |

**Answer and decision.** Weighting does not merely leave the offset unchanged: it
consistently shifts the posterior younger, creating or worsening downward bias.
H2b is therefore struck as the source of the production bias. Retain the uniform
conditional as the documented default; keep the weighted mode as a diagnostic.

**Two caveats on T1, and the second is unresolved.**

**(i) T1 changed two things at once, so it could not attribute its own
difference.** The same commit replaced the 16-node trapezoid with exact analytic
integration in *both* arms, so neither was the pre-T1 behaviour. Bundling them was
recommended in this document and that was a mistake for attribution. **Resolved by
T4**, which separated them: the quadrature change is worth only -2.7 to -18.9
generations on real data, so the bundling turned out not to matter for the
conclusion.

**(ii) The weighted column reproduces `bias_ideas.md` to within one generation.
RESOLVED — and it is validation, not a problem.**

| Ne | `bias_ideas.md` | T1 "weighted" |
|---:|---:|---:|
| 10,000 | -836 / 890 | -836 / 889 |
| 50,000 | -564 / 653 | -563 / 653 |
| 100,000 | -611 / 705 | -612 / 705 |

The explanation is that `bias_ideas.md` was itself den-weighted: its numbers came
from `/tmp/dnAging-bigsims-estimate.py`, which implemented denominator weighting
independently and never called production `phi_lookup`. Two unrelated
implementations of the same weighted marginalisation — dense trapezoid with
recomputed denominators, versus format5 `log_den` with analytic per-knot
integration — agreeing to within one generation on six quantities is a strong
cross-check of both. See "State of the bias" for the provenance.

This document previously recorded a suspicion that the arms were transposed and
that "H2b struck" might be inverted. That suspicion was **wrong**: T4 ran the
actual pre-session code with uniform marginalisation and got +250 / RMSE 484,
confirming the labels. H2b is struck, not provisional.

### T2. What a nonzero eps costs on error-free data

**Question.** On perfect simulated data the true eps is **zero**: no genotyping
error, true ARGs, and the hard wall never fires (0 carriers in 41,395
opportunities on the `betabinom` branch). So `--epsilon 0.01` is not a
data-quality setting here, it is a **misspecification**, and the question is
simply how much it costs.

**Answer: +151 generations in a synthetic test.** With eps profiled, T_hat came
out at one grid step from truth at every injected error level; with eps wrongly
fixed at 0.01 on error-free data, T_hat moved from 2067 to 2218 against a truth
of 2000 (see D1 for the full table). That is a substantial fraction of the offset
under diagnosis, and it is removed by a one-line change.

**Action: set eps as small as is numerically safe and re-measure.** Do this
before any structural hypothesis, because it is nearly free and it takes H0 and
H1 off the table for simulated data. Keep 1e-6 / 1e-3 / 0.01 as a reported
sensitivity, not a calibration.

**Why the cost is not a uniform shift** — this is the mechanism, and it is why
the effect is a bias rather than a variance loss.

**Measured, at n=26, Ne=20,000, mutation age 4,000 generations** (one point only —
sweep before trusting the magnitudes). All figures in nats.

phi for d0=1 spans 4.9e-3 to 4.3e-2, so eps=0.01 lands *inside* the range of phi
itself, which is what makes the effect asymmetric rather than a uniform shift:

|      | small-T half | (truth) | large-T half | (truth) |
|------|-------------:|--------:|-------------:|--------:|
| d0=1, eps=0.01  | 0.134 | 0.168 | 1.117 | 1.997 |
| d0=1, eps=1e-3  | 0.164 | 0.168 | 1.837 | 1.997 |
| d0=13, eps=0.01 | 0.902 | 0.951 | 2.054 | 2.999 |
| d0=13, eps=1e-3 | 0.946 | 0.951 | 2.845 | 2.999 |

So eps=0.01 keeps 80% of the small-T signal but only 56% of the large-T signal.
It degrades discrimination *among old ages* specifically.

**eps acts almost entirely on carried sites.** Raising it from 1e-6 to 0.01
changes the carried-site signal by 42% (2.191 -> 1.272 at d0=1) and the absent-site
signal by **1%** (0.0397 -> 0.0393). Absent sites also carry 10-70x less
information each (0.04-0.37 nats against 1.2-3.0 carried), so the two channels are
balanced by *counts*, not per-site weight. This is the quantitative content of
bias H1, and it is why H0 and H1 are one hypothesis.

**Deferred to real data.** How eps would be estimated when it is *not* known —
joint profiling, the violation-count check, radiocarbon anchoring, genotype
likelihoods — is worked out and parked under "Deferred: real data" below. None of
it is needed to solve the simulated bias, where the truth is eps = 0.

### T3. Predicted-probability calibration (bias H4) — COMPLETE

At each simulation's true sample age, sites were binned by their **predicted**
ancient carriage probability, not by true frequency. The corrected uniform
lookup is well calibrated across the Ne sweep:

| Ne | all sites: predicted / observed | d0=1: predicted / observed |
|---:|---:|---:|
| 10K | 0.17229 / 0.17584 | 0.02606 / 0.02773 |
| 50K | 0.22656 / 0.22733 | 0.03476 / 0.03473 |
| 100K | 0.23404 / 0.23430 | 0.03559 / 0.03556 |
| 200K | 0.23844 / 0.23873 | 0.03647 / 0.03626 |

The small 10K discrepancy is not repeated across Ne; singleton calibration is
essentially exact at 50K--200K. H4 is retired as the explanation for the
systematic age bias. Full decile results are in
`results/frequency_calibration.tsv`; the reusable harness is
`scripts/check_frequency_calibration.py`.

### T4. Disentangle what actually fixed the bias — COMPLETE

**Result: none of the conditional-side changes explain the bias, and the old
baseline was a different estimator.** Ten replicates per Ne, true ARGs, eps = 0,
identical sites and table throughout; only the code path varies. Site grouping
copied verbatim from `scripts/compare_t1_marginalisation.py` so the numbers are
directly comparable to T1's.

| Ne | orig+float64 | orig+exact | +bug fixes | +analytic | +weighted |
|---:|---:|---:|---:|---:|---:|
| 10,000 | +251 / 484 | +250 / 484 | +250 / 484 | +248 / 479 | -836 / 889 |
| 50,000 | — | +434 / 628 | — | +423 / 613 | -563 / 653 |
| 100,000 | — | +651 / 749 | — | +632 / 726 | -612 / 705 |

Per-step attribution, mean shift in generations:

| Ne | precision | phi_lookup bug fixes | analytic quadrature | den weighting |
|---:|---:|---:|---:|---:|
| 10,000 | -0.2 | +0.0 | -2.7 | **-1083.7** |
| 50,000 | — | — | -11.3 | **-986.5** |
| 100,000 | — | — | -18.9 | **-1243.5** |

**Conclusions.**

1. **Precision, the phi_lookup bug fixes, and the quadrature fix are all
   negligible here** — under 19 generations combined, against a bias of
   250-651. All three are correct and worth keeping; none is the bias. H5 and 5b
   are closed on this evidence.
2. **Den weighting is the only large effect, and it is harmful:** roughly -1,000
   to -1,240 generations, flipping bias negative and worsening RMSE at 10K and
   50K. H2b is struck, now confirmed by an independent implementation.
3. **`bias_ideas.md` was never measuring the production estimator.** Its numbers
   are the den-weighted harness's, reproduced here to within one generation. The
   uniform estimator gives +250 / +434 / +651, not -836 / -564 / -611. See "State
   of the bias".
4. **The bias to explain is therefore positive, growing sub-linearly with Ne, and
   not attributable to anything in the age-marginalisation machinery.** Since eps
   is 0 in these runs and the ARGs are true, the remaining candidates are the
   conditional itself (H3: information in the ARG beyond d0 and the age interval)
   and the sampling model in eq. (4)/(7).

**Caveat.** 200K was started and dropped as expensive and uninformative once the
pattern was established at three Ne values; 500K was never run. The A_float64 and
B_bugfixed arms were only run at Ne=10,000, since both were exactly zero there.

### T5. Does the ARG carry information beyond d0 and the age interval? (bias H3)

The last structural hypothesis standing. Test designs carried over from
`bias_ideas.md` before that file was removed, plus what T4 now implies.

**Why it survived.** Exact-mutation-time msprime inference is now essentially
unbiased, so the diffusion conditional given `(d0, t_i)` is no longer the generic
suspect. The remaining candidate is specifically what changes when exact `t_i`
is replaced by its containing edge: eq. (4) reduces the tree's dependence on x to
the descendant count, while the mutation time is treated as unknown across the
edge. The paired run gives the predicted separation: +10 generations with the
point age and +364 with the true edge interval. The remaining task is to identify
which information discarded by the interval conditional produces that shift.

**One of the three candidates in that transition is now excluded.** Replacing a
point time by an edge admits three distinct failures: the tree carries
information beyond `(d0, interval)`; the measure on candidate times within the
edge is wrong; or the marginalisation is composed in the wrong order relative to
the likelihood. The third is now proven exact (T6), and the second was tested and
rejected as an explanation of the bias (T1, H2b — though note it *is* worth
hundreds of generations, so it is influential without being the culprit). If the
paired edge run is biased, the surviving candidate is the information the edge
retains, which is what tests 2-4 below address.

**Tests, in order of directness.**

1. **Exact-time versus edge-interval inference. COMPLETE.** On the same 100
   Ne=50K simulations, exact mutation times give +10 generation MAP bias and true
   edge intervals give +364. Results are in
   `msprime_exact_time_ne50k/edge_interval_uniform_results/`.
2. **Exact sharing state.** In simulations that contain the ancient lineage,
   compute whether that lineage actually descends from the mutation, and compare
   with the model's p. This is a direct residual, not a calibration curve, and it
   is available because the simulations know the truth.
3. **Matched-edge comparison.** Compare empirical ancient-sharing frequencies
   among edges matched on (d0, age interval) but differing in other tree
   features. Any systematic difference is exactly the information eq. (4)
   discards.
4. **Residual structure.** Test whether residuals cluster by terminal versus
   internal edges, branch length, tree height, or descendant topology. Terminal
   versus internal is the sharpest single split, since a singleton on a terminal
   branch and a singleton on a short internal branch have very different clade
   shapes.

**What would confirm it.** A residual that depends on a tree feature *at fixed*
(d0, age interval), with the sign and magnitude to produce a 250-651 generation
bias growing sub-linearly with Ne. Note the scaling requirement: whatever this is
must reproduce the observed Ne dependence, which is neither generation-constant
nor tau-constant. That is a strong constraint and worth checking against any
candidate before running a full sweep.

**Related, and cheaper:** characterise the bias scaling first (next steps step 3).
If the Ne dependence turns out to match something simple, it may identify the
mechanism without a matched-edge study.

### T6. Does the within-edge marginalisation commute with the likelihood? — COMPLETE

**Question, as raised externally.** Eq. (10b) marginalises the mutation age into
the *frequency* and only then forms the site likelihood. The estimand is the
other order — marginalise the *likelihood* over the age. Written out, the worry
is that the code computes `prod_i integral p(a_i | t_i) p(t_i | edge) dt_i` where
the estimand is `integral prod_i p(a_i | t_i) p(t_i | edge) d(all t_i)`, i.e.
that a product of integrals has been substituted for an integral of a product.

**Answer: the two are identical here, and the reason is worth stating.** The
swap is valid iff the `t_i` are independent under the age prior and each factor
depends only on its own `t_i`. Given one ARG draw the tree is fixed, and under
infinite sites a mutation's placement time is uniform on its own edge
independently across sites — so the integral belongs inside the product, exactly.
The latent that is *not* site-local is the draw index, which fixes every site's
age at once; that one sits outside the site product (T7), and must.

A second, narrower commutation is also load-bearing: the branch integral is
applied to the frequency table rather than to `ell`. That is exact because `ell`
is affine in the tabulated quantities — `qA = eps + (1-2 eps) phi` is linear in
`phi`, and the diploid dosage probabilities are linear in `(phi, phi2)` jointly.
The diploid case only works because `phi_lookup` averages the second-moment
plane over the same edge with the same weights; averaging one plane and not the
other would be wrong, since `P(dosage)` is quadratic in `r`.

**Result** (`tests/test_branch_marginalisation_commutes.py`, 49 tests pass).
Averaging `ell` over point ages along the edge equals `ell` of the averaged
`phi` to `rtol=1e-12` — machine precision, since it is pure algebra — for
haploid carried and absent sites, for all three diploid dosages, and under both
`uniform` and `weighted` measures. Separately, the production knot-splitting
integrator matches an independent 40,001-node trapezoid over production
*point-age* lookups to `rtol=2e-3`; that reference resolves the `T >= t_i`
boundary at node spacing while production resolves it as an exact integration
limit, so it is a genuine check of the integrator rather than a restatement of
it. The low-clip convention is checked too: normalising by the covered width
instead of the true edge width would inflate affected sites, and it does not.

**What this rules out.** The order of the within-edge marginalisation as a bias
candidate, and any future change that quietly breaks affinity in the tabulated
moments — a frequency-dependent eps, a clip applied before the average, or a
likelihood needing a third moment would all now fail these tests.

### T7. Is the draw mixture actually integrating the ARG posterior? — COMPLETE

**Question.** Eq. (11) is composed in the right order, but `sum_i log ell_ig` is
O(n_sites), so the mixture weights `w_g = exp(sum_i log ell_ig - max)` can
collapse onto one draw. If the effective sample size is 1, the answer is
conditional on the modal draw and the correct ordering buys nothing beyond
choosing it. My prior expectation was ESS ~ 1; that was wrong.

**Method, with no rerun required.** `--save-epsilon-data` already stores per-site
per-draw `phi_alt` over the sample-age grid together with the ancient calls, so
the per-draw log-likelihoods can be reconstructed exactly outside the production
run. `scripts/draw_mixture_ess.py` does that and reports both the ESS and the
cost of the swapped order. Applied to the real 10-draw maize output
(`logan_try/results/genome_parts_10draw_eps001`, 10 chromosomes, 15,058 sites, 50
ancient samples, eps = 0.01); outputs in `results/draw_mixture_ess/`.

**Result 1: the mixture does not collapse at this scale.** ESS at the MAP is a
median 7.82 of 10 draws (min 1.16, max 9.97 over chromosome x sample); nothing
falls below 1.05. It does decline with site count, monotonically:

| sites | mean ESS | median | p05 |
|---|---:|---:|---:|
| 1 | 10.00 | 10.00 | 9.99 |
| 100 | 9.70 | 9.93 | 8.69 |
| 300 | 8.81 | 9.41 | 4.97 |
| 1000 | 7.67 | 8.34 | 3.25 |

Whole chromosomes here are 834-2,803 sites. Extrapolating the decline to the
125K-site or genome scale is a projection, not a measurement — the collapse rate
depends on how much the draws actually differ — but the direction is established
and the ESS should be reported per run rather than assumed.

**Result 2: the swapped order is cheap here, and errs the other way.** Replacing
eq. (11) with the per-site draw average (exactly the old behaviour, and for
pseudo-haploid calls an exact algebraic alternative rather than an approximation
of one) gives, across the 50 samples: total variation between posteriors median
0.0245, p95 0.0864, max 0.1211; `|MAP shift|` median 15.0 generations — one grid
cell of 15.04, i.e. nothing — p95 53, max 346; and 90% interval width ratio
swapped/correct median **1.019**, wider in 84% of samples. So the wrong order
**overstates** uncertainty slightly. REVIEW.md asserted it understates it, with
no derivation; the flattening argument was right and the assertion was not.

**Caveat on scope.** These two numbers are real-data measurements at G=10 with
eps = 0.01, not the perfect-simulation setting this document is scoped to. The
script takes any `epsilon_calibration_data.npz`, so pointing it at a simulation
run with `--save-epsilon-data` would give the matched-scope version, and at the
125K-site simarg output would settle the ESS extrapolation directly.

## Next steps, in order

Perfect simulated data only. Real-data work is parked under "Deferred".

The hypothesis list has collapsed. H0, H1, H2a, H2b, H4, H5, H5b, H6 and H7 are
now all closed or retired, and none of them was the bias. What survives:

1. **Decompose edge conditioning.** The paired run is biased: exact `t_i` gives
   +10 generations and the true edge interval gives +364. Compare uniform
   versus denominator-weighted mutation-time measures, exact sharing residuals,
   and matched edges differing in within-clade topology. Do this on these same 100
   simulations before adding another simulation suite.
2. **Extend the re-measured baseline only if needed.** The old uniform baseline is
   +250 / +434 / +651 at Ne = 10K / 50K / 100K. The new exact-time result provides
   a much cleaner zero-bias reference. Widen the candidate-age grid before any
   200K or 500K run.
3. **Characterise the scaling.** The bias is +250 / +434 / +651 across a 10x range
   of Ne: neither constant in generations nor in diffusion time. Identifying what
   it *is* proportional to would point at the mechanism. Worth doing before more
   hypothesis testing, because it is cheap and discriminating.
4. **Watch the variance.** The exact-time run's ordinary 95% coverage is 0.68
   because linked sites are treated as independent. Use a block bootstrap for
   interval calibration; do not confuse this precision problem with point bias.
   In the old baseline, bias/RMSE is 0.52 at Ne=10,000, so ten replicates no
   longer resolve the bias cleanly there. More replicates, or larger regions,
   before drawing fine conclusions at small Ne.
5. **H8 (composite-likelihood dependence) last, as before.** Test designs, carried
   over from `bias_ideas.md`: thin sites by genetic distance and compare point
   estimates; use one mutation per tree or per recombination block; compare
   ordinary posterior intervals against block-bootstrap uncertainty.

**Retired without further work** (carried over from `bias_ideas.md`, which is now
removed): the 20-generation T grid is far too fine to explain errors of hundreds
of generations; the mutation-age cutoff removes few sites in these simulations;
ARG inference error cannot apply because true ARGs are used; and genotype error
cannot apply because these runs have none.

## Deferred: real data

Parked by decision: the estimator must be unbiased on error-free simulations with
true ARGs before any of this matters. Kept because it is worked out and verified,
and because it is the reason eps is not simply "a number we do not know".

### D1. Estimating eps when it is genuinely unknown

**Do not sweep eps. Profile it.** A sweep implies we could then *set* eps to the
right value, and we cannot: the genotyping error rate of real aDNA is not known a
priori. That objection dissolves rather than blocking, because **eps is sharply
identifiable jointly with T**, so it never has to be assumed.

Measured on a synthetic set built from the real conditional (n=26, Ne=20,000,
T_true=2000, 12 age classes x 26 d0 classes, ~15,000 independent sites, carriage
drawn with a KNOWN injected error rate):

| injected eps | eps_hat | T_hat | T error | T_hat with eps fixed at 0.01 | error |
|---:|---:|---:|---:|---:|---:|
| 0      | 1e-6 (grid floor) | 2067 | +67 | 2218 | +218 |
| 1e-4   | 1e-4  | 2067 | +67 | 2067 | +67 |
| 1e-3   | 1e-3  | 2067 | +67 | 2067 | +67 |
| 1e-2   | 1e-2  | 2067 | +67 | 2067 | +67 |

The injected rate is recovered exactly at every level, and the residual +67 is one
T-grid step (151 generations), i.e. essentially unbiased. **Profiling eps makes
T_hat insensitive to it** — 2067 at every injected level — whereas wrongly fixing
eps=0.01 on error-free data costs +151 generations. The profile is peaked, not
ridged: at eps_true=1e-3 the penalty is -56 nats at 1e-6 and -37 nats at 1e-2.

**The real limit is the effective number of INDEPENDENT sites, not the site
count.** Curvature scales with it, and identification fails outright once it is
small:

| independent sites | penalty at eps=1e-6 | penalty at eps=0.01 | eps_hat |
|---:|---:|---:|---:|
| 14,976 | -56.1 | -37.0 | 1e-3 (correct) |
| 2,988  | -20.4 |  -5.4 | 3e-3 |
| 756    |  -6.9 |  -1.4 | 3e-3 |
| 144    |   0.0 |  -0.9 | **1e-6 — identification lost** |
| 24     |   0.0 |  -0.2 | **lost** |

Crossover is around 500-1000 effective independent units. Linked sites are worth
somewhere between one site and their number, so **10 Mb may be borderline while
genome scale is comfortable** — the block bootstrap's 100 blocks over 10 Mb sits
uncomfortably near the row where identification fails. This must be measured on
real linked data before relying on a profiled eps; the table says what to look for
(is the profile peaked, and by how many nats).

**Consequences.**

- Replace "choose eps" with "profile eps jointly with T", and report the profile
  shape as a diagnostic. If it is flat on real data, eps is unknowable *for that
  dataset* and the estimate must be reported as conditional on it.
- **This demotes the split further.** Profiling a scalar eps already removes the
  T-sensitivity, so eps_e/eps_m is only needed if a profiled scalar leaves residual
  bias. Test the scalar profile first.
- **eps_hat is an omnibus failure rate, not a sequencing error rate.** It will
  absorb ARG error, polarity error and model misspecification along with genotyping
  error. That is arguably what the likelihood wants, but it means the fitted value
  must not be reported as, or sanity-checked against, a damage-based error estimate.

**Caveat on all of the above.** This is a well-specified simulation: the model
generating the data is exactly the model being fit, sites are drawn independently,
and phi is exact with no ARG uncertainty. Identifiability under correct
specification does not imply it under misspecification.

**Practical routes, ranked by what they deliver.**

1. **Radiocarbon-dated samples.** Fix T at the known age and fit eps. The only
   route that identifies eps cleanly. Transfers only across comparable library
   prep, coverage and ARG quality.
2. **Joint profile with block-bootstrap uncertainty.** Works, but **never read
   confidence off the profile depth**: under a composite likelihood the score is
   unbiased, so eps_hat stays consistent, while the curvature is inflated by
   pseudo-replication — the sum runs over all N sites while the information is
   worth roughly N/k. Bootstrap the whole (T_hat, eps_hat) fit over genomic
   blocks, which the pipeline already does for T.
3. **The violation-count check — cheap, and can falsify an assumed eps today.**
   Sites whose mutation postdates the sample have p = 0 exactly, so r = eps and
   every observed carriage there is pure failure:
   eps_hat(T) = carried among impossible sites / impossible sites.
   Verified unbiased, but its precision is set by the expected event count
   ~ eps * n_impossible: at eps=1e-3 and 6,240 impossible sites you expect ~6
   events, so ~40% relative precision — a bound, not an estimate, at 10 Mb. It is
   also monotone in T (0.0016 at T=856 rising to 0.0267 at T=6000), so it cannot
   stand alone. **Its value is falsification:** 9 violations among 6,240
   impossible sites is incompatible with eps=0.01, which predicts 62. Even zero
   violations is informative — by the rule of three that bounds eps below
   4.8e-4, ruling out 0.01 outright.
4. **Replicate libraries or duplicated individuals.** Discordance gives the
   genotyping component directly, cleanly separated from ARG error. Gives the
   floor.
5. **Per-site genotype likelihoods — the principled replacement.** Propagate
   per-site, per-allele likelihoods from the BAM (base quality, position-specific
   damage, coverage) instead of a scalar. MATH.md eq. (2)/(3) generalises cleanly:
   replace the symmetric eps-flip with a per-site pair P(obs | derived),
   P(obs | ancestral). This also fixes something a scalar structurally cannot —
   damage is strand- and base-specific, and reference bias is asymmetric again.
6. **Damage/quality models alone — a lower bound only**, since eps_hat is an
   omnibus failure rate absorbing ARG and polarity error. The *gap* is the useful
   part: eps_hat >> eps_damage quantifies ARG or model error.

### D2. The error-term split (former T2)

Splitting eps into eps_e (observation) and eps_m (model failure) with a
T-independent background b_i = d0/n. **Twice demoted:** it gains nothing on signal
(2.021 nats against 2.027 for eps_e alone, because the flat term dominates), and
profiling a scalar eps already removes the T-sensitivity. Its only distinctive
value is making the wall site-dependent — at eps_e=1e-5 a flat 11.51-nat wall
becomes 7.84 at d0=1 against 4.60 at d0=26 — so a violation at a common allele
(likely an ARG error) is cheap while one at a singleton stays expensive. Build it
only if a profiled scalar leaves residual bias on real data.

### D3. Other real-data exposures, not yet worked out

- Reference bias and missing-to-ref: excluded by decision, controlled upstream.
  Note the exposure is not symmetric — it under-calls derived alleles in both the
  panel count and the ancient call, so errors point the same way rather than
  cancelling, and it breaks the symmetric-eps form.
- Contamination.
- **Inferred rather than true ARGs.** The largest untested gap on either branch;
  everything to date used true ARGs.
- Store coordinate-convention validation (see open bugs).
