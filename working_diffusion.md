# Working notes: the diffusion approach (`main`)

Status as of 2026-09-08. Live working document: **status lines are the point**, so
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

The bias is **not** resolved. It is also not yet cleanly characterised, because
every bias measurement to date was taken through a table that we now know carried
percent-level errors and silent zeros. **Re-baselining is the blocking task.**

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

Two datasets disagree about the sign, and this is unresolved.

**(a) `bias_ideas.md`, production pipeline, 26 haplotypes, 10 Mb, true ARG.**
Underestimates, growing at high Ne:

| Ne      | bias  | RMSE | implied SD | bias in tau |
|--------:|------:|-----:|-----------:|------------:|
| 10,000  |  -836 |  890 |        305 |    -0.04180 |
| 50,000  |  -564 |  653 |        329 |    -0.00564 |
| 100,000 |  -611 |  705 |        352 |    -0.00305 |
| 200,000 | -1093 | 1219 |        540 |    -0.00273 |
| 500,000 | -1360 | 1599 |        841 |    -0.00136 |

Two readings not drawn out in `bias_ideas.md` itself:

- **Bias is 88% of RMSE in every row.** A deterministic offset, not noise. Few
  simulations per condition suffice, and the cause is structural or numerical
  rather than stochastic.
- **The bias is near-constant in *generations* (2.4x spread over a 50x range of
  Ne) and wildly non-constant in *diffusion time* (31x).** This points away from
  every drift-scale mechanism: anything living in tau would be flat in tau and
  fan out 50-fold in generations. Generation-scale candidates are the eps floor,
  grid resolution, and the age prior.

**(b) My 300-sim reimplementation on true trees at eps = 1e-6** gives bias **+538**
(50-digit) — opposite sign. Candidate reasons, untested: eps (0.01 vs 1e-6),
different Ne range, and that (b) is a reimplementation of the math rather than the
production pipeline (no normalizeTEs store here, so the eps handling and ARG-draw
mixture in (a) are untested by (b)). **Reconciling these two is a prerequisite for
believing any bias number.**

Contrast with the beta-binomial, whose offset behaves the opposite way: flat in
tau (1.40x across Ne) and fanning 3.3x in generations. Different mechanisms.

## Hypotheses for the bias

Numbering follows `bias_ideas.md`; **H0 is new**, and the priority order below
differs from that file's.

| # | Hypothesis | Status |
|---|---|---|
| **H0** | **The eps floor, acting asymmetrically in T.** r = eps + (1-2 eps) p, and phi *decreases* with T, so the floor bites hardest where phi is smallest — at large T. It therefore compresses the likelihood's ability to discriminate among old ages rather than shifting it uniformly. Measured below: it removes 44% of the large-T signal at d0=1 while leaving 80% of the small-T signal, and it acts almost entirely on **carried** sites (a 42% change there against 1% on absent sites). An eps-driven offset is also ~Ne-independent in generations, matching (a). | **OPEN, highest priority.** Codex measured a 200-generation MAP shift between eps = 1e-6 and 0.01 on a synthetic 2-site fixture. Worth ~1/3 of the offset in the beta-binomial work, the only mechanism of six that survived there. Not in `bias_ideas.md`. Because its effect is concentrated on carried sites, **H0 and H1 are quantitatively the same object** — see TODO test T2. |
| 5b | **Table numerics.** `bias_ideas.md` notes in passing that "a handful of entire very-young moment-table rows were NaN ... at Ne = 500,000". That was the tip of it: float64 discarded 11% of entries at n=12 and 59% at n=40, wrote 18% of a test build as silent zeros, and was already 0.8% wrong where it passed its own guard. | **FIXED** (d75b9b9, 54b9113, 33b23d6). Plausibly the main driver of the high-Ne rows in (a): larger Ne compresses a fixed generation grid into smaller tau, where the intermediate-d0 denominator underflows. **All of (a) must be re-measured.** |
| 1 | **Leverage from rare carried alleles.** log p moves fast when p is small, so a few carried rare sites can outweigh many singleton absences. | OPEN, and probably the same object as H0. Note the *selective* filter test in `bias_ideas.md` cannot distinguish "model wrong about rare carried alleles" from "deleting terms of one sign moves the estimate": removing carried-but-not-absent singletons deletes the log p terms and keeps the log(1-p) terms, so an upward swing is guaranteed. The swings are 3–6x the bias being explained, which is the tell. Run the **symmetric** arms instead. |
| 2a | **Double conditioning on d0.** Once the ARG edge is observed, d0 is determined, so reweighting candidate mutation ages by P(d0 \| t) may condition on the modern count twice. | **DEMOTED** from "leading structural hypothesis". The `betabinom` branch is essentially this fix carried to its limit, and it *loses* at matched precision (calibrated RMSE 1407 vs 974). The exact test in `bias_ideas.md` is cheap and still worth running; the reasoning is sound, but the empirical direction is against it. |
| **2b** | **Marginalisation ORDER over the edge: integral of ratios vs ratio of integrals.** Distinct from 2a. `phi_lookup` averages the conditional uniformly along the branch (an integral of ratios); the alternative weights candidate ages by P(d0 \| t_i), giving a ratio of integrals. | **TESTED AND REJECTED as the bias explanation (T1).** In matched 10 Mb infinite-sites simulations, denominator weighting moved estimates strongly younger and generally increased RMSE. Across the complete 10K/50K/100K sets its bias was -836/-563/-612 generations, versus +248/+423/+632 for uniform; RMSE was 889/653/705 versus 479/613/726. Eight completed 200K replicates agreed (weighted bias -1270, RMSE 1326; uniform +480, 742). The `betabinom` calibration result does not transfer because its weight conditions on k observed at T, whereas this one conditions on d0 observed at the present. Retain uniform as the default. |
| 3 | Insufficient conditioning on the ARG beyond d0 and age interval. | OPEN, untested. |
| 4 | Wrong diffusion conditioning / boundary behaviour. | OPEN. **Caution:** "biased upward at low true frequency" is what an *unbiased* posterior mean does — shrinkage toward the prior. Calibration must be binned on the *predicted* value, E[p_true \| p_hat], not on the true one. Bin on truth and you will reproduce the artifact. I withdrew a claim of my own for this reason. |
| 5 | Edge quadrature and interpolation. | **PARTLY CONFIRMED, OPEN.** See "fixed 16-node quadrature" under open bugs — codex constructed a 30x overestimate. Directionally this pulls toward *younger* ages, matching (a). |
| 6 | Ne scaling / haploid-diploid convention mismatch. | **RETIRED.** A factor-of-two convention error would produce a clean factor-of-two displacement in diffusion time and an error in generations proportional to Ne. The observed offset is roughly generation-scale across Ne and has neither signature. |
| 7 | Modern-polymorphism ascertainment mismatch. | **RETIRED as an explanation for the simulation bias.** The simulated data are generated and analyzed under the same modern-polymorphism ascertainment, yet the bias remains. Ascertainment differences may still matter when transferring the method to real data, but they cannot cause the bias under diagnosis here. |
| 8 | Composite-likelihood dependence. | OPEN, keep last. Dependence inflates precision without biasing calibrated marginals. My block-bootstrap SDs came out *conservative* (947 estimated vs 717 actual scatter), so on that evidence coverage failures here are bias, not underestimated variance. |

Retired in `bias_ideas.md` already and not revisited: T-grid resolution, ARG
inference error (true ARGs used), genotype error (none simulated).

## Open bugs

| Severity | Bug | Notes |
|---|---|---|
| Medium | **Fixed 16-node edge quadrature overestimates near the existence boundary.** Codex's case: age rows [100,1000,10000], T=999, frequencies [0,0.5,0.5], edge [100,1000] gives 0.0166667 against an analytic 0.000555435 — **30x**. Only the endpoint samples the narrow nonzero region and trapezoid gives it half a panel's width. | Fix is to split the integral at each T and table-age knot and integrate the log-linear pieces analytically. Changes the scheme MATH.md documents. Overlaps bias H5 and pulls the right direction to matter. |
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
- **eps = 0.01 remains the default**, unchanged pending H0. This is a judgement
  call by omission and the one I would revisit first.
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

**(i) T1 changed two things at once, so it cannot attribute the difference.** The
same commit replaced the 16-node trapezoid with exact analytic integration split
at every age knot — in *both* arms. Neither arm is therefore the pre-T1
behaviour, and the uniform-vs-weighted contrast is measured on top of a
quadrature change worth up to 30x on the regression case. Bundling the two was
recommended in this document and that was a mistake for attribution.

The more important number is buried by the framing above. Against
`bias_ideas.md`'s -836 / RMSE 890 at Ne=10,000, the current uniform arm gives
**+248 / RMSE 479**: the bias flips sign and shrinks, and RMSE nearly halves.
That is a real improvement, but its cause is **not isolated** — it is the joint
effect of the table precision fixes (d75b9b9, 54b9113, 33b23d6), the five
inference bug fixes (bc5268d), and T1's quadrature change. Disentangling these is
now the most valuable single run available: see T4.

**(ii) The weighted column reproduces `bias_ideas.md` to within one generation,
which is unexplained.**

| Ne | `bias_ideas.md` (pre-T1 production) | T1 "weighted" |
|---:|---:|---:|
| 10,000 | -836 / 890 | -836 / 889 |
| 50,000 | -564 / 653 | -563 / 653 |
| 100,000 | -611 / 705 | -612 / 705 |

Six quantities matching across three conditions. `bias_ideas.md` documents the
pre-existing pipeline, which **could not** compute a den-weighted marginal — the
denominator was not in the table until T1 added it. So the arm reproducing the
old numbers ought to be the *uniform* one. The code's branches read correctly
(`uniform` divides by branch length; `weighted` accumulates den-weighted
numerator and denominator), so mislabelling in the implementation is not the
explanation. Remaining possibilities: the weighted arm ran through a stale table
or code path, the numbers were carried over rather than recomputed, or it is
coincidence. **Until this is resolved, treat "H2b struck" as provisional** — if
the arms are effectively transposed, the same data would instead say weighting
moves -836 to +248 and nearly halves RMSE, i.e. H2b confirmed as a major
contributor. The conclusion is currently inverted by an unexplained match.

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

### T4. Disentangle what actually fixed the bias

**Why.** Between `bias_ideas.md` (-836 / RMSE 890 at Ne=10,000) and T1's uniform
arm (+248 / RMSE 479), three independent sets of changes landed: table precision
(d75b9b9, 54b9113, 33b23d6), five inference bug fixes (bc5268d), and T1's exact
edge integration. The bias flipped sign and RMSE nearly halved, and **we do not
know which change did it.** That is now the most valuable thing to establish,
because it decides where the remaining effort goes.

**Method.** Re-run the Ne sweep at each of these, uniform marginalisation
throughout, on identical SNPs and ARGs:

1. Pre-T1 code with the OLD float64 table — reproduces `bias_ideas.md`, confirming
   the comparison is like-for-like. If it does not reproduce, stop: something else
   differs and the whole sweep is not comparable.
2. Pre-T1 code with the NEW exact table — isolates the precision fixes.
3. Current code (exact table + bug fixes + analytic integration) — the +248 arm.
4. Current code with the 16-node trapezoid restored — isolates T1's quadrature
   change from the bug fixes.

**Also settles caveat (ii) above**, since step 1 establishes which arm reproduces
the old numbers and whether the weighted column's match is real.

**Expected value.** If the precision fixes did it, the bias story is largely over
and H1/H3/H4 can be closed. If the quadrature fix did it, H5 is confirmed and the
remaining offset is whatever survives. Either way the hypothesis list collapses.

## Next steps, in order

Perfect simulated data only. Real-data work is parked under "Deferred".

1. **T4: disentangle what already fixed the bias.** Between `bias_ideas.md`
   (-836 / RMSE 890) and T1's uniform arm (+248 / RMSE 479) the bias flipped sign
   and RMSE nearly halved, and the cause is not isolated. Most of the hypothesis
   list may already be closed. Do this first — it decides where everything else
   goes.
2. **Resolve T1 caveat (ii)**: why the weighted column reproduces
   `bias_ideas.md` to within one generation. Until it is explained, "H2b struck"
   is provisional and could be inverted. T4 step 1 settles it.
3. **T2: set eps small and re-measure.** Nearly free, and takes H0 and H1 off the
   table for simulated data, where the truth is eps = 0.
4. **Re-baseline everything** on the fixed table. Every number in
   `bias_ideas.md` predates the precision and inference fixes.
5. Then whatever offset survives: H3, H4 (binned on *predicted* frequency, not
   true), H8.

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
