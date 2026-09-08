# Working notes: the diffusion approach (`main`)

Status as of 2026-09-08. Live working document: **status lines are the point**, so
update them rather than appending. Companion: [working_betabinom.md](working_betabinom.md)
for the alternative likelihood. Spec is [MATH.md](MATH.md); code is
[precompute_freq_trajectory_moments.py](precompute_freq_trajectory_moments.py)
(table build) and [posterior_sample_age_infer.py](posterior_sample_age_infer.py)
(inference).

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
| **H0** | **The eps floor.** `--epsilon` defaults to 0.01 and enters as r = eps + (1-2 eps) p, so every site with true p < 0.01 is modelled at ~0.01 — exactly the rare carried sites H1 identifies as high-leverage. An eps-driven offset is ~Ne-independent in generations, matching (a). | **OPEN, highest priority.** Codex measured a 200-generation MAP shift between eps = 1e-6 and 0.01 on a synthetic 2-site fixture. Worth ~1/3 of the offset in the beta-binomial work, the only mechanism of six that survived there. Not in `bias_ideas.md`. |
| 5b | **Table numerics.** `bias_ideas.md` notes in passing that "a handful of entire very-young moment-table rows were NaN ... at Ne = 500,000". That was the tip of it: float64 discarded 11% of entries at n=12 and 59% at n=40, wrote 18% of a test build as silent zeros, and was already 0.8% wrong where it passed its own guard. | **FIXED** (d75b9b9, 54b9113, 33b23d6). Plausibly the main driver of the high-Ne rows in (a): larger Ne compresses a fixed generation grid into smaller tau, where the intermediate-d0 denominator underflows. **All of (a) must be re-measured.** |
| 1 | **Leverage from rare carried alleles.** log p moves fast when p is small, so a few carried rare sites can outweigh many singleton absences. | OPEN, and probably the same object as H0. Note the *selective* filter test in `bias_ideas.md` cannot distinguish "model wrong about rare carried alleles" from "deleting terms of one sign moves the estimate": removing carried-but-not-absent singletons deletes the log p terms and keeps the log(1-p) terms, so an upward swing is guaranteed. The swings are 3–6x the bias being explained, which is the tell. Run the **symmetric** arms instead. |
| 2a | **Double conditioning on d0.** Once the ARG edge is observed, d0 is determined, so reweighting candidate mutation ages by P(d0 \| t) may condition on the modern count twice. | **DEMOTED** from "leading structural hypothesis". The `betabinom` branch is essentially this fix carried to its limit, and it *loses* at matched precision (calibrated RMSE 1407 vs 974). The exact test in `bias_ideas.md` is cheap and still worth running; the reasoning is sound, but the empirical direction is against it. |
| **2b** | **Marginalisation ORDER over the edge: integral of ratios vs ratio of integrals.** Distinct from 2a, and much better supported. `phi_lookup` averages the conditional uniformly along the branch (`np.linspace(lo, hi, n_quad)`, trapezoidal), which is an integral of ratios. The weighted form needs P(d0 \| t_i), which is not in the table. | **OPEN, second priority after H0.** The `betabinom` branch **tested exactly this question** and the effect is large: uniform averaging gives observed/predicted of 0.09, 0.17, 0.31 across bins of the fraction of edge above T, against 1.12, 1.01, 0.92 den-weighted — up to **11x** miscalibration, and up to 2x inflation of p on long edges (MATH2.md sections 4 and 8). `main` currently uses the form that failed there. Caveats: the weights are analogous but not identical (k observed *at* T versus d0 observed at the present, much later), and a uniform inflation of p would push *older* whereas the observed bias is younger — so if this matters it acts through how the inflation varies with T, not its level. |
| 3 | Insufficient conditioning on the ARG beyond d0 and age interval. | OPEN, untested. |
| 4 | Wrong diffusion conditioning / boundary behaviour. | OPEN. **Caution:** "biased upward at low true frequency" is what an *unbiased* posterior mean does — shrinkage toward the prior. Calibration must be binned on the *predicted* value, E[p_true \| p_hat], not on the true one. Bin on truth and you will reproduce the artifact. I withdrew a claim of my own for this reason. |
| 5 | Edge quadrature and interpolation. | **PARTLY CONFIRMED, OPEN.** See "fixed 16-node quadrature" under open bugs — codex constructed a 30x overestimate. Directionally this pulls toward *younger* ages, matching (a). |
| 6 | Ne scaling / haploid-diploid convention mismatch. | **Cheapest to retire; recommend retiring.** A factor-of-two error gives a clean 2x in tau, not this pattern. Checkable analytically. |
| 7 | Modern-polymorphism ascertainment mismatch. | OPEN, overlaps H2. |
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
- **MATH.md eq. (10b) specifies an integral of ratios**, and the code implements
  that faithfully — it cannot do otherwise, since the denominator is not in the
  table. Whether it is *right* is a modelling question: uniform-along-edge is
  correct as the conditional age distribution given the draw, but wrong as a prior
  that still needs the sampling weight applied. **This is bias H2b and needs a
  decision, not a patch**; storing the numerator and denominator separately is the
  precondition for even testing the alternative.

  The same question appears in three places, which is worth seeing as one issue:
  here as eq. (10b); on the `betabinom` branch as the prior-versus-posterior
  straddling weight w (MATH2.md eq. 8 against 8a, where the honest reason for the
  prior is that the posterior is not well defined — below T the allele does not
  exist, so there is no binomial factor to weigh against); and in `bias_ideas.md`
  as hypothesis 2. On the `betabinom` branch the *order* question is settled in
  favour of the weighted form and verified; the *weight for the straddling
  fraction* remains approximate there. Here neither is settled.

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

### T1. Marginalisation order: uniform average vs den-weighted (bias H2b)

**Question.** MATH.md eq. (10b) averages the conditional uniformly along the
mutation's edge (an integral of ratios). The alternative weights each candidate
age by the probability of the observed panel count, P(d0 | t_i), giving a ratio
of integrals. On the `betabinom` branch the uniform form miscalibrates by up to
11x; `main` uses the uniform form and has never been able to compute the other.

**Why it is now cheap.** den is already computed and discarded, and it does *not*
depend on T -- in `ExactMomentEngine.grid` it is built outside the `for it, tT`
loop, so it is a function of (n, d0, age) alone. Since num = phi * den, the
weighted form is a **den-weighted average of the phi already in the table**:

    p(T) = integral[ den(a) * phi(a,T) da ] / integral[ den(a) da ]

**Prerequisite, already met.** den underflowing to zero was exactly the
silent-zero bug fixed in 33b23d6 (576 of 3120 entries). Weighting by den before
that fix would have been unreliable, which is part of why this is newly
approachable rather than long-neglected.

**Changes.**

1. Return den from `ExactMomentEngine.grid` and write it as an extra table plane
   of shape (n_panel, n_sample, n_age) -- **no T axis**, so about 1/300 the size
   of `table` at the default n_t=300. Store **log den**: it is essentially
   P(d0 | t_i) and will underflow float32 for rare configurations.
2. In `phi_lookup`, multiply the trapezoidal node weights by den at each node.
   Do this **together with** the quadrature fix (open bugs, and bias H5) -- both
   live in the same few lines, and measuring one through the other's error would
   waste the run.
3. Expose `--marginalise {uniform,weighted}` so both are runnable from one table.
   Commits to neither reading of eq. (10b).

**No separate w factor is needed.** The existence dilution falls out: for
placements with a <= T, phi = 0 contributes nothing to the numerator while den
stays in the denominator. That is the diffusion analogue of betabinom's explicit
w -- see MATH2.md eq. (8).

**Comparison.** Both settings over `~/Projects/mutrates/simbatch300` (harness
exists from the head-to-head). Report bias, RMSE, slope, r, coverage, and
parameter dependence in generations *and* in diffusion time, since the two
approaches' offsets behave oppositely under that reparameterisation.

**What counts as an answer.** The weighted form reducing the offset is the
hypothesis. Two outcomes are informative and one is a trap:

- Offset shrinks materially -> H2b confirmed, and eq. (10b) needs rewriting.
- Offset unchanged -> H2b struck for `main`, and the betabinom calibration result
  does not transfer. Record *why*: the weights are analogous but not identical
  (k observed at T against d0 observed at the present, much later).
- **Trap:** a uniform inflation of p pushes estimates *older*, while the observed
  bias is *younger*. Do not read "wrong direction" as "no effect" -- if this
  matters it acts through how the inflation varies with T, not through its level.
  Compare the whole p(T) curve, not just the MAP.

**Decision left open.** Which form becomes the default in MATH.md is a modelling
call, not something the run settles by itself: uniform-along-edge is defensible as
the conditional age distribution given the draw, and indefensible as a prior still
awaiting its sampling weight.

## Next steps, in order

1. **Rebuild the tables with the fixed engine and re-measure the `bias_ideas.md`
   table.** Blocks everything else; every existing bias number came through a
   contaminated table.
2. **Reconcile the sign disagreement** between the production pipeline and the
   300-sim reimplementation. Most likely eps.
3. **H0: eps sweep** at 1e-6, 1e-3, 0.01 on the same simulations.
4. **H2b: the marginalisation order — see TODO test T1**, which is specified.
   Currently untestable rather than tested, and the fix is small: den is already
   computed and thrown away. Bundle with the quadrature fix below.
5. Fix the edge quadrature (bias H5, and an open bug in its own right). Same few
   lines as T1; do them together.
6. H1's symmetric arms; likelihood pull by d0 and carried/absent state.
7. H4 calibration, binned on predicted rather than true frequency.
8. Retire H6 analytically. Then H3, H7, H8.
