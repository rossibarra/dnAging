# Working notes: the beta-binomial approach (`betabinom` branch)

Status as of 2026-09-08. Live working document: **status lines are the point**, so
update them rather than appending. Companion: [working_diffusion.md](working_diffusion.md)
for the approach on `main`.

**Code lives on the `betabinom` branch, not here** — `betabinom/` and `MATH2.md`
are only on that branch. This file sits on `main` so the two approaches can be
read side by side, which is what the branch decision requires. Spec is `MATH2.md`
on that branch (18 commits, all unpushed).

## Where this stands

The likelihood is correct and well validated in its parts, and it was the better
performer for most of this work — but that was measured against a **float64**
implementation of the diffusion, which loses ~0.3 digits per chromosome and so
degrades as the panel grows. At matched precision the diffusion wins:

| calibrated (LOO, in tau, 1/panel covariate) | bias | med abs err | RMSE | 95% cov |
|---|---:|---:|---:|---:|
| beta-binomial | +114 | 758 | 1407 | 0.87 |
| diffusion, 50 digits | +80 | **460** | **974** | 0.83 |

Paired over 300 simulations the beta-binomial wins only 106/300 (mean abs-err
difference +345 +- 51 generations against it, t = 6.8). Uncalibrated it is worse
still: RMSE 2289 against 1182, r 0.874 against 0.946.

**Recommendation: this is not the way forward.** Its value now is as an
independent check on the diffusion — it conditions on the ARG at time T rather
than on the present-day panel count, so the two share no numerics and few
assumptions, and agreement between them is real evidence. Keep the branch, keep
the validation suite, do not build on it.

## The idea, in one paragraph

Given a sample at time T, the probability a pseudo-haploid carries a mutation at
population frequency x is x. Conditional on the ARG at T, a mutation whose child
node is older than T sits above k of n_T lineages, and the frequency posterior is
Beta — so the carrier probability follows directly from the tree rather than from
a diffusion trajectory. A mutation whose parent is younger than T has p = 0. An
edge straddling T is integrated over the placement of the mutation along it,
weighted by the same conditional. Three cases, no matrix exponentials, no
present-day binomial.

## State of the bias

Measured on the 300-simulation benchmark (`~/Projects/mutrates/simbatch300`; Ne
10k–100k, mu 1e-9–1e-8, panel 10–40, age 100–10,000 generations, 10 Mb, true ARGs,
eps = 1e-6):

```
n = 300    bias +1573    median +1404    SD 1666    RMSE 2289
MAP = 1.048 * age + 1318      r = 0.874
```

**The offset is not a constant number of generations.** This corrects `MATH2.md`
section 8, which reports "a constant additive offset of about +1340 generations"
from 20 fixed-parameter simulations. Across varying Ne it is not constant:

| Ne bin              |   n | bias (gen) | bias (tau) |
|---------------------|----:|-----------:|-----------:|
| [1.0e4, 3.1e4)      |  75 |        780 |     0.0188 |
| [3.1e4, 5.4e4)      |  75 |       1137 |     0.0134 |
| [5.4e4, 7.7e4)      |  75 |       1767 |     0.0141 |
| [7.7e4, 1.0e5)      |  75 |       2607 |     0.0145 |
| all                 | 300 |       1573 |     0.0152 |

Spread across Ne: **3.34x in generations, 1.40x in diffusion time.** In
generations Ne is overwhelmingly significant (+271 per 10,000, t = 8.24); in tau
that dependence vanishes (t = -1.29) and R2 falls from 0.201 to 0.031.

**So the offset is ~0.0152 of coalescent time — a fixed amount of drift, not a
fixed number of generations.** That is a better-behaved artifact than it looked,
and correctable given Ne, which the pipeline already assumes known. Applying the
correction in tau rather than in generations: LOO RMSE 1520 against 1669, median
abs err 857 against 971.

A small residual panel-size effect survives (-25.8 generations per haplotype,
t = -2.64; -0.00023 in tau, t = -2.21).

Note this is the **opposite** behaviour to the bias in `bias_ideas.md` for the
diffusion pipeline, which is near-constant in generations (2.4x) and fans 31x in
tau. Different mechanisms; do not pool the two.

## Hypotheses for the offset

| Hypothesis | Status |
|---|---|
| **eps floor.** For a rare carried site, r = eps + (1-2 eps) p is dominated by eps. | **CONFIRMED contributor, ~1/3 of the offset.** The only survivor of the candidates tested. Raising eps from 1e-6 to 1e-3 on error-free data moved RMSE 1169 to 1562 and bias +886 to higher. |
| Clade shape | **STRUCK** by test. |
| The straddling-edge term | **STRUCK** by test. |
| A second binomial on clade size (the "phi2 / nleaf" idea) | **STRUCK, and it was wrong.** It multiplied two binomials for the same sample; observed/predicted degraded from 1.002 to 1.079. |
| Age-marginalisation order (ratio of integrals vs integral of ratios) | **STRUCK *for this branch*, and the evidence matters for `main`.** This branch uses the den-weighted ratio-of-integrals form, which is what makes calibration flat; the uniform average gives observed/predicted of 0.09, 0.17, 0.31 across bins of the fraction of edge above T against 1.12, 1.01, 0.92 with it, i.e. up to **11x** miscalibration, and inflates p by up to 2x on long edges. Struck here means "we do it right and checked", **not** "it does not matter" — `main` uses the uniform average, so this is live evidence for H2b in [working_diffusion.md](working_diffusion.md) — specified there as TODO test T1 — not a closed question. |
| Panel ascertainment | **STRUCK** by test. |
| Multiallelic exclusion | **STRUCK** by test. |
| Residual: whatever remains of the ~0.015 tau offset after eps | **OPEN.** Unexplained. Note the diffusion at 50 digits carries a same-signed tau offset of ~0.005 (2.12x spread across Ne), so part of this may be common to both likelihoods rather than specific to this one. |

## What is validated

Established against msprime with true ARGs:

- **p = 0 for mutations postdating the sample is exact**: 0 carriers in 41,395
  opportunities.
- **The large-age limit k/(n_T + 1) is accurate to 0.15%** for edges above T. The
  unpolarised alternative k/n_T is wrong by 9.6% — the +1 matters and comes from
  the polarised 1/x prior.
- **Calibration is flat across strata**: by predicted p, by fraction of edge above
  T, by clade size, and on a log scale down to p ~ 1e-3. Observed/predicted 0.996
  for straddling edges, 0.999 for edges above T; worst single bin 1.2 SD.
- **The den-weighting is what makes that true.** Averaging uniformly instead gives
  ratios of 0.09, 0.17, 0.31 across bins of the fraction of edge above T, against
  1.12, 1.01, 0.92 with it.
- **Slope 1.002** (correlation 0.960) over 20 fixed-parameter simulations at
  Ne = 50,000 with true ages 454–9,293 generations.
- **Blind test** on `bigsims`: 8,709 against a true 8,587, +1.4%.
- **Pseudohaploid matches single-haplotype** MAPs to within one grid step (252
  generations mean abs difference, 4.7% of mean true age).

## Judgement calls

- **Polarised prior**, giving posterior Beta(k, n_T - k + 1) and mean k/(n_T + 1)
  in the large-age limit. The unpolarised alternative gives k/n_T and is wrong by
  9.6%, so this is validated rather than merely chosen. Note this is *not* a
  polarity decision: see below.
- **Two approximations are knowingly retained** (MATH2.md sections 3 and 4):
  dropping the structured-coalescent term, and using the prior rather than the
  posterior in the straddling weight w. Tests cannot resolve either, which is not
  the same as their being exact — and this likelihood has already proved sensitive
  to sub-1% structure, so neither should be assumed harmless.
- **mpmath at dps 45–60** throughout, because the unified conditional's alternating
  sums lose ~0.3n digits. This was the right call and is precisely what the
  diffusion side lacked until today — the comparison was unfair for that reason.
- **eps = 1e-6** for all benchmarking here, against MATH.md's default of 0.01.
  Any cross-comparison of the two approaches must match eps first.
- **Composite likelihood over linked sites, with a block bootstrap** over
  contiguous genomic windows instead of analytic intervals. Block-bootstrap SDs
  came out slightly **conservative** (947 estimated against 717 actual scatter),
  so every coverage failure here is bias, not underestimated variance.
- **Pseudohaploid is treated as one Bernoulli draw** per site, so the likelihood
  uses only the first moment. Reference bias and missing-to-ref are explicitly out
  of scope, controlled upstream.
- The pseudohaploid SD is **0.66x** the single-haplotype SD (860 to 568
  generations). That tightening is not a gain: the offset is unchanged, so
  bias/SD rises from ~1.0 to ~1.8 and coverage gets *worse*. Anything that cuts
  variance without addressing the offset makes the estimator more confidently
  wrong.

## Polarity: this approach largely does not need it

Worth stating explicitly, because it is a genuine structural difference from
`main` and easy to get backwards.

**Polarity is implied by the ARG, not supplied to the model.** A mutation is an
event on an edge, so "derived" means precisely "descends from that edge's child
node", and k is read off the topology:

```python
m = site.mutations[0]                        # the mutation's edge
nl = tree.num_samples(m.node)                # leaves below it
kv = nl - searchsorted(cache[m.node], GRID)  # lineages at T descending from m.node
```

Nothing in `betabinom/` reads REF/ALT or an ancestral state — there is no polarity
table and no ALT-ancestral branch. The `len(site.mutations) != 1` filter drops
recurrent and back mutations, so the identification is exact. The prior is not
literally f(x) ∝ 1/x either: the moments come from the neutral recursion started
at a single copy a generations before T, which is polarised by construction; 1/x
and k/(n_T+1) are only its large-a limit.

Contrast `main`, which conditions on d0, a **count from the VCF**. It must know
whether ALT is derived to index the conditional, and polarity can flip between ARG
draws — hence MATH.md's polarity table, per-draw ALT-ancestral transform, and the
d0 -> n - d0 complement. That machinery exists because the conditioning variable
is a called allele count rather than a tree.

Polarity still reaches this approach in two places, both **outside** the model:

1. **The ancient sample's call.** `f[9].startswith('1')` means "carries ALT", so
   `carried` assumes ALT is derived — true by msprime convention, but needing an
   ancestral state on real data. A flip inverts the *observation*, not the model:
   a more localised failure than picking the wrong row of a conditional.
2. **Upstream ARG inference.** With inferred rather than true ARGs, polarised data
   goes into building the ARG, so mispolarisation would corrupt mutation placement
   and topology. The exposure moves upstream rather than vanishing, and is
   untested here — see open question 2.

**Consequence for the branch decision.** The two approaches fail *differently*
under polarity error, so agreement between them is evidence against polarity error
as a cause of any shared offset. This is part of why the branch is worth keeping as
an independent check even though it is not the way forward.

## Struck, with reasons — do not revisit

- **phi2 / nleaf double-binomial conditioning.** Multiplied two binomials for the
  same sample. Degraded observed/predicted 1.002 to 1.079.
- **Private-mutation clock.** Inapplicable: SNPs in the ancient sample are
  ascertained as present in some modern sample. Verified concretely — 1,166 of
  18,714 carried sites are in neither the modern VCF nor the tree.
- **"Calibration is unstable between datasets."** Withdrawn; it was a post-hoc
  comparison, and pooling 20 simulations gives slope 1.002.
- **"The frequency channel is nearly empty."** Withdrawn; 0.007 nats/site read as
  no signal, but x 15,489 sites is ~108 nats.
- **"Rare-end over-prediction causes the age bias."** Withdrawn; shrinkage toward
  the prior is expected of a posterior mean, not evidence of miscalibration.
  Calibration must be binned on the predicted value, not the true one.

## Open questions

1. What is the residual tau offset after eps is accounted for, and is it shared
   with the diffusion? Both carry same-signed tau offsets (~0.015 here, ~0.005
   there), which would suggest a common cause rather than two coincidences.
2. Does the correction hold under **inferred** rather than true ARGs? Everything
   here used true ARGs. This is the largest untested gap.
3. Does it hold under non-constant Ne? All benchmarking was constant-Ne.
4. The residual panel-size effect (t = -2.6) is unexplained.

## If this branch is revived

Two things would have to happen first: apply the offset correction **in diffusion
time** rather than in generations (a constant-generation correction is invalid, see
above), and re-run the head-to-head against the 50-digit diffusion rather than the
float64 one, since that is the comparison that actually decided against it.
