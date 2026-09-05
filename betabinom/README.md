# betabinom: an ARG-conditioned likelihood for ancient-sample age

Prototype and validation for an alternative to the diffusion likelihood in
[MATH.md](../MATH.md). Nothing here is wired into the pipeline; it is simulation
work establishing whether the reformulation is correct and whether it recovers a
known age.

## The model

For a pseudohaploid ancient sample at time `T`, the probability it carries a
mutation is the population frequency `X(T)`. Read two numbers off the
**modern-only** marginal tree at `T`:

* `n_T` — ancestral lineages of the panel at `T`
* `k`   — how many of them descend from the mutation

Given `X(T)=x`, the `n_T` lineages are an exchangeable sample from the population
at `T`, so `P(k of n_T | x) ∝ x^k (1-x)^(n_T-k)`. Conditioning on the tree makes
the present-day panel count `d_0` redundant — the tree already determines it — so
the only prior needed is the frequency density of a mutation of age `a`, whose
moments `M_m` the neutral recursion supplies. Expanding `(1-x)^(n_T-k)`:

    phi(k, n_T, a) = sum_j C(n_T-k, j) (-1)^j M_{k+1+j}(tau_a)
                   / sum_j C(n_T-k, j) (-1)^j M_{k+j}(tau_a)

One formula covers every case. A mutation whose edge lies entirely above `T` has
`k >= 1` and a large age; a mutation on an edge straddling `T` has `k = 1` and a
small age; a mutation whose edge is entirely below `T` gets `p = 0`. Mutation
position is uniform on its edge, so `p` integrates `phi` over the part of the edge
above `T`, normalised by the full edge length.

This is cheaper than the diffusion route: moments to order `n_T+1` (about the same
size as the current table) but **one** matrix exponential per age instead of three
per entry, and no two-time bridge — conditioning at `T` on the tree replaces
conditioning at the present on `d_0`.

## What simulation established

Verified against msprime (`Ne = 100,000`, 26 modern + 1 ancient, true ARG):

* **`p = 0` for mutations younger than the sample is exact** — 0 carriers in
  41,395 opportunities.
* **`k/(n_T+1)` is right for edges above `T`**, to 0.15%. It is the
  *age-marginalised* limit of `phi`; the `tau -> inf` limit is `(k+1)/(n_T+2)`.
* **`k/n_T` is wrong** by 9.6% — that is the `Beta(a1,a2)` (unpolarised) prior.
  The correct posterior is `Beta(a1, a2+1)`, i.e. prior `f(x) ∝ 1/x`.
* **Uniform placement alone is not enough for straddling edges.** A mutation born
  just above `T` sits at frequency `~1/(2Ne)`, not `1/(n_T+1)`. Before `phi`, the
  observed/predicted ratio across bins of "fraction of edge above T" ran
  0.00 / 0.09 / 0.17 / 0.31 / 1.07; with `phi` it is flat at ~1.

## What still fails

The age posterior is biased old and does not cover truth (`post3.py`):

    seed  T_true   sites    MAP         95% CI  covers
      11    2500  15,489   7658   7155-7798         no
      22    2500  15,393   3631   2920-3770         no
      33    2500  14,991   4376   4230-4580         no
      44    1000  15,182   1095    953-1240        YES
      55    6000  15,410   8254   7980-8403         no

Findings about why:

1. **The data is informative; the model is not accurate enough.** The expected
   log-likelihood penalty for the wrong `T`, if the model were exact, is ~100 nats
   at `T = 7658` (`power.py`). Misspecification supplies a comparable gain in the
   same direction, so it outbids the signal.
2. **The composite likelihood understates uncertainty by ~30x.** Its intervals
   cover the truth **0 times in 12 seeds**. Resampling 100 contiguous genomic
   blocks (`bootstrap.py`) widens them to a median of 5,334 generations and lifts
   coverage to 8/12. No interval from the raw composite curvature should be
   reported.
3. **Clade shape is NOT the missing information.** The natural hypothesis was that
   `(k, n_T, a)` discards the mutant clade's coalescent history -- fast coalescence
   within the clade implying a lower frequency. Tested directly on 862,343
   observations (`clade.py`), stratified within `(k, n_T)` cells and binned by
   leaves-per-lineage: obs/pred is flat (1.001, 0.995, 1.012, 1.038, 0.990) with no
   monotone trend. The hypothesis is wrong; `(k, n_T, a)` is closer to sufficient
   than expected, and the residual lies elsewhere.
4. **There is a real upward bias, and more sequence does not remove it.** Across 12
   seeds at 2 Mb the MAP median is 4,376 against a truth of 2,500, with 75% of
   estimates above truth. Going to 8 Mb tightens nothing toward the truth -- the
   median moves to 6,017 with 4/4 above (`converge.py`). So the residual is not
   only sampling noise: bias survives a 4x increase in data, and since variance
   shrinks as `1/sqrt(L)` while bias does not, a genome-scale run would converge
   confidently on the wrong age.

### The bias source, located

The likelihood *formula* is correct (Bernoulli in `r = eps + (1-2eps) p`; the tree
factor is common to all `T` and cancels), and the *numerics* are correct
(`numcheck.py`: production quadrature and table interpolation match a 400-node,
240-age reference to 0.2% even at the smallest `p`). The error is in `phi` itself,
in two places the earlier calibration checks could not see:

* **`phi` ignores `nleaf` on straddling edges** (`straddle_nleaf.py`). A mutant
  lineage at `T` that founds a four-leaf clade below `T` sits at a higher frequency
  than one founding a singleton, but both are assigned the same `p`. obs/pred rises
  monotonically with clade size: 0.993 (nleaf=1), 1.046 (2), 1.167 (3-4), 1.182
  (5-8). `clade.py` missed this because it skipped straddling sites, which is
  exactly where small `p` lives.
* **26% over-prediction at `p` in [0.003, 0.01]** (`smallp.py`), about 4 SE and
  worth ~47 excess predicted events against a ~100 nat total signal. The old
  `[0, 0.2)` bin lumped `p=0.003` together with `p=0.19` and averaged it away.

### The nleaf fix was attempted and does not work

`phi2.py` conditions a straddling site on both "1 of n_T carries it at T" and
"nleaf of n carry it today", by multiplying the two binomial factors. This is
**wrong**: the n_T lineages at T are the ancestors of those same 26 samples, so
the two factors are the same observation counted twice. It inflates p uniformly
(1.45x at nleaf=1, where phi was already right) and makes calibration worse --
`phi2_check.py` gives obs/pred 1.079 against phi's 1.002, and 1.127 at nleaf=1
against phi's 1.010.

`nleaf` is not a second sample. It records how fast the mutant lineage branched
between T and the present, and under the structured coalescent that rate is
proportional to 1/X(u) across the whole interval (0,T). Using it correctly needs
the frequency trajectory, not just X(T), which is a much larger change than an
extra table dimension.

The lesson for any further work here: aggregate and coarse-binned calibration is
not evidence. Both real defects were invisible until the bins were made fine on a
log scale and the straddling case was tested on its own.

## Independent test on supplied simulations

Ten true-ARG simulations (`run_archive.py`), independently generated: Ne = 50,000,
10 Mb, 26 modern haplotypes, one diploid ancient individual, true ages 454-9,293.
Ne and the true ages are read from the tree-sequence provenance and sample node
times; `simplify()` to the modern samples strips the ancient lineage, so the true
age never enters the likelihood. One ancient haplotype is used, since that is the
pseudohaploid case that has been validated.

    n=10   Pearson r = 0.970
    bias   mean +1062   median +713   9/10 above truth
    regression:  MAP = 1.032 * T_true + 892
    bootstrap CI width: median 3,466 generations
    coverage: composite 1/10   block-bootstrap 8/10

Three points:

* **The bias is essentially additive** -- slope 1.032, intercept +892. The model
  tracks `T` correctly and sits about 900 generations too old, rather than
  mis-scaling time. This is why the Ne = 100,000 runs above look so much worse: at
  a truth of 2,500 an offset that size is the whole answer, while at T = 8,000 it
  is a 10% error.
* **This regime carries more signal.** tau_T reaches 0.093 here against 0.0125 in
  the Ne = 100,000 runs, and 10 Mb gives five times the blocks.
* **The composite interval is unusable, now confirmed on independent data** --
  1/10 coverage here, 0/12 there. The block bootstrap gives 8/10, about right for
  a nominal 95% given the residual bias, but the intervals are thousands of
  generations wide.

### The additive bias is not the straddling term

`no_straddle.py` retains only sites whose edge never straddles any `T` in the grid,
so the observation set is identical at every `T`. Dropping the straddling term
makes the bias worse, not better:

                       slope   intercept       r   mean bias   coverage
      full model       1.032        +892   0.970       +1062       8/10
      no straddling    1.100       +2187   0.944       +2725       5/10

The test is confounded -- the filter removes 60% of sites and not at random, since
it strips every recent edge -- but the direction settles it: a term responsible for
a +892 bias cannot produce +2187 when deleted.

This matches the earlier per-case decomposition. Above the true age the `old`
channel runs obs/pred > 1 (under-predicting, 1.008 -> 1.031) while straddling runs
< 1 (over-predicting, ~0.91). They push in **opposite** directions, so removing
straddling leaves the `old` drift unopposed. The additive bias therefore lives in
`phi(k, n_T, a)` for edges *above* `T` -- the case measured at 0.15% and considered
settled -- and its slow monotone drift with `T` is why the displacement is roughly
constant rather than proportional.

### The mutation age was being marginalised the wrong way

The code evaluated `phi` at each quadrature node and averaged the *ratios*
uniformly over the edge. That is wrong: the age posterior is not uniform once the
ARG is conditioned on. Writing num(a) and den(a) for the two alternating sums,

    f(a | ARG)  proportional to  f(a) * INT f_a(x) x^k (1-x)^(n_T-k) dx  =  f(a) den(a)

so the marginal is INT num / INT den -- a ratio of integrals, not an integral of
ratios. `den(a)` is P(k of n_T | a); it is near-flat while the allele is rare and
decays as the allele drifts to loss or fixation, so it up-weights young ages. The
error is large and grows with edge span: the den-weighted value is 0.99 of the
uniform average for a 200-generation edge but 0.44 for a 100,000-generation one.

`phid.py` exposes num and den; `run_corrected.py` uses them. It is a genuine
correctness fix and it changes the answer almost not at all:

    CORRECTED  MAP = 1.023*T +1117   r=0.962   mean bias +1239   coverage 8/10
    previous   MAP = 1.032*T  +892   r=0.970   mean bias +1062   coverage 8/10

Six of ten MAPs are bit-identical. The reason is that the correction is a smooth
function of edge span applied at every `T` alike, so it shifts the level of `p`
without much changing the likelihood's *gradient* in `T`, and only the gradient
sets the peak.

### eps matters, and it was misspecified throughout

The supplied simulations have no genotype error -- the ancient sample never carries
a mutation that postdates it (0 violations in 75,827 sites) -- so eps > 0 can only
ever cost likelihood. Every run in this branch before now used eps = 1e-3.

eps is a floor on r, so it converts the hard constraint "the sample cannot carry a
mutation younger than itself" into a finite penalty of log(eps) per violation. The
optimiser can then buy a larger T by paying that finite cost, which is why it has a
T-gradient and the age-marginalisation fix does not.

Same model, only eps changed (`run_corrected.py` takes eps as argv[2]):

                    RMSE    MAE    bias   slope   intercept      r   coverage
      eps=1e-3      1562   1244   +1239   1.023       +1117   0.962      8/10
      eps=1e-6      1169    892    +886   1.026        +748   0.975      8/10

A 25% reduction in RMSE and 28% in bias from the error rate alone. `plot_archive.py`
draws `age_vs_truth.png`: both settings sit above the diagonal at every true age
with slope ~1.02, so the residual is a near-constant offset rather than a timescale
error, and eps accounts for about a third of it.

Note for the real application: on aDNA eps cannot be set to 1e-6, since damage is
real. The gain seen here is only available if eps is estimated rather than fixed --
MATH.md currently fixes it at 0.01 by fiat.

### Blind test on a second set (Archive2)

Ten more true-ARG simulations, set up the realistic way: the `.trees` holds only
the 26 modern haplotypes, the ancient sample arrives as a separate haploid VCF
listing carried sites, and nothing records Ne or the true ages. VCF POS is
`int(tree position) + 1`. Ne was recovered from Watterson at 49,170-50,095
assuming mu = 1e-8, later confirmed as 50,000.

Predictions were made before seeing any true age (`run_archive2.py`), then scored
(`plot_archive2.py` -> `archive2_blind.png`, truths in `ages_archive2.tsv`):

                            RMSE    MAE    bias     SD      r   coverage
      raw MAP               1940   1817   +1817    717  0.963      5/10
      Archive1-corrected    1125    893    +893    720  0.963      9/10

The correction is `T = (MAP - 748)/1.026`, fitted on Archive1, so applying it here
is genuinely out of sample. It cuts RMSE 42% and lifts coverage from 5/10 to 9/10 --
the offset is real and partly transferable.

But it under-corrects, and the calibration is **not stable between datasets**:

      Archive1:  MAP = 1.026*T  +748
      Archive2:  MAP = 0.928*T +2238

Archive2's raw bias is more than double Archive1's despite matching Ne, sequence
length, sample size and age range. So the offset cannot be removed by a universal
constant. The per-simulation errors are U-shaped rather than flat -- worst at the
youngest sample (+3369 on a truth of 1265, an error larger than the age itself) and
rising again at the oldest -- because a young sample has few mutations postdating
it, and the hard constraint that carries most of the information is weakest exactly
where it is needed most.

### The intervals are the right size, in the wrong place

Re-running Archive2 with Ne fixed at exactly 50,000 changes nothing -- every MAP is
identical to the Watterson-estimated run -- so the 2% rounding was irrelevant.
Reporting +/-1 block-bootstrap SD instead of percentiles separates the two failure
modes (`plot_sd.py` -> `archive2_sd.png`):

                    RMSE   bias   scatter SD   median bootstrap SD   |z|>2
      raw           1940  +1817          717                   947    4/10
      corrected     1125   +893          720                   923    1/10

The bootstrap estimates a per-simulation SD of ~947 against an actual scatter of
717, so it is slightly conservative -- it is not underestimating the noise. Every
coverage failure is therefore bias, not variance: 4 of 10 raw estimates sit beyond
2 SD of truth, falling to 1 of 10 once most of the offset is removed.

## Files

| file | purpose |
|---|---|
| `phi.py` | the `phi(k, n_T, a)` table, built in mpmath (dps=60) to survive the alternating-sum cancellation |
| `betabinom_check.py` | validates `k/(n_T+1)` against `k/n_T` for edges above `T` |
| `straddle.py` | calibration split into the three cases; exposes the straddling error |
| `calib2.py` | calibration gate for `phi` — run before any posterior |
| `post3.py` | age posterior across seeds |
| `power.py` | identifiability budget: expected log-likelihood penalty for the wrong `T` |
| `decomp2.py` | per-case obs/pred at several `T`, for locating residual bias |
| `scaling.py` | posterior versus sequence length |
| `clade.py` | tests whether mutant-clade shape carries the residual (it does not) |
| `bootstrap.py` | block bootstrap over genomic windows for honest intervals |
| `bias_vs_var.py` | many seeds at one true age: separates bias from variance |
| `converge.py` | does the estimate approach truth as sequence length grows |
| `smallp.py` | calibration on a log-`p` scale, where the coarse bins hid a 26% error |
| `numcheck.py` | separates numerical error from model error (numerics are clean) |
| `straddle_nleaf.py` | shows straddling `p` depends on clade size, which `phi` ignores |
| `phi2.py` | attempted nleaf conditioning via a second binomial -- incorrect, kept as a record |
| `phi2_check.py` | demonstrates phi2 is worse than phi |
| `run_archive.py` | runs the likelihood over a directory of supplied `.trees` simulations |
| `no_straddle.py` | tests whether the additive bias comes from the straddling term (it does not) |
| `phid.py` | `phi` with num and den exposed, so the age integral can be weighted correctly |
| `run_corrected.py` | archive run with the corrected age marginalisation; takes eps as argv[2] |
| `plot_archive.py` | estimated vs true age for both eps settings -> `age_vs_truth.png` |
| `run_archive2.py` | blind run on a panel-only ARG plus a separate haploid ancient VCF |
| `plot_archive2.py` | scores the blind predictions -> `archive2_blind.png` |
| `ages_archive2.tsv` | true ages for Archive2, supplied after the predictions were made |
| `plot_sd.py` | estimated vs true with +/-1 bootstrap SD bars -> `archive2_sd.png` |

Requires `msprime` and `mpmath` (both in `environment.yml`).
