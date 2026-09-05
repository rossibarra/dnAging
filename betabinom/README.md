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

The bias source is still unidentified. Calibration is flat in every stratum tested
(by `p`, by fraction of edge above `T`, by clade compression, and in aggregate at
the true age it is the *best* of any `T`), yet the likelihood peaks several
thousand generations too old. A model whose per-site mean is right everywhere
tested but whose likelihood peaks in the wrong place is misspecified in something
not yet measured.

Reported intervals are also not trustworthy at any stage: they treat linked sites
as independent, so even an exactly specified model would give intervals far too
narrow. Coverage needs a block bootstrap.

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

Requires `msprime` and `mpmath` (both in `environment.yml`).
