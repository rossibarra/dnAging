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

Three findings about why:

1. **The data is informative; the model is not accurate enough.** The expected
   log-likelihood penalty for the wrong `T`, if the model were exact, is ~100 nats
   at `T = 7658` (`power.py`). But misspecification supplies ~129 nats of spurious
   gain in the same direction, so it outbids the signal.
2. **More sites will not fix it.** Signal and systematic bias both scale with site
   count, so the peak location is invariant while the width shrinks as
   `1/sqrt(N)`. Going 0.5 -> 4 Mb leaves the MAP wrong and only narrows the
   interval (`scaling.py`).
3. **`(k, n_T, a)` is not a sufficient statistic for the tree.** At the true age
   the aggregate calibration is now *best* (obs/pred 1.0019, versus 0.9954 and
   1.0091 at neighbouring `T`), yet the likelihood peaks elsewhere. A model whose
   mean is right but whose likelihood peaks in the wrong place is misspecified in
   the spread of `p`, not its level. Two sites with the same `(k, n_T, a)` get the
   same `p`, but a mutant clade that coalesces quickly implies a lower frequency
   than one that does not.

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

Requires `msprime` and `mpmath` (both in `environment.yml`).
