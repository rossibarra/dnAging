# The math: dating an ancient sample from the derived alleles it carries

This note derives, from first principles, the model implemented in this folder.
Every symbol is defined in the glossary at the end; display equations are numbered
for reference.

---

## 1. Setup and goal

We have one ancestral recombination graph (ARG) inferred by SINGER (Deng, Nielsen
& Song 2025) on **a panel of 26 haplotypes** (`--n-sample`; at a site, $n$ is
however many of them are *called* there — [section 5](#5-exact-computation-via-the-neutral-moment-recursion)), summarised per site per
posterior draw as a **SNP age interval store** of mutation-age intervals, built by
the normalizeTEs pipeline (github.com/rossibarra/normalizeTEs). We have a set of
**ancient samples** (in one multi-sample VCF) genotyped at SNPs that were
**ascertained as polymorphic in a discovery panel that *contains* those ARG
haplotypes** — that panel supplies only the *site positions* and is used nowhere
else in the model, and the containment is what makes ascertainment exactly
ignorable ([NOTES.md](./NOTES.md)).

For one ancient sample we want the posterior over its age $T$ (generations before
present):

$$
p(T \mid \text{data}) \propto p(T) \mathcal{L}(T), \tag{1}
$$

with $p(T)$ a prior (uniform by default) and $\mathcal L(T)$ the likelihood over
ascertained sites $i$. Following the Poisson random field (PRF) approach, the model
uses a **composite likelihood** that treats site contributions as conditionally
independent given an ARG posterior draw. The site product is therefore taken within a draw,
$\mathcal L^{(g)}(T)=\prod_i \ell^{(g)}_i(T)$, and the ARG posterior is
marginalised as a *mixture over draws* rather than site by site ([section 6](#6-draws-polarity-chromosomes), eq. 11).

---

## 2. The per-site likelihood (with ploidy)

At site $i$ the ancient sample contributes $c_i$ **called alleles**, of which $a_i$
are the derived allele. Each called allele is a lineage drawn from the population
at time $T$, derived with probability $X_i(T)$ = the **derived-allele population
frequency** at time $T$ ([section 3](#3-what-we-compute-the-age-conditioned-expected-frequency)). Crucially $X_i(T)$ is a *random variable*: the
trajectory is unobserved, and we know only its conditional law given the site's
present count and mutation age. Adding a symmetric **per-allele** genotype-error
probability $\varepsilon$ (aDNA damage and sequencing/genotyping error in the
ancient allele calls), the per-allele probability of *observing* the derived
state, **given** the frequency, is

$$
r_i(T) = (1-\varepsilon)X_i(T) + \varepsilon\bigl(1-X_i(T)\bigr)
       = \varepsilon + (1-2\varepsilon)X_i(T). \tag{2}
$$

Here $\varepsilon$ is fixed rather than estimated, defaults to $0.01$, and must
satisfy $0\le\varepsilon<0.5$. It describes error in each ancient-VCF allele call
and nothing else. In particular it does **not** stand in for polarity uncertainty —
which allele is derived is resolved per draw from the polarity table ([section 6](#6-draws-polarity-chromosomes)) — nor for
recurrent mutation, which alters the evolutionary process rather than the
observation of it and is outside this model.

*Given* $X_i(T)$, the sample's called alleles at a site are independent draws
(Hardy–Weinberg — i.e. no recent inbreeding within the individual), so the observed
derived count is $\text{Binomial}(c_i, r_i)$. The frequency itself must then be
marginalised out, so the per-site likelihood is an **expectation** over the
conditional trajectory law (dropping the $T$-independent binomial coefficient):

$$
\ell_i(T) = \mathbb{E}\left[r_i(T)^{a_i}\bigl(1-r_i(T)\bigr)^{c_i-a_i}\right]. \tag{3}
$$

Because (2) is *affine* in $X_i(T)$, this expectation needs only the first two
conditional moments of the frequency, $\bar p_i(T)=\mathbb{E}[X_i(T)]$ and
$\bar p^{(2)}_i(T)=\mathbb{E}[X_i(T)^2]$ ([section 3](#3-what-we-compute-the-age-conditioned-expected-frequency)):

$$
\mathbb{E}[r_i] = \varepsilon + (1-2\varepsilon)\bar p_i,
\qquad
\mathbb{E}[r_i^2] = \varepsilon^2 + 2\varepsilon(1-2\varepsilon)\bar p_i +
                    (1-2\varepsilon)^2\bar p^{(2)}_i. \tag{3a}
$$

**Ploidy.** $c_i$ is set by a `--ploidy` flag matching the ancient genotype calls:

- **Haploid / pseudo-haploid** ($c_i = 1$): one allele per called site, $a_i\in\{0,1\}$.
  A homozygous call is collapsed to a single observation (crucial for pseudo-haploid
  aDNA written as `0/0` or `1/1`, so it is not counted twice). Equation (3) is
  *linear* in $r_i$, so it reduces to $\ell_i = \mathbb{E}[r_i]$ if the allele is
  derived and $1-\mathbb{E}[r_i]$ if ancestral: **the first moment alone suffices.**
- **Diploid** ($c_i = 2$): the true genotype, $a_i\in\{0,1,2\}$. Equation (3) is now
  *quadratic* in $r_i$, so the conditional second moment is required — see below.

For the diploid case the three genotype probabilities are

$$
P(a_i=2) = \mathbb{E}[r_i^2],\quad
P(a_i=1) = 2\bigl(\mathbb{E}[r_i]-\mathbb{E}[r_i^2]\bigr),\quad
P(a_i=0) = 1-2\mathbb{E}[r_i]+\mathbb{E}[r_i^2], \tag{3b}
$$

which sum to one and require the conditional **second** moment through (3a).
Substituting the plug-in mean $\bar p_i$ into a Hardy–Weinberg genotype
likelihood — i.e. using $\bar q^2,2\bar q(1-\bar q),(1-\bar q)^2$ with
$\bar q = \varepsilon+(1-2\varepsilon)\bar p_i$ — is **wrong**, because the
genotype probabilities are nonlinear in the latent frequency and in general
$\mathbb{E}[X^2]\neq\mathbb{E}[X]^2$; it understates *each* homozygote
probability by $(1-2\varepsilon)^2\mathrm{Var}(X_i(T))$ and overstates the
heterozygote probability by twice that. Using the full three-genotype form
also keeps the heterozygote-vs-homozygote information rather than collapsing to
presence/absence. A diploid site with only **one** allele called falls back to the
$c_i=1$ form, which needs the first moment only.

A missing/uncalled site has $c_i = 0$ and contributes $\log\ell_i = 0$.

**Ascertainment.** A site absent from a sample's VCF is missing data, not an
observation of "does not carry." Conditioning on the ascertained set does not bias
$\hat T$ **for this site set**, because the discovery panel *contains* the ARG panel:
we condition on $1 \le d_0 \le n-1$, i.e. both alleles seen among the $n$ ARG
haplotypes, which implies the site is polymorphic in the discovery panel — so the
ascertainment event has probability 1 given what we already condition on, and
$E[p_T \mid d_0, t_i, A] = E[p_T \mid d_0, t_i]$ holds identically. This is a
property of *these* data, not of ascertainment in general: a discovery panel
applying a frequency cutoff, or one not containing the ARG panel, would require the
ascertainment factor to be modelled. See [NOTES.md](./NOTES.md) for the condition
and what breaks it. Ascertainment does still cost power — and most where the
lower-bound signal lives (young, rare-in-discovery alleles are under-ascertained).

Everything now hinges on the conditional law of one quantity: $X_i(T)$, the
derived-allele frequency at time $T$ — and, by (3a), only on its first two moments.

---

## 3. What we compute: the age-conditioned expected frequency

We use the derived-allele frequency **conditioned on the two things known robustly
per site**: its present count in the panel, $d_0$, and its mutation age, $t_i$
(from the store). Conditioning on age encodes **survival to the present** — an
allele that is old and still segregating was, on average, at *higher* frequency in
the past than a young allele of the same present frequency (it had to be, to
persist):

$$
\bar p(T \mid d_0, t_i) = \mathbb{E}\left[X(T)\middle| d_0, t_i\right],
\qquad
\bar p^{(2)}(T \mid d_0, t_i) = \mathbb{E}\left[X(T)^2\middle| d_0, t_i\right], \tag{4}
$$

where $X(\cdot)$ is the neutral population-frequency trajectory. The second moment
is tabulated alongside the first (it costs one extra moment order, [section 5](#5-exact-computation-via-the-neutral-moment-recursion)) and is what
makes the diploid likelihood (3b) correct. In simulation the first moment
captures the survival bias. At a present count of $d_0=1$ in $n=26$ (3.85%) and a
constant $N_e=100{,}000$ — the size simulated in
[ANCIENT_TEST.md](./ANCIENT_TEST.md) — an allele of age 20,000 generations sits at
3.5% halfway back to its origin, while one of age 120,000 generations sits at
15.3%: more than four times higher at the same present frequency, purely from
having had to survive longer. This is a
classical object: Griffiths (2003, eq. 27) gives the joint density of a mutation's
population frequency and its age, conditioned on the sample count: $b$ copies out
of $n$ genes, which is exactly our conditioning with $b=d_0$. What we need is one step further: the
frequency at an intermediate time $T$ rather than at the present, which we obtain
from the same diffusion via its moments ([section 5](#5-exact-computation-via-the-neutral-moment-recursion)).

---

## 4. Time-varying $N_e$: the diffusion-time change

Under the neutral Wright–Fisher diffusion $X$ has **no drift**; its only dependence
on population size is the *rate* of drift, with infinitesimal variance
$X(1-X)/(2N_e(t))$ per generation. A time-varying $N_e(t)$ is therefore absorbed
exactly by rescaling to **diffusion time**

$$
\tau(t) = \int_0^{t} \frac{dt'}{2N_e(t')}. \tag{5}
$$

In $\tau$-time the process is the *standard*, parameter-free neutral diffusion
(Ewens 2009 lecture notes, eqs. 268–269 and 273: unit diffusion time is $2N$
generations, with drift $a(x)=0$ and variance $b(x)=x(1-x)$ — so every
constant-size result applies verbatim;
for time-varying $N_e$ this same rescaling — measuring time in units of
$\int^t dt'/2N_e$ — is the one used by Griffiths (2003, eq. 49), whose transformed
time $\int_0^t \nu(u)du$ with $\nu = N_e(0)/N_e(u)$ is $2N_e(0)\tau(t)$: the same
time change, in units of the present size rather than of the standard diffusion. The method
requires a **piecewise-constant** $N_e(t)$ — a step function of the diploid
effective size over time windows; any demographic inference expressed in that form
works. Here we use the ~50 log-spaced windows inferred by ARGtest's
`coalescence_ne_plots_from_ts.py`, where $N_e = 1/(2\text{ rate})$. Note $N_e(t)$
enters *only* through the integral (5), which for a step function is a cumulative
sum over the windows (each window contributes its width divided by $2N_e$), so each
site/draw needs just $\tau_T=\tau(T)$ and $\tau_i=\tau(t_i)$.

---

## 5. Exact computation via the neutral moment recursion

The standard neutral diffusion has (backward) generator
$\mathcal{L}f = \tfrac{x(1-x)}{2}f''(x)$ — the drift-free case
the drift-free case $\mu=0$ with $\sigma^2(x)=x(1-x)$, taken from the general
diffusion generator of Griffiths (2003, eq. 1); see also Ewens (2009 lecture notes,
eq. 218) and Kimura (1955). Applying it to $f(x)=x^k$
gives $\mathcal{L}x^k = \tfrac{k(k-1)}{2}\bigl(x^{k-1}-x^k\bigr)$, so the
moments $M_k(\tau) = \mathbb{E}[X(\tau)^k]$ obey a **closed** linear system:

$$
\frac{dM_k}{d\tau} = \frac{k(k-1)}{2}\bigl(M_{k-1} - M_k\bigr),
\qquad k = 1, 2, \dots \tag{6}
$$

$M_1$ is conserved (the frequency is a martingale). Writing $B$ for the
lower-bidiagonal generator of (6), the moment vector propagates by a matrix
exponential, $M(\tau) = e^{B\tau}M(0)$; a new mutation starts from a single copy,
$M_k(0) = \varepsilon_0^{k}$ with $\varepsilon_0 = 1/(2N_e(t_i))$.

**Sampling to the observed count.** We observe not $x_0$ but a **count**
$d_0 \sim \text{Binom}(n, x_0)$ among the $n$ called panel haplotypes. Because $n$
can vary with panel missingness, the implementation tabulates each
$n\in\{n_{\min},\ldots,26\}$ separately and drops sites below `--min-n` (default
20). The binomial pmf is a
polynomial in $x_0$ of degree $n$,

$$
P(d_0 \mid n, x_0) = \binom{n}{d_0}\sum_{m=d_0}^{n}
\binom{n-d_0}{m-d_0}(-1)^{m-d_0}x_0^{m}, \tag{7}
$$

so conditioning on $d_0$ needs frequency moments only up to order $n$ — **moments
to order $n+2$, from $(n+3)\times(n+3)$ matrix exponentials** ($n+1$ for the first
conditional moment, one further order for the second, eq. 9a) — each one tiny for
the modest $n$ of an ARG panel. The implementation evaluates **three** per table
entry ($e^{Bu_1}$, $e^{B\tau_T}$ and $e^{B\tau_i}$; eqs. 8–9), which is where the
build's cost actually sits — see the end of this section.

**Trajectory (bridge).** Let $u_1 = \tau_i - \tau_T$ be the diffusion time from
origin to sample age. With the conditional-moment map $C(\Delta)=e^{B\Delta}$
(so $\mathbb{E}[X(\Delta)^m\mid X_0=y]=\sum_j C(\Delta)_{m,j}y^{j}$), the joint
moments across the two times are

$$
\mathbb{E}\left[X_TX_{\text{pres}}^{m}\right]
= \sum_j C(\tau_T)_{m,j} M(u_1)_{j+1}, \tag{8}
$$

and the age-conditioned expected frequency is the sampling-weighted ratio

$$
\bar p(T\mid d_0,t_i) =
\frac{\displaystyle\sum_{m=d_0}^{n} \binom{n-d_0}{m-d_0}(-1)^{m-d_0}
        \mathbb{E}[X_T X_{\text{pres}}^{m}]}
     {\displaystyle\sum_{m=d_0}^{n} \binom{n-d_0}{m-d_0}(-1)^{m-d_0}
        M_{\text{pres},m}},
\qquad M_{\text{pres}} = e^{B\tau_i}M(0). \tag{9}
$$

**The conditional second moment** required by the diploid likelihood (3a–3b) comes
from the *same* contraction shifted one index. The $(m,j)$ entry of $C(\tau_T)$ is
the coefficient of $x^j$ in the conditional moment
$\mathbb{E}[X^m_{\mathrm{pres}}\mid X_T=x]$, so each extra factor of $X_T$ raises that
power by one:

$$
\mathbb{E}\left[X_T^2X_{\text{pres}}^{m}\right]
= \sum_j C(\tau_T)_{m,j} M(u_1)_{j+2},
\qquad
\bar p^{(2)}(T\mid d_0,t_i) =
\frac{\sum_{m} \binom{n-d_0}{m-d_0}(-1)^{m-d_0}
        \mathbb{E}[X_T^2 X_{\text{pres}}^{m}]}
     {\sum_{m} \binom{n-d_0}{m-d_0}(-1)^{m-d_0} M_{\text{pres},m}}, \tag{9a}
$$

with the **same denominator** as (9) — the sampling weight does not change, only the
functional being averaged. The cost is one extra moment order in $M(u_1)$. The
result necessarily obeys $\bar p^{2} \le \bar p^{(2)} \le \bar p$ (Cauchy–Schwarz,
and $X^2\le X$ on $[0,1]$).

For $T \ge t_i$ (sample older than the mutation) the allele does not yet exist, so
$\bar p = \bar p^{(2)} = 0$. This is exact and was validated against a forward
Wright–Fisher Monte Carlo — agreement to MC noise for **both** moments, including
the rare $d_0=2$ bin (0.22 vs 0.23 for $\bar p$; 0.066 vs 0.067 for
$\bar p^{(2)}$), the $T\ge t_i$ boundary, and Kimura's constant-size limit
(`validate_moments_vs_mc.py`). The tables $\bar p(T\mid d_0,t_i)$ and
$\bar p^{(2)}(T\mid d_0,t_i)$ are built over $(n, d_0, t_i, T)$ — a separate plane
for each called-panel size $n$, as above — once, and shared across all sites and
all samples.

**Exact in the model, numerical in the table.** Equations (6)–(10b) are exact
within the neutral diffusion, but what inference reads is a *tabulation* of them:
entries are stored in `float32`, the mutation-age axis is a finite log-spaced grid
recovered by interpolation in $\log t_i$, the branch integral (10b) is a 16-node
trapezoidal quadrature, mutation-age mass beyond $\tau_i=3$ is truncated (below),
and the conditioning sums lose precision as described next. Where this note calls
something exact, it means exact *as an identity of the model*, never that the
pipeline's arithmetic is.

**Numerics / scaling with $n$.** The alternating-sign conditioning sums in (7), (9)
and (9a) are prone to catastrophic cancellation: the binomial weights grow like
$2^{n}$ while the moments are small, so the sum loses on the order of $0.3n$
decimal digits. The loss also grows with $\tau_i$: once the allele is almost surely
lost or fixed, the denominator of (9)/(9a) — the sampling weight of an intermediate
$d_0$ — underflows to cancellation noise, and *both* moments become meaningless
(clipping them into the mathematically valid range would hide the failure rather
than recover the information). The implementation therefore measures
$A=\sum|\mathrm{term}|/|\sum\mathrm{term}|$ for each conditioning sum and writes
`NaN` once only roughly 1–2 significant decimal digits remain. Inference treats a
`NaN` in **any** ARG draw as disqualifying the whole site: the mixture in eq. (11)
is defined over all $M$ draws, so a single draw cannot simply be dropped ([section 6](#6-draws-polarity-chromosomes)).

At $n=26$ the relative error against an 80-digit `mpmath` reference
(`tests/_reference.py`, `dps=80`) is $\sim10^{-3}$ for
$\tau_i\lesssim3$, reaches a few percent by $\tau_i\approx 5$ to $6$, and the values are
noise by $\tau_i\approx10$. As an additional operational safeguard,
`--mutation-age-max` defaults to a cutoff of $\tau_i=3$ in diffusion units and
discards mutation-age mass beyond it. The table stores the diffusion time
$\tau(t)=\int_0^t ds/(2N_e(s))$ of every generation-age row, and inference recovers
the cutoff *in generations* by **linear interpolation between those log-spaced
rows**; it does not re-evaluate $N_e(t)$. For a piecewise-constant $N_e$, $\tau(t)$
is piecewise linear in $t$, so the inverse would be exact whenever the bracketing
pair of rows lies inside a single demographic window — but a window boundary
falling between two neighbouring rows leaves the recovered cutoff age approximate,
by up to the spacing of the age grid (100 log-spaced rows by default, i.e. ~17% in
$t$ per step). At constant $N_e=100{,}000$, $\tau=3$ is about 600,000 generations.
Table construction requires `--age-max` to extend beyond $\tau=3$, and inference
rejects insufficient table coverage. This is a numerical cutoff, not an assertion
that all older mutations are biologically uninformative. Double precision is
therefore comfortable for panels up to a few tens of haplotypes (the ARG regime
here) at moderate $\tau_i$, but the method degrades and eventually breaks beyond
roughly $n \approx 40$ to $50$ — there one must switch to extended precision
(e.g. `mpmath`) or reformulate the conditioning in a numerically stable basis
(orthogonal-polynomial / spectral moments rather than the raw power moments).

**Cost of the build.** A *single* matrix exponential is cheap ($O((n+3)^3)$ on a
matrix this small), but the number of them is not negligible: `Emoments()`
recomputes all three per $(n,d_0,t_i,T)$ entry, and the build loops over every
panel size $n$, every $d_0\in\{1,\ldots,n\}$, every mutation age and every sample
age (entries with $T\ge t_i$ return $0$ before doing any work). The default grid
($n=20\ldots26$, $d_0$ up to $n$, 100 log-spaced ages, 300 sample ages) issues
$\approx7.8\times10^{6}$ exponentials of dimension 23–29, of order an hour on one
core — that is the dominant term in precomputation, not a rounding error. All three
depend only on $(n,\tau_i,\tau_T)$ and **not** on $d_0$, so hoisting them out of the
$d_0$ loop would cut the count by a factor of $n$, to $\approx3.4\times10^{5}$.
That reuse is **not** implemented; it is the obvious optimisation should the grid
ever be refined.

---

## 6. Draws, polarity, chromosomes

**Polarity / ALT convention.** Which allele is "derived" can flip between ARG
draws, so we track the **ALT-allele frequency** $\varphi_{\text{alt}}$
consistently. With $c_{\text{alt}}$ = ALT count in the panel and, for draw $g$, the
ancestral base from the polarity table,

$$
\varphi_{\text{alt}}^{(g)}(T) =
\begin{cases}
\bar p(T \mid d_0=c_{\text{alt}}, t_i^{(g)}), & \text{ALT derived in draw }g,\cr
1 - \bar p(T \mid d_0=n-c_{\text{alt}}, t_i^{(g)}), & \text{ALT ancestral in draw }g,
\end{cases}
\tag{10}
$$

and, for the diploid likelihood, the matching **second** moment of the same ALT
frequency — where the ALT-ancestral branch must transform *both* moments, since
$\mathbb{E}[(1-X)^2] = 1-2\mathbb{E}[X]+\mathbb{E}[X^2]$:

$$
\varphi^{(2,g)}_{\text{alt}}(T) =
\begin{cases}
\bar p^{(2)}(T \mid c_{\text{alt}}, t_i^{(g)}), & \text{ALT derived},\cr
1 - 2\bar p(T \mid n-c_{\text{alt}}, t_i^{(g)}) +
  \bar p^{(2)}(T \mid n-c_{\text{alt}}, t_i^{(g)}), & \text{ALT ancestral}.
\end{cases}
\tag{10a}
$$

Because $c_{\text{alt}}$ comes from the panel VCF while $a_i, c_i$ come from the
ancient VCF, the two files' REF/ALT representations must agree: the implementation
harmonises them per site (using $n-c_{\text{alt}}$ where they are exactly swapped)
and **skips** any site whose alleles cannot be matched, rather than joining on
position alone.

Here $t_i^{(g)}$ is the mutation age in draw $g$. **This age is uncertain, and the
uncertainty is marginalised at two levels.** *Within* a draw the age is uniform on
its branch $[\text{below}_g,\text{above}_g]$ (from the store; infinite-sites model),
so the value entering eq. (10) is the branch integral

$$
\bar p\big(T \mid d_0,[\text{below}_g,\text{above}_g]\big)
= \frac{1}{\text{above}_g-\text{below}_g}
  \int_{\text{below}_g}^{\text{above}_g} \bar p(T\mid d_0, t)dt .
\tag{10b}
$$

**The lower limit may fall outside the table, and this costs nothing.**
`below` can lie below the youngest tabulated mutation age $t_0$ — it
is $0$ for a branch reaching the present. The integral is then evaluated only over
the covered part $[t_{\min},\text{above}_g]$, but normalised by the **true** branch
length the true branch length. That is *exact* for $T \ge t_0$:
every uncovered age satisfies $t < t_0 \le T$, so the mutation postdates the
sample, $\bar p(T\mid d_0,t)=0$, and the uncovered region contributes nothing to
the numerator of (10b). Normalising by the covered width instead would inflate the
site by

$$
\frac{\text{above}_g-\text{below}_g}{\text{above}_g-t_{\min}}.
$$

The correction is **asymmetric by construction**: it is not applied at the upper
limit, where uncovered ages are *older* than the sample and genuinely contribute.
That end is bounded instead by the `--mutation-age-max` cutoff, which caps
`above` and rejects tables that do not reach it. For $T < t_0$ the
identity fails — uncovered ages in $(T,t_0)$ do contribute — so the posterior
below $t_0$ is biased low; see NOTES.md.

Exactly one mapped branch interval is required per site and ARG draw. If any draw
supplies multiple intervals, the mutation is treated as multiply mapped and the
entire site is excluded rather than assigning weights among mappings that are not
trusted. The requirement is symmetric: **every** one of the $M$ draws must supply
an interval and pass the polarity, age-cutoff, table-coverage and numerical checks,
or the site is dropped. Because eq. (11) below is an equal-weight mixture over all
$M$ draws, keeping a site that survives in only a subset would silently replace the
estimand with an expectation conditional on that subset, and the missing terms
cannot be imputed (a missing polarity call is not evidence of zero ALT frequency).
`run.json` counts such sites under `sites_incomplete_draws`, with per-cause draw
tallies alongside.

*Between* draws, averaging over the $M$ posterior draws integrates the remaining
posterior on $t_i$ — each draw places the mutation on a different branch — so
together with (10b) this marginalises the full ARG posterior over the mutation age,
within-branch and between-draw. But that average is taken **at the level of the
likelihood, not of the frequency.** One ARG draw is a single chromosome-wide
genealogy: it fixes the mutation age at *every* site at once, so the sites' ages are
dependent, and that dependence is carried by the draw index. Conditional on a draw,
the PRF composite likelihood *models* the remaining per-site observations as
independent. So (10) and (10a) enter (3a) **within** a draw —
$\bar p_i\to\varphi^{(g)}(T)$ and $\bar p^{(2)}\to\varphi^{(2,g)}(T)$ — writing
$\varphi^{(g)}$ for $\varphi^{(g)}_{\text{alt}}$ from here on — with $a_i,c_i$ the
sample's ALT dosage and ploidy at the site (ALT carriage is polarity-independent) —
giving a per-draw per-site likelihood $\ell^{(g)}_i(T)$, and the draws are averaged
only *after* the site product:

$$
\mathcal{L}_c(T) = \frac{1}{M}\sum_{g=1}^{M}\ \prod_i\ \ell^{(g)}_i(T),
\qquad
\ell^{(g)}_i(T) = \mathbb{E}\left[r^{(g)}_i(T)^{a_i}
      \bigl(1-r^{(g)}_i(T)\bigr)^{c_i-a_i}\right], \tag{11}
$$

where the expectation is the one in (3), evaluated through (3a) with draw $g$'s
moments $\varphi^{(g)}$ and $\varphi^{(2,g)}$. The
implementation accumulates this in logs, keeping a chromosome's $M$ per-draw
log-likelihoods separate until every site has been multiplied in and only then
$\log\mathcal L_c = \mathrm{logsumexp}_g \sum_i \log\ell^{(g)}_i - \log M$.

**Across-site composite likelihood.** Equation (11) is a PRF-style composite
likelihood, not a claim that linked SNPs are literally independent. This is the
standard scalable approximation behind frequency-spectrum inference (Sawyer &
Hartl 1992; Gutenkunst et al. 2009): linkage changes the joint variance while the
one-site marginal model remains the same. Local LD can therefore make likelihood
curvature overstate the amount of independent information, so intervals derived
from the normalized composite likelihood are nominal rather than guaranteed to
have exact frequentist coverage. For the intended maize application, however,
millions of SNPs are distributed across more than 10 Morgans and hence many
recombining genealogical blocks. We retain all quality-controlled sites rather
than impose arbitrary physical thinning; block bootstrap or simulation calibration
is the appropriate route if coverage itself is a target.

Averaging the *frequency* per site instead — using
$\bar\varphi_{\text{alt}} = \frac1M\sum_g \varphi^{(g)}_{\text{alt}}$ in (3a) and
then taking $\prod_i$ — is a **different and wrong** quantity: it discards the
across-site coupling an ARG draw carries, treating each site's age as if it were
redrawn independently at every site. The site product is nonlinear in the per-draw
likelihood, so the two disagree as soon as a chromosome has more than one site (and
under `--ploidy 2` already at a single site, since (3b) is nonlinear in the
frequency; under `--ploidy 1` a lone site is the one case where they coincide). Two
sites whose per-draw frequencies share the same draw *means* but differ in their
*pairing across draws* have different likelihoods under (11) and identical ones
under the per-site average (`tests/test_draw_marginalization.py`). Each draw's
$\varphi^{(g)}$ and $\varphi^{(2,g)}$ come from the same
conditional law, so $(\varphi^{(g)})^2\le\varphi^{(2,g)}\le\varphi^{(g)}$ holds
draw by draw, and the `--ploidy 1` path never touches the second moment.

**Chromosomes.** ARG draws are sampled independently per chromosome, and the $C$
chromosomes are unlinked — hence conditionally independent given $T$ — so the joint
likelihood factorises and the per-sample log-likelihood is the sum of the
per-chromosome terms:

$$
\log \mathcal{L}(T) = \sum_{c=1}^{C} \log \mathcal{L}_{c}(T), \tag{12}
$$

where $\mathcal L_c$ is the within-chromosome likelihood of eq. (11), over that
chromosome's sites and its own ARG draws. Each chromosome carries its own
independent draw mixture, which is why the marginalisation in (11) is per
chromosome and only the resulting log-marginals are summed here. Each run computes
one $\log\mathcal L_c$; the merge step sums them over the $C$ chromosomes. (Here
$c$ indexes chromosomes, distinct from the per-site allele count $c_i$.) Sites
monomorphic in the panel ($c_{\text{alt}}\in\{0,n\}$) carry no trajectory and are
skipped.

---

## 7. Assumptions and caveats

1. **Neutrality.** Eq. (4) is the *neutral* trajectory — right for the vast
   majority of sites; strongly selected sites are not neutral. Guard with an
   approximately-neutral `--include-positions` set if selection is a concern.
2. **Ascertainment caps power** ([section 2](#2-the-per-site-likelihood-with-ploidy)): expect broad posteriors on array-ascertained
   SNPs.
3. **The ancient sample must not be in the ARG/ascertainment panels.** Otherwise
   its own alleles are among the $n$ that produced $d_0$, and the independence
   eq. (3) assumes between the panel count and the ancient genotype breaks — the
   sample would partly be conditioned on itself.
4. **Hardy–Weinberg** within a diploid individual assumes no recent inbreeding —
   the two alleles are independent *given* the population frequency, which is why
   (3b) needs $\mathbb{E}[r^2]$ rather than $\mathbb{E}[r]^2$; for
   haploid/pseudo-haploid data this is moot ($c_i=1$).
5. **Effective-size convention.** The curve is a diploid effective size
   with $N_e = 1/(2\text{ rate})$,
   so a single new copy has frequency $1/(2N_e)$ and drift variance $X(1-X)/(2N_e)$.
6. **Only two moments.** Eq. (3a) is exact for $c_i\le 2$. A likelihood over three
   or more alleles drawn from the same latent frequency would need moments of $X$ up
   to that order; the table carries only the first two.
7. **Panel missingness is allele-blind**, hence ignorable given the called count
   ([section 5](#5-exact-computation-via-the-neutral-moment-recursion)). Eq. (7) treats the $n$ called panel haplotypes as a binomial sample of the
   population, i.e. as an allele-blind subset of the 26: choosing the moment plane
   for the site's exact $n$ makes the sampling distribution right for that $n$, but
   it conditions on *how many* haplotypes were called, not on *which*.
   Allele-dependent missingness (reference bias making one allele harder to call)
   therefore biases $d_0$ and is
   **not** removed by the $n$-specific plane. See [NOTES.md](./NOTES.md) for the
   condition and the mitigation.

---

## Glossary of variables

| symbol | meaning |
|---|---|
| $T$ | age of the ancient sample, generations before present (inferred) |
| $t_i$ | age of the mutation at site $i$ (interval $[\text{below},\text{above}]$ per draw, from the store) |
| $n$ | number of ARG-panel haplotypes **called at the site** (`--min-n` $\le n\le 26$; the full panel size only where the panel is fully called — [section 5](#5-exact-computation-via-the-neutral-moment-recursion)) |
| $C$ | number of chromosomes (independent given $T$) |
| $d_0$ | present count of the derived allele among the $n$ *called* panel haplotypes |
| $c_{\text{alt}}$ | present count of the ALT allele among the $n$ *called* panel haplotypes (from the panel VCF) |
| $x_0$ | present *population* derived frequency, unobserved, with $d_0\sim\text{Binom}(n,x_0)$ |
| $a_i$ | number of derived (ALT) alleles the ancient sample shows at site $i$ |
| $c_i$ | number of called alleles at site $i$ (ploidy: 1 haploid/pseudo-haploid, 2 diploid, 0 missing) |
| $X_i(T),X(T)$ | derived-allele population frequency at time $T$ (a random variable) |
| $\bar p(T\mid d_0,t_i)$ | expected frequency at $T$ given present count and age — the table value (eq. 4, 9) |
| $\bar p^{(2)}(T\mid d_0,t_i)$ | expected *squared* frequency — the second table plane (eq. 4, 9a) |
| $r_i(T)$ | per-allele probability of *observing* derived given the frequency, with error $\varepsilon$ (eq. 2) |
| $\varphi_{\text{alt}}^{(g)}$ | ALT-allele frequency at $T$ in draw $g$ (eq. 10) |
| $\varphi^{(2,g)}_{\text{alt}}$ | the corresponding second moment (eq. 10a) |
| $N_e(t)$ | diploid effective population size at time $t$, equal to $1/(2\text{ rate})$ |
| $\tau(t)$ | diffusion time, $\int_0^t dt'/(2N_e(t'))$ |
| $\tau_T,\tau_i$ | $\tau(T)$ and $\tau(t_i)$ |
| $\varepsilon_0$ | frequency of a single new copy at origin, $1/(2N_e(t_i))$ |
| $M_k(\tau)$ | $k$-th moment $\mathbb{E}[X(\tau)^k]$ of the neutral diffusion |
| $B$ | generator of the closed moment recursion (eq. 6) |
| $C(\Delta)=e^{B\Delta}$ | conditional-moment map over diffusion-time $\Delta$ |
| $\varepsilon$ | symmetric per-allele genotype-error probability |
| $\ell_i(T),\mathcal L(T)$ | per-site and total likelihood |
| $\ell^{(g)}_i(T)$ | per-site likelihood in draw $g$ — (3) with draw $g$'s moments (eq. 11) |
| $\mathcal L_c(T)$ | within-chromosome likelihood: draw mixture of site products (eq. 11, 12) |
| $g,M$ | ARG posterior draw index and number of draws |

---

## References

- Kimura, M. (1955). *Solution of a process of random genetic drift with a
  continuous model.* PNAS 41:144–150.
- Sawyer, S.A. & Hartl, D.L. (1992). *Population genetics of polymorphism and
  divergence.* Genetics 132:1161–1176.
  [doi:10.1093/genetics/132.4.1161](https://doi.org/10.1093/genetics/132.4.1161)
- Gutenkunst, R.N., Hernandez, R.D., Williamson, S.H. & Bustamante, C.D. (2009).
  *Inferring the joint demographic history of multiple populations from
  multidimensional SNP frequency data.* PLoS Genetics 5:e1000695.
  [doi:10.1371/journal.pgen.1000695](https://doi.org/10.1371/journal.pgen.1000695)
- Griffiths, R.C. (2003). *The frequency spectrum of a mutation, and its age, in a
  general diffusion model.* Theoretical Population Biology 64:241–251.
- Ewens, W.J. (2009). *Mathematical Population Genetics: Introduction to the
  Stochastic Theory* (lecture notes, Guanajuato, March 2009; abstracted from Ewens
  2004, Springer).
  [PDF](http://www.rgwinther.com/Ewens2009MathematicalPopulationGeneticsTheGuanajuatoLectures.pdf)
- Deng, Y., Nielsen, R. & Song, Y.S. (2025). *Robust and accurate Bayesian
  inference of genome-wide genealogies for hundreds of genomes.* Nature Genetics
  57:2124–2135 (doi:10.1038/s41588-025-02317-9); software:
  github.com/popgenmethods/SINGER.
- bis101 & Ross-Ibarra, J. (2026). *RILAB/argtest: ARGtest v1.10.* Zenodo
  (doi:10.5281/zenodo.21865676).
- normalizeTEs pipeline — github.com/rossibarra/normalizeTEs.
