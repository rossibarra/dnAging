# The math: dating an ancient sample from the derived alleles it carries

This note derives, from first principles, the model implemented in this folder.
Every symbol is defined in the glossary at the end; display equations are numbered
for reference.

---

## 1. Setup and goal

We have one ancestral recombination graph (ARG) inferred by SINGER (Deng, Nielsen
& Song 2025) on **a panel of $n$ haplotypes**, summarised per site per
posterior draw as a **SNP age interval store** of mutation-age intervals, built by
the normalizeTEs pipeline (github.com/rossibarra/normalizeTEs). We have a set of
**ancient samples** (in one multi-sample VCF) genotyped at SNPs that were
**ascertained as polymorphic in a discovery panel that *contains* the $n$ ARG
haplotypes** — that panel supplies only the *site positions* and is used nowhere
else in the model, and the containment is what makes ascertainment exactly
ignorable ([NOTES.md](./NOTES.md)).

For one ancient sample we want the posterior over its age $T$ (generations before
present):

$$
p(T \mid \text{data}) \;\propto\; p(T)\, \mathcal{L}(T), \tag{1}
$$

with $p(T)$ a prior (uniform by default) and $\mathcal L(T)=\prod_i \ell_i(T)$ the
likelihood over ascertained sites $i$.

---

## 2. The per-site likelihood (with ploidy)

At site $i$ the ancient sample contributes $c_i$ **called alleles**, of which $a_i$
are the derived allele. Each called allele is a lineage drawn from the population
at time $T$, derived with probability $X_i(T)$ = the **derived-allele population
frequency** at time $T$ (§3). Crucially $X_i(T)$ is a *random variable*: the
trajectory is unobserved, and we know only its conditional law given the site's
present count and mutation age. Adding a symmetric **per-allele** genotype-error
probability $\varepsilon$ (aDNA damage, sequencing/genotyping error, recurrent
mutation, mis-polarisation), the per-allele probability of *observing* the derived
state, **given** the frequency, is

$$
r_i(T) = (1-\varepsilon)\,X_i(T) + \varepsilon\,\bigl(1-X_i(T)\bigr)
       = \varepsilon + (1-2\varepsilon)\,X_i(T). \tag{2}
$$

*Given* $X_i(T)$, the sample's called alleles at a site are independent draws
(Hardy–Weinberg — i.e. no recent inbreeding within the individual), so the observed
derived count is $\text{Binomial}(c_i, r_i)$. The frequency itself must then be
marginalised out, so the per-site likelihood is an **expectation** over the
conditional trajectory law (dropping the $T$-independent binomial coefficient):

$$
\ell_i(T) = \mathbb{E}\!\left[\,r_i(T)^{\,a_i}\bigl(1-r_i(T)\bigr)^{\,c_i-a_i}\,\right]. \tag{3}
$$

Because (2) is *affine* in $X_i(T)$, this expectation needs only the first two
conditional moments of the frequency, $\bar p_i(T)=\mathbb{E}[X_i(T)]$ and
$\bar p^{(2)}_i(T)=\mathbb{E}[X_i(T)^2]$ (§3):

$$
\mathbb{E}[r_i] = \varepsilon + (1-2\varepsilon)\,\bar p_i,
\qquad
\mathbb{E}[r_i^2] = \varepsilon^2 + 2\varepsilon(1-2\varepsilon)\,\bar p_i
                    + (1-2\varepsilon)^2\,\bar p^{(2)}_i. \tag{3a}
$$

**Ploidy.** $c_i$ is set by a `--ploidy` flag matching the ancient genotype calls:

- **Haploid / pseudo-haploid** ($c_i = 1$): one allele per called site, $a_i\in\{0,1\}$.
  A homozygous call is collapsed to a single observation (crucial for pseudo-haploid
  aDNA written as `0/0` or `1/1`, so it is not counted twice). Equation (3) is
  *linear* in $r_i$, so it reduces to $\ell_i = \mathbb{E}[r_i]$ if the allele is
  derived and $1-\mathbb{E}[r_i]$ if ancestral: **the first moment alone suffices.**
- **Diploid** ($c_i = 2$): the true genotype, $a_i\in\{0,1,2\}$. Equation (3) is now
  *quadratic* in $r_i$, so the three genotype probabilities are

$$
P(a_i=2) = \mathbb{E}[r_i^2],\quad
P(a_i=1) = 2\bigl(\mathbb{E}[r_i]-\mathbb{E}[r_i^2]\bigr),\quad
P(a_i=0) = 1-2\,\mathbb{E}[r_i]+\mathbb{E}[r_i^2], \tag{3b}
$$

  which sum to one and require the conditional **second** moment through (3a).
  Substituting the plug-in mean $\bar p_i$ into a Hardy–Weinberg genotype
  likelihood — i.e. using $\bar q^2,\,2\bar q(1-\bar q),\,(1-\bar q)^2$ with
  $\bar q = \varepsilon+(1-2\varepsilon)\bar p_i$ — is **wrong**, because the
  genotype probabilities are nonlinear in the latent frequency and in general
  $\mathbb{E}[X^2]\neq\mathbb{E}[X]^2$; it understates homozygote probabilities and
  overstates heterozygote probabilities by
  $(1-2\varepsilon)^2\operatorname{Var}(X_i(T))$. Using the full three-genotype form
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
\bar p(T \mid d_0, t_i) = \mathbb{E}\!\left[\,X(T)\,\middle|\, d_0,\, t_i\,\right],
\qquad
\bar p^{(2)}(T \mid d_0, t_i) = \mathbb{E}\!\left[\,X(T)^2\,\middle|\, d_0,\, t_i\,\right], \tag{4}
$$

where $X(\cdot)$ is the neutral population-frequency trajectory. The second moment
is tabulated alongside the first (it costs one extra moment order, §5) and is what
makes the diploid likelihood (3b) correct. In simulation the first moment
captures the survival bias: at the same present frequency (~5%), an allele of age
200 gen sits near 3.6% halfway back, one of age 1200 gen near 14.8%. This is a
classical object: the joint density of a mutation's frequency and age is given by
Griffiths (2003, eq. 27). We need the frequency conditional on age, which we obtain
from the same diffusion via its moments (§5).

---

## 4. Time-varying $N_e$: the diffusion-time change

Under the neutral Wright–Fisher diffusion $X$ has **no drift**; its only dependence
on population size is the *rate* of drift, with infinitesimal variance
$X(1-X)/(2N_e(t))$ per generation. A time-varying $N_e(t)$ is therefore absorbed
exactly by rescaling to **diffusion time**

$$
\tau(t) = \int_0^{t} \frac{dt'}{2\,N_e(t')}. \tag{5}
$$

In $\tau$-time the process is the *standard*, parameter-free neutral diffusion
(Ewens 2009 lecture notes, eq. 211), so all constant-$N_e$ results apply verbatim;
for time-varying $N_e$ this same rescaling — measuring time in units of
$\int^t dt'/2N_e$ — is the one used by Griffiths (2003, eq. 49). The method
requires a **piecewise-constant** $N_e(t)$ — a step function of the diploid
effective size over time windows; any demographic inference expressed in that form
works. Here we use the ~50 log-spaced windows inferred by ARGtest's
`coalescence_ne_plots_from_ts.py` (where $N_e = 1/(2\,\text{rate})$). $N_e(t)$
enters *only* through the integral (5), which for a step function is a cumulative
sum over the windows (each window contributes its width divided by $2N_e$), so each
site/draw needs just $\tau_T=\tau(T)$ and $\tau_i=\tau(t_i)$.

---

## 5. Exact computation via the neutral moment recursion

The standard neutral diffusion has (backward) generator
$\mathcal{L}f = \tfrac{x(1-x)}{2}\,f''(x)$ — the drift-free case
($\mu=0$, $\sigma^2(x)=x(1-x)$) of the general diffusion generator (Griffiths 2003,
eq. 1; Ewens 2009 lecture notes, eq. 218; Kimura 1955). Applying it to $f(x)=x^k$
gives $\mathcal{L}x^k = \tfrac{k(k-1)}{2}\bigl(x^{k-1}-x^k\bigr)$, so the
moments $M_k(\tau) = \mathbb{E}[X(\tau)^k]$ obey a **closed** linear system:

$$
\frac{dM_k}{d\tau} = \frac{k(k-1)}{2}\bigl(M_{k-1} - M_k\bigr),
\qquad k = 1, 2, \dots \tag{6}
$$

$M_1$ is conserved (the frequency is a martingale). Writing $B$ for the
lower-bidiagonal generator of (6), the moment vector propagates by a matrix
exponential, $M(\tau) = e^{B\tau}\,M(0)$; a new mutation starts from a single copy,
$M_k(0) = \varepsilon_0^{\,k}$ with $\varepsilon_0 = 1/(2N_e(t_i))$.

**Sampling to the observed count.** We observe not $x_0$ but a **count**
$d_0 \sim \text{Binom}(n, x_0)$ in the $n$-haplotype panel. The binomial pmf is a
polynomial in $x_0$ of degree $n$,

$$
P(d_0 \mid n, x_0) = \binom{n}{d_0}\sum_{m=d_0}^{n}
\binom{n-d_0}{m-d_0}(-1)^{m-d_0}\,x_0^{\,m}, \tag{7}
$$

so conditioning on $d_0$ needs frequency moments only up to order $n$ — **moments
to order $n+2$, one $(n+3)\times(n+3)$ matrix exponential** ($n+1$ for the first
conditional moment, one further order for the second, eq. 9a) — tiny for the modest
$n$ of an ARG panel.

**Trajectory (bridge).** Let $u_1 = \tau_i - \tau_T$ be the diffusion time from
origin to sample age. With the conditional-moment map $C(\Delta)=e^{B\Delta}$
(so $\mathbb{E}[X(\Delta)^m\mid X_0=y]=\sum_j C(\Delta)_{m,j}\,y^{\,j}$), the joint
moments across the two times are

$$
\mathbb{E}\!\left[X_T\,X_{\text{pres}}^{\,m}\right]
= \sum_j C(\tau_T)_{m,j}\; M(u_1)_{\,j+1}, \tag{8}
$$

and the age-conditioned expected frequency is the sampling-weighted ratio

$$
\bar p(T\mid d_0,t_i) =
\frac{\displaystyle\sum_{m=d_0}^{n} \binom{n-d_0}{m-d_0}(-1)^{m-d_0}\,
        \mathbb{E}[X_T X_{\text{pres}}^{\,m}]}
     {\displaystyle\sum_{m=d_0}^{n} \binom{n-d_0}{m-d_0}(-1)^{m-d_0}\,
        M_{\text{pres},\,m}},
\qquad M_{\text{pres}} = e^{B\tau_i}M(0). \tag{9}
$$

**The conditional second moment** required by the diploid likelihood (3a–3b) comes
from the *same* contraction shifted one index. Since $C(\tau_T)_{m,j}$ is the
coefficient of $x^{\,j}$ in $\mathbb{E}[X_{\text{pres}}^{\,m}\mid X_T=x]$, each extra
factor of $X_T$ raises that power by one, so

$$
\mathbb{E}\!\left[X_T^2\,X_{\text{pres}}^{\,m}\right]
= \sum_j C(\tau_T)_{m,j}\; M(u_1)_{\,j+2},
\qquad
\bar p^{(2)}(T\mid d_0,t_i) =
\frac{\sum_{m} \binom{n-d_0}{m-d_0}(-1)^{m-d_0}\,
        \mathbb{E}[X_T^2 X_{\text{pres}}^{\,m}]}
     {\sum_{m} \binom{n-d_0}{m-d_0}(-1)^{m-d_0}\, M_{\text{pres},\,m}}, \tag{9a}
$$

with the **same denominator** as (9) — the sampling weight does not change, only the
functional being averaged. The cost is one extra moment order in $M(u_1)$. The
result necessarily obeys $\bar p^{\,2} \le \bar p^{(2)} \le \bar p$ (Cauchy–Schwarz,
and $X^2\le X$ on $[0,1]$).

For $T \ge t_i$ (sample older than the mutation) the allele does not yet exist, so
$\bar p = \bar p^{(2)} = 0$. This is exact and was validated against a forward
Wright–Fisher Monte Carlo — agreement to MC noise for **both** moments, including
the rare $d_0=2$ bin (0.22 vs 0.23 for $\bar p$; 0.066 vs 0.067 for
$\bar p^{(2)}$), the $T\ge t_i$ boundary, and Kimura's constant-$N_e$ limit
(`validate_moments_vs_mc.py`). The tables $\bar p(T\mid d_0,t_i)$ and
$\bar p^{(2)}(T\mid d_0,t_i)$ over $(d_0, t_i, T)$ are built once and shared across
all sites and all samples.

**Numerics / scaling with $n$.** The alternating-sign conditioning sums in (7), (9)
and (9a) are prone to catastrophic cancellation: the binomial weights grow like
$2^{n}$ while the moments are small, so the sum loses on the order of $0.3\,n$
decimal digits. The loss also grows with $\tau_i$: once the allele is almost surely
lost or fixed, the denominator of (9)/(9a) — the sampling weight of an intermediate
$d_0$ — underflows to cancellation noise, and *both* moments become meaningless
(clipping them to $[0,1]$ and to $[\bar p^{\,2},\bar p]$ keeps them usable but does
not recover the information). At $n=26$ the relative error against a 60-digit
reference is $\sim10^{-3}$ for $\tau_i\lesssim3$, reaches a few percent by
$\tau_i\approx5$–$6$, and the values are noise by $\tau_i\approx10$ — so the
`--age-max` end of the mutation-age grid should be read with that in mind. Double
precision is therefore comfortable for panels up to a few tens of haplotypes (the
ARG regime here) at moderate $\tau_i$, but the method degrades and eventually breaks
beyond roughly $n \approx 40$–$50$; there one must switch to extended precision
(e.g. `mpmath`) or reformulate the conditioning in a numerically stable basis
(orthogonal-polynomial / spectral moments rather than the raw power moments). The
matrix-exponential cost itself is negligible ($(n+3)^3$).

---

## 6. Draws, polarity, chromosomes

**Polarity / ALT convention.** Which allele is "derived" can flip between ARG
draws, so we track the **ALT-allele frequency** $\varphi_{\text{alt}}$
consistently. With $c_{\text{alt}}$ = ALT count in the panel and, for draw $g$, the
ancestral base from the polarity table,

$$
\varphi_{\text{alt}}^{(g)}(T) =
\begin{cases}
\bar p(T \mid d_0=c_{\text{alt}},\, t_i^{(g)}), & \text{ALT derived in draw }g,\\[3pt]
1 - \bar p(T \mid d_0=n-c_{\text{alt}},\, t_i^{(g)}), & \text{ALT ancestral in draw }g,
\end{cases}
\tag{10}
$$

and, for the diploid likelihood, the matching **second** moment of the same ALT
frequency — where the ALT-ancestral branch must transform *both* moments, since
$\mathbb{E}[(1-X)^2] = 1-2\,\mathbb{E}[X]+\mathbb{E}[X^2]$:

$$
\varphi^{(2,g)}_{\text{alt}}(T) =
\begin{cases}
\bar p^{(2)}(T \mid c_{\text{alt}},\, t_i^{(g)}), & \text{ALT derived},\\[3pt]
1 - 2\,\bar p(T \mid n-c_{\text{alt}},\, t_i^{(g)})
  + \bar p^{(2)}(T \mid n-c_{\text{alt}},\, t_i^{(g)}), & \text{ALT ancestral}.
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
\bar p\big(T \mid d_0,\,[\text{below}_g,\text{above}_g]\big)
= \frac{1}{\text{above}_g-\text{below}_g}
  \int_{\text{below}_g}^{\text{above}_g}\! \bar p(T\mid d_0, t)\,dt .
$$

*Between* draws, averaging over the $M$ posterior draws (eq. 11) integrates the
remaining posterior on $t_i$ — each draw places the mutation on a different branch.
So $\bar\varphi_{\text{alt}}$ marginalises the full ARG posterior over the mutation
age, within-branch and between-draw. Averaging over draws,

$$
\bar\varphi_{\text{alt}}(T) = \frac{1}{M}\sum_g \varphi_{\text{alt}}^{(g)}(T),
\qquad
\bar\varphi^{(2)}_{\text{alt}}(T) = \frac{1}{M}\sum_g \varphi^{(2,g)}_{\text{alt}}(T), \tag{11}
$$

and use $\bar p_i\!\to\!\bar\varphi_{\text{alt}}(T)$ and
$\bar p^{(2)}_i\!\to\!\bar\varphi^{(2)}_{\text{alt}}(T)$ in (3a), with $a_i,c_i$
the sample's ALT dosage and ploidy at the site (ALT carriage is
polarity-independent). Both averages are over the *same* mixture, so the mixture
preserves $\bar\varphi^{\,2}\le\bar\varphi^{(2)}\le\bar\varphi$, and the
`--ploidy 1` path never touches the second moment.

**Chromosomes.** ARG draws are sampled independently per chromosome, and the $C$
chromosomes are unlinked — hence conditionally independent given $T$ — so the joint
likelihood factorises and the per-sample log-likelihood is the sum of the
per-chromosome terms:

$$
\log \mathcal{L}(T) = \sum_{c=1}^{C} \log \mathcal{L}_{c}(T), \tag{12}
$$

where $\mathcal{L}_{c}$ is the within-chromosome likelihood (eq. 3 over that
chromosome's sites, with its own ARG draws averaged in eq. 11). Each run computes
one $\log\mathcal{L}_{c}$; the merge step sums them over the $C$ chromosomes. (Here
$c$ indexes chromosomes, distinct from the per-site allele count $c_i$.) Sites
monomorphic in the panel ($c_{\text{alt}}\in\{0,n\}$) carry no trajectory and are
skipped.

---

## 7. Assumptions and caveats

1. **Neutrality.** Eq. (4) is the *neutral* trajectory — right for the vast
   majority of sites; strongly selected sites are not neutral. Guard with an
   approximately-neutral `--include-positions` set if selection is a concern.
2. **Ascertainment caps power** (§2): expect broad posteriors on array-ascertained
   SNPs.
3. **The ancient sample must not be in the ARG/ascertainment panels**, or the
   independence behind the ancestral-allele term breaks.
4. **Hardy–Weinberg** within a diploid individual assumes no recent inbreeding —
   the two alleles are independent *given* the population frequency, which is why
   (3b) needs $\mathbb{E}[r^2]$ rather than $\mathbb{E}[r]^2$; for
   haploid/pseudo-haploid data this is moot ($c_i=1$).
5. **Effective-size convention.** The curve is a diploid effective size
   ($N_e=1/(2\,\text{rate})$),
   so a single new copy has frequency $1/(2N_e)$ and drift variance $X(1-X)/(2N_e)$.
6. **Only two moments.** Eq. (3a) is exact for $c_i\le 2$. A likelihood over three
   or more alleles drawn from the same latent frequency would need moments of $X$ up
   to that order; the table carries only the first two.

---

## Glossary of variables

| symbol | meaning |
|---|---|
| $T$ | age of the ancient sample, generations before present (inferred) |
| $t_i$ | age of the mutation at site $i$ (interval $[\text{below},\text{above}]$ per draw, from the store) |
| $n$ | number of ARG-panel haplotypes |
| $C$ | number of chromosomes (independent given $T$) |
| $d_0$ | present count of the derived allele among the $n$ panel haplotypes |
| $c_{\text{alt}}$ | present count of the ALT allele among the $n$ panel haplotypes (from the panel VCF) |
| $x_0$ | present *population* derived frequency (unobserved; $d_0\sim\text{Binom}(n,x_0)$) |
| $a_i$ | number of derived (ALT) alleles the ancient sample shows at site $i$ |
| $c_i$ | number of called alleles at site $i$ (ploidy: 1 haploid/pseudo-haploid, 2 diploid, 0 missing) |
| $X_i(T),\,X(T)$ | derived-allele population frequency at time $T$ (a random variable) |
| $\bar p(T\mid d_0,t_i)$ | expected frequency at $T$ given present count and age — the table value (eq. 4, 9) |
| $\bar p^{(2)}(T\mid d_0,t_i)$ | expected *squared* frequency — the second table plane (eq. 4, 9a) |
| $r_i(T)$ | per-allele probability of *observing* derived given the frequency, with error $\varepsilon$ (eq. 2) |
| $\varphi_{\text{alt}}^{(g)},\,\bar\varphi_{\text{alt}}$ | ALT-allele frequency at $T$ in draw $g$, and its average over draws |
| $\varphi^{(2,g)}_{\text{alt}},\,\bar\varphi^{(2)}_{\text{alt}}$ | the corresponding second moments (eq. 10a, 11) |
| $N_e(t)$ | diploid effective population size at time $t$ ($N_e=1/(2\,\text{rate})$) |
| $\tau(t)$ | diffusion time, $\int_0^t dt'/(2N_e(t'))$ |
| $\tau_T,\tau_i$ | $\tau(T)$ and $\tau(t_i)$ |
| $\varepsilon_0$ | frequency of a single new copy at origin, $1/(2N_e(t_i))$ |
| $M_k(\tau)$ | $k$-th moment $\mathbb{E}[X(\tau)^k]$ of the neutral diffusion |
| $B$ | generator of the closed moment recursion (eq. 6) |
| $C(\Delta)=e^{B\Delta}$ | conditional-moment map over diffusion-time $\Delta$ |
| $\varepsilon$ | symmetric per-allele genotype-error probability |
| $\ell_i(T),\,\mathcal L(T)$ | per-site and total likelihood |
| $g,\,M$ | ARG posterior draw index and number of draws |

---

## References

- Kimura, M. (1955). *Solution of a process of random genetic drift with a
  continuous model.* PNAS 41:144–150.
- Griffiths, R.C. (2003). *The frequency spectrum of a mutation, and its age, in a
  general diffusion model.* Theoretical Population Biology 64:241–251.
- Ewens, W.J. (2009). *Mathematical Population Genetics: Introduction to the
  Stochastic Theory* (lecture notes, Guanajuato; abstracted from Ewens 2004,
  Springer).
- Deng, Y., Nielsen, R. & Song, Y.S. (2025). *Robust and accurate Bayesian
  inference of genome-wide genealogies for hundreds of genomes.* Nature Genetics
  (doi:10.1038/s41588-025-02317-9); software: github.com/popgenmethods/SINGER.
- bis101 & Ross-Ibarra, J. (2026). *RILAB/argtest: ARGtest v1.10.* Zenodo
  (doi:10.5281/zenodo.21865676).
- normalizeTEs pipeline — github.com/rossibarra/normalizeTEs.
