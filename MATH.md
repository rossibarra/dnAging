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
**ascertained in a panel different from the one the ARG was built on** — that
ascertainment panel supplies only the *site positions* and is used nowhere else in
the model.

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
at time $T$, derived with probability $p_i(T)$ = the **derived-allele population
frequency** at time $T$ (§3). Adding a symmetric **per-allele** genotype-error
probability $\varepsilon$ (aDNA damage, sequencing/genotyping error, recurrent
mutation, mis-polarisation), the per-allele probability of *observing* the derived
state is

$$
q_i(T) = (1-\varepsilon)\,p_i(T) + \varepsilon\,\bigl(1-p_i(T)\bigr). \tag{2}
$$

Treating the sample's called alleles at a site as independent draws (Hardy–Weinberg
— i.e. no recent inbreeding within the individual), the observed derived count is
$\text{Binomial}(c_i, q_i)$, so (dropping the $T$-independent binomial coefficient)

$$
\ell_i(T) = q_i(T)^{\,a_i}\,\bigl(1-q_i(T)\bigr)^{\,c_i-a_i},
\qquad
\log\ell_i(T) = a_i \log q_i(T) + (c_i-a_i)\log\!\bigl(1-q_i(T)\bigr). \tag{3}
$$

**Ploidy.** $c_i$ is set by a `--ploidy` flag matching the ancient genotype calls:

- **Haploid / pseudo-haploid** ($c_i = 1$): one allele per called site, $a_i\in\{0,1\}$.
  A homozygous call is collapsed to a single observation (crucial for pseudo-haploid
  aDNA written as `0/0` or `1/1`, so it is not counted twice). Equation (3) reduces
  to $\ell_i = q_i$ if the allele is derived, $1-q_i$ if ancestral.
- **Diploid** ($c_i = 2$): the true genotype, $a_i\in\{0,1,2\}$, under Hardy–Weinberg:
  $\ell_i \propto (1-q_i)^2,\; 2q_i(1-q_i),\; q_i^2$ for a homozygous-ancestral,
  heterozygous, or homozygous-derived call respectively. Equivalently, the
  probability an individual **carries at least one** derived allele is
  $1-(1-p_i)^2$ and of being homozygous ancestral is $(1-p_i)^2$; using the full
  three-genotype form keeps the heterozygote-vs-homozygote information rather than
  collapsing to presence/absence.

A missing/uncalled site has $c_i = 0$ and contributes $\log\ell_i = 0$.

**Ascertainment.** A site absent from a sample's VCF is missing data, not an
observation of "does not carry." Because ascertainment depends on a separate panel
(not the ARG panel) and is independent of $T$ and of the ancient genotypes,
conditioning on the ascertained set does not bias $\hat T$; it only costs power —
and most where the lower-bound signal lives (young, rare-in-discovery alleles are
under-ascertained).

Everything now hinges on one quantity: $p_i(T)$, the derived-allele frequency at
time $T$.

---

## 3. What we compute: the age-conditioned expected frequency

We use the derived-allele frequency **conditioned on the two things known robustly
per site**: its present count in the panel, $d_0$, and its mutation age, $t_i$
(from the store). Conditioning on age encodes **survival to the present** — an
allele that is old and still segregating was, on average, at *higher* frequency in
the past than a young allele of the same present frequency (it had to be, to
persist):

$$
p_i(T) \;=\; \bar p(T \mid d_0, t_i) \;=\; \mathbb{E}\!\left[\,X(T)\,\middle|\, d_0,\, t_i\,\right], \tag{4}
$$

where $X(\cdot)$ is the neutral population-frequency trajectory. In simulation this
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
to order $n+1$, one $(n+2)\times(n+2)$ matrix exponential** — tiny for the modest
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

For $T \ge t_i$ (sample older than the mutation) the allele does not yet exist, so
$\bar p = 0$. This is exact and was validated against a forward Wright–Fisher Monte
Carlo — agreement to MC noise, including the rare $d_0=2$ bin (0.22 vs 0.23), the
$T\ge t_i$ boundary, and Kimura's constant-$N_e$ limit. The table
$\bar p(T\mid d_0,t_i)$ over $(d_0, t_i, T)$ is built once and shared across all
sites and all samples.

**Numerics / scaling with $n$.** The alternating-sign conditioning sums in (7) and
(9) are prone to catastrophic cancellation: the binomial weights grow like $2^{n}$
while the moments are small, so the sum loses on the order of $0.3\,n$ decimal
digits. Double precision is therefore comfortable for panels up to a few tens of
haplotypes (the ARG regime here), but the method degrades and eventually breaks
beyond roughly $n \approx 40$–$50$; there one must switch to extended precision
(e.g. `mpmath`) or reformulate the conditioning in a numerically stable basis
(orthogonal-polynomial / spectral moments rather than the raw power moments). The
matrix-exponential cost itself is negligible ($(n+2)^3$).

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

where $t_i^{(g)}$ is the mutation age in draw $g$. **This age is uncertain, and the
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
\bar\varphi_{\text{alt}}(T) = \frac{1}{M}\sum_g \varphi_{\text{alt}}^{(g)}(T), \tag{11}
$$

and use $p_i(T)\!\to\!\bar\varphi_{\text{alt}}(T)$ in $q_i$ (eq. 2), with $a_i,c_i$
the sample's ALT dosage and ploidy at the site (ALT carriage is
polarity-independent).

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
4. **Hardy–Weinberg** within a diploid individual (eq. 3) assumes no recent
   inbreeding; for haploid/pseudo-haploid data this is moot ($c_i=1$).
5. **Effective-size convention.** The curve is a diploid effective size
   ($N_e=1/(2\,\text{rate})$),
   so a single new copy has frequency $1/(2N_e)$ and drift variance $X(1-X)/(2N_e)$.

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
| $p_i(T),\,X(T)$ | derived-allele population frequency at time $T$ |
| $\bar p(T\mid d_0,t_i)$ | expected frequency at $T$ given present count and age — the table value (eq. 4, 9) |
| $q_i(T)$ | per-allele probability of *observing* derived, with error $\varepsilon$ (eq. 2) |
| $\varphi_{\text{alt}}^{(g)},\,\bar\varphi_{\text{alt}}$ | ALT-allele frequency at $T$ in draw $g$, and its average over draws |
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
