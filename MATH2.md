# The math, second approach: dating an ancient sample from the ARG directly

This note derives the likelihood implemented in [betabinom/](./betabinom/). It is an
alternative to [MATH.md](./MATH.md), not a revision of it. The two differ in what
they condition on:

| | MATH.md | this note |
|---|---|---|
| conditions on | present count $d_0$ and mutation age $t_i$ | the marginal tree at time $T$ |
| needs | $\mathbb{E}[p_T]$ from a two-time diffusion bridge | one moment vector per mutation age |
| cost | three matrix exponentials per table entry | one per age |
| demography | enters through $N_e(t)$ and the diffusion clock | enters through the tree, plus the clock |

Every symbol is defined in the glossary at the end; display equations are numbered
for reference. Section 7 states the assumptions, section 8 what is and is not
validated.

---

## 1. Setup and goal

We have an ancestral recombination graph on a panel of $n$ modern haplotypes, and
one **ancient sample** drawn from the same population at an unknown time $T$
generations before present. The ancient sample is *not* in the panel and does not
contribute to the ARG.

At each site the ARG gives, per marginal tree, the edge carrying the mutation and
the times of that edge's endpoints. Cutting the tree at $T$ gives two numbers that
are the whole basis of the method:

- $n_T$, the number of panel ancestral lineages alive at $T$
- $k$, how many of those descend from the mutation

We want the posterior over $T$:

$$
p(T \mid \text{data}) \propto p(T) \mathcal{L}(T). \tag{1}
$$

The essential simplification is that **conditioning on the tree makes the
present-day panel count redundant**. Given the ARG, the panel's allele
configuration is determined, so nothing is left to condition on except the
population frequency at $T$, which is what the ancient sample actually samples.

---

## 2. The per-site likelihood

The ancient calls are **pseudo-haploid**: one allele per called site. Let $X(T)$ be
the derived-allele population frequency at time $T$. A single lineage drawn from
that population carries the derived allele with probability $X(T)$. With a
symmetric per-allele error probability $\varepsilon$, the probability of
*observing* the derived state given the frequency is

$$
r(T) = \varepsilon + (1 - 2\varepsilon) X(T). \tag{2}
$$

Because (2) is affine in $X(T)$, the per-site likelihood needs only the **first
moment** of the frequency. Writing $p_i(T) = \mathbb{E}[X_i(T) \mid \text{ARG}]$
and $g_i$ for the observed allele,

$$
\ell_i(T) = \begin{cases}
\varepsilon + (1-2\varepsilon) p_i(T), & g_i \text{ derived} \cr
1 - \varepsilon - (1-2\varepsilon) p_i(T), & g_i \text{ ancestral}.
\end{cases} \tag{3}
$$

No second moment is required anywhere. A true-diploid path would need
$\mathbb{E}[X^2]$, which the same machinery supplies but which is not used here.

**Pseudo-haploid sampling needs no change to any of this.** A real ancient sample
is a diploid individual from which one allele is drawn at random per locus. That
individual carries two lineages, each an exchangeable draw from the population at
$T$, so given $X(T)=x$ each allele is Bernoulli in $x$; a random choice between
them is a mixture of two Bernoulli variables with the same parameter, which is
again Bernoulli in $x$. The marginal at every site is therefore identical to a
true haploid lineage and (2) and (3) stand unaltered.

It also means **Hardy-Weinberg is not assumed here**, unlike the diploid path of
MATH.md. The two alleles of one individual are never used at the same locus, so
whether they are independent given the frequency never arises, and recent
inbreeding in the ancient individual is harmless.

What does differ is dependence *across* sites. A single haploid genome is one
chromosome, so neighbouring sites share an ancestry; pseudo-haploid switches
chromosome at random between sites, so roughly half of adjacent pairs are drawn
from two independent lineages. Real pseudo-haploid data therefore carries **less**
linkage than the haploid genomes this model was validated on. Since the estimate
depends only on the per-site marginals, that is a variance effect and not a bias,
and the block bootstrap of section 6 absorbs it because it resamples the observed
data rather than assuming a correlation structure.

Everything therefore reduces to one quantity: $p_i(T)$, the expected
derived-allele frequency at $T$ given the tree.

---

## 3. The conditional frequency

Fix a site and a candidate $T$. Suppose the mutation is present at $T$ and carried
by $k$ of the $n_T$ panel lineages, and that it arose $a$ generations before $T$.

**Sampling.** The $n_T$ panel ancestral lineages at $T$ are an exchangeable sample
from the population alive at $T$, so given $X(T) = x$,

$$
P(k \mid n_T, x) = \binom{n_T}{k} x^k (1-x)^{n_T - k}. \tag{4}
$$

Equation (4) keeps only the dependence of the tree on $x$ through the *count*. It
is not the whole dependence: under the structured coalescent the mutant class
coalesces at rate proportional to $1/x$, so the branch lengths within that class
also carry information about the frequency (Griffiths 2003, section 4.3). That term
is dropped. Section 8 records that it was tested and found undetectable.

**Prior.** Before seeing the tree, $X(T)$ is the frequency of a mutation of age $a$
that began as a single copy. Call its density $f_a$; its moments
$M_m(\tau_a) = \mathbb{E}[X^m]$ come from the neutral recursion of section 5.

**Posterior.** Combining, the conditional density is proportional to
$x^k (1-x)^{n_T-k} f_a(x)$. Expanding the binomial,
$(1-x)^{n_T-k} = \sum_j \binom{n_T-k}{j} (-1)^j x^j$, both the numerator and the
normaliser become alternating sums of moments:

$$
\mathrm{num}(k, n_T, a) = \sum_{j=0}^{n_T-k} \binom{n_T-k}{j} (-1)^j M_{k+1+j}(\tau_a),
\qquad
\mathrm{den}(k, n_T, a) = \sum_{j=0}^{n_T-k} \binom{n_T-k}{j} (-1)^j M_{k+j}(\tau_a),
\tag{5}
$$

$$
\varphi(k, n_T, a) = \mathbb{E}\bigl[X(T) \mid k, n_T, a\bigr]
  = \frac{\mathrm{num}(k, n_T, a)}{\mathrm{den}(k, n_T, a)}. \tag{6}
$$

$\mathrm{den}$ has a meaning of its own that matters in section 4: it is
proportional to $P(k \mid n_T, a)$, the probability of the observed lineage count
given the mutation's age.

**Limits.** For $a \to 0$, all moments collapse to powers of the initial frequency
and $\varphi \to 1/(2N_e)$: a mutation born an instant before $T$ is one copy in
the whole population. For large $a$ the age dependence saturates and $\varphi$
approaches the beta-binomial predictive $k/(n_T+1)$, which is the posterior mean of
$\mathrm{Beta}(k, n_T-k+1)$ under the polarised neutral prior $f(x) \propto 1/x$.
That limit is accurate to 0.15% for edges comfortably older than $T$, but it is
wrong by an order of magnitude for young mutations, which is exactly why (6) is
needed rather than the closed form.

---

## 4. Marginalising the mutation age

Under the infinite-sites model a mutation falls uniformly along its edge, between
the child node at $t_c$ and the parent at $t_p$. Three cases arise.

**The edge is entirely older than $T$**, i.e. $t_c > T$. The mutation certainly
exists at $T$, and $k$ is the number of lineages at $T$ descending from the child
node, which does not depend on where on the edge the mutation sits.

**The edge straddles $T$**, i.e. $t_c \le T < t_p$. Only the portion above $T$ puts
the mutation in existence at $T$, and there it sits on exactly **one** lineage, so
$k = 1$. The portion below $T$ contributes zero.

**The edge is entirely younger than $T$**, i.e. $t_p \le T$. The mutation postdates
the sample and $p_i(T) = 0$ exactly.

The age must be marginalised with the correct weight. Conditioning on the ARG makes
the age posterior **non-uniform**, because an age that makes the observed lineage
count more likely is itself more likely:

$$
f(a \mid \text{ARG}) \propto f(a) \int f_a(x) x^k (1-x)^{n_T-k} dx
  = f(a) \cdot \mathrm{den}(k, n_T, a). \tag{7}
$$

So the marginal is a **ratio of integrals, not an integral of ratios**:

$$
p_i(T) = w \cdot
\frac{\int_{\max(T, t_c)}^{t_p} \mathrm{num}(k, n_T, t - T) dt}
     {\int_{\max(T, t_c)}^{t_p} \mathrm{den}(k, n_T, t - T) dt},
\qquad
w = \frac{t_p - \max(T, t_c)}{t_p - t_c}. \tag{8}
$$

The factor $w$ is the fraction of the edge lying above $T$; it equals 1 for an edge
entirely older than $T$ and handles the straddling case. Averaging $\varphi$
uniformly instead over-weights young ages and inflates $p$ by up to a factor of two
on long edges.

**Equation 8 is an approximation in one respect.** $w$ is the *prior* probability
that the mutation lies above $T$, whereas the reweighting argument of (7) calls for
the posterior,

$$
P(\text{above} \mid \text{ARG}) =
\frac{\int_{\text{above}} f \cdot \mathrm{den}}
     {\int_{\text{above}} f \cdot \mathrm{den} + W_{\text{below}}}, \tag{8a}
$$

and $W_{\text{below}}$ is not commensurable with $\mathrm{den}(a)$ in any obvious
way: below $T$ the allele does not exist, so there is no binomial factor at $T$ to
weigh against. The prior is used because the posterior is not well defined without
settling that comparison. Section 8 records the calibration evidence that the
difference is below the level simulation can resolve.

A note on what $w$ does and does not do. For a straddling edge, the placements
below $T$ make carriage impossible, so they reduce $p$; they do **not** zero the
site's likelihood. A site removes a value of $T$ outright only when every placement
is impossible *and* the sample carries the allele, which is the third case above.
That hard wall is where most of the information about $T$ lives, and $\varepsilon$
controls how hard it is: a violation costs $\log\varepsilon$, not $-\infty$.

---

## 5. The moments

The moments in (5) are those of the standard neutral diffusion. Measuring time in
units of $2N_e$ generations removes $N_e$ from the process, so with

$$
\tau(t) = \int_0^{t} \frac{dt'}{2 N_e(t')} \tag{9}
$$

the moments obey the closed linear system of the drift-free diffusion,

$$
\frac{dM_m}{d\tau} = \frac{m(m-1)}{2}\bigl(M_{m-1} - M_m\bigr),
\qquad m = 1, 2, \dots \tag{10}
$$

Writing $B$ for its lower-bidiagonal generator, $M(\tau) = e^{B\tau} M(0)$, and a
new mutation starts from a single copy, so $M_m(0) = \bigl(1/(2N_e)\bigr)^m$.

Equation (5) needs moments to order $n_T + 1$, so the matrix is about the same size
as MATH.md's, but **one** exponential per age suffices: there is no two-time bridge
and no joint moment $\mathbb{E}[X_T X_{\mathrm{pres}}^m]$, because conditioning at
$T$ on the tree has replaced conditioning at the present on $d_0$.

**Numerics.** The alternating sums in (5) cancel catastrophically, losing on the
order of $0.3 n_T$ decimal digits, and in float64 they return negative values or
NaN at large $\tau$. The table is therefore built in `mpmath` at 60 decimal digits
from an exact moment map and only then rounded to float64. Quadrature and
interpolation were checked against a 400-node, 240-age reference and agree to 0.2%
even at the smallest $p$.

---

## 6. The full likelihood

Sites are combined as a composite likelihood:

$$
\log \mathcal{L}(T) = \sum_i \log \ell_i(T). \tag{11}
$$

This treats linked sites as independent, which they are not. The consequence is
severe and measured: intervals taken from the curvature of (11) cover the truth
0/12, 1/10 and 5/10 across three independent simulation sets. **Any interval must
come from a block bootstrap** over contiguous genomic windows, which resamples
whole blocks and so preserves linkage within them. Block-bootstrap standard
deviations were checked against the actual scatter of estimates and are slightly
conservative, so they are honest about the noise even though (11) is not.

Where the ARG is a posterior sample rather than a single tree, draws are
marginalised as a mixture at the chromosome level, exactly as in
[MATH.md](./MATH.md): the site product is taken within a draw and the draws are
averaged only afterwards.

---

## 7. Assumptions and conditioning

**What is conditioned on.**

1. The marginal tree at each site, through $n_T$ and $k$ only.
2. The mutation's edge, with uniform placement along it.
3. The demography, through $\tau$ in (9) and through $1/(2N_e)$ as the initial
   frequency.

**What is deliberately not conditioned on.**

4. The present-day panel count $d_0$. Given the tree it is redundant; it is a
   deterministic function of the same information.
5. The mutant clade's internal coalescent times. Faster coalescence within the
   clade implies a lower frequency, so this is information in principle. It was
   tested on 862,343 observations and found to carry no detectable signal once
   $(k, n_T, a)$ is fixed.
6. The ancient lineage's own branch. Mutations private to the ancient sample are
   informative in simulation but **unavailable in practice**, because ancient
   samples are genotyped only at sites ascertained as polymorphic in a modern
   discovery panel, so a private mutation has no site to be called at.

**Assumptions.**

7. **Neutrality.** Equations (4) to (6) are the neutral diffusion and the neutral
   coalescent. Strongly selected sites violate both.
8. **Exchangeability of the panel lineages.** Equation (4) treats the $n_T$
   ancestral lineages at $T$ as a random sample from the population at $T$. This
   is a property of the neutral coalescent and fails under structure or selection.
9. **The ancient sample is an independent draw from the same population at $T$**,
   and is absent from the panel and from the ARG. If it were in the panel it would
   partly be conditioned on itself.
10. **Pseudo-haploid calls.** One allele per site, so only the first moment enters.
    Collapsing a homozygous diploid call to one observation is required, or each
    site is counted twice.
11. **Infinite sites**, so each site has one mutation and placement along an edge
    is uniform. Recurrent mutation is outside the model.
12. **Known demography.** $N_e(t)$ is taken as given. It enters both the clock (9)
    and, through the tree, the values of $n_T$.
13. **Known error rate.** $\varepsilon$ is fixed, not estimated. Section 8 records
    how much this matters.
14. **Upstream control of the aDNA-specific artefacts.** Three properties of real
    pseudo-haploid data are outside this model and are assumed to be handled before
    it sees the data. Each biases the sample *older* if left uncorrected, which is
    the same direction as the unexplained offset in section 8, so none of them is
    a safe omission.
    - *Reference bias*: reads carrying the non-reference base map less readily, so
      the drawn allele skews toward REF and the sample appears to carry fewer
      derived alleles. $\varepsilon$ cannot represent this, being symmetric by
      construction.
    - *Missing calls scored as non-carriage*: absence of a site means uncalled, not
      that the sample lacks the derived allele.
    - *Deamination*: C-to-T damage is strand- and context-dependent, so it is not
      the symmetric per-allele error that $\varepsilon$ models. End-trimming or
      restriction to transversions belongs upstream.

---

## 8. What is validated, and what is not

**Established by simulation**, against msprime with true ARGs:

- $p = 0$ for mutations postdating the sample is exact: 0 carriers in 41,395
  opportunities.
- The large-age limit $k/(n_T+1)$ is accurate to 0.15% for edges above $T$. The
  unpolarised alternative $k/n_T$ is wrong by 9.6%.
- Calibration of (8) is flat across strata: by predicted $p$, by fraction of edge
  above $T$, by clade size, and on a log scale down to $p \approx 10^{-3}$.
  Overall observed/predicted is 0.996 for straddling edges and 0.999 for edges
  above $T$, with the worst single bin at 1.2 standard deviations. The
  den-weighting of (8) is what makes this true: averaging uniformly instead gives
  ratios of 0.09, 0.17 and 0.31 across bins of the fraction of edge above $T$,
  against 1.12, 1.01 and 0.92 with it.
- The two approximations flagged in sections 3 and 4 -- dropping the
  structured-coalescent term in (4), and using the prior rather than the posterior
  in $w$ -- are therefore good to within what these tests can resolve. That is not
  the same as being exact, and this likelihood has already proved sensitive to
  sub-1% structure, so neither should be assumed harmless.
- Over 20 simulations at $N_e = 50{,}000$, 10 Mb, 26 haplotypes, with true ages
  from 454 to 9,293 generations, the estimate regresses on truth with **slope
  1.002** and correlation 0.960.

**Not established.**

- The fit carries a **constant additive offset of about +1340 generations**
  (95% CI 965 to 1738). Its origin is unknown. Six candidate mechanisms were
  tested and eliminated: clade shape, the straddling term, a second binomial on
  the clade size, the age-marginalisation order, panel ascertainment, and the
  multiallelic exclusion.
- Corrected leave-one-out, the estimator is unbiased with RMSE 967 and about 12%
  median relative error. But the correction is **fitted, not derived**, and was
  calibrated on true ARGs with known $N_e$ and $\mu$ at $\varepsilon = 10^{-6}$.
- **Sensitivity to $\varepsilon$ is strong.** Raising it from $10^{-6}$ to
  $10^{-3}$ on error-free data moved RMSE from 1169 to 1562 and bias from +886 to
  +1239; earlier runs moved by thousands of generations. Since $\varepsilon$ softens
  the hard constraint of section 4, and that constraint carries most of the
  information, a fixed $\varepsilon$ is a real modelling choice rather than a
  detail.
- Everything here uses **true** ARGs. Inferred trees will add error on top, and
  none of it has been measured.

---

## Glossary of variables

| symbol | meaning |
|---|---|
| $T$ | age of the ancient sample, generations before present (inferred) |
| $n$ | number of haplotypes in the modern panel |
| $n_T$ | panel ancestral lineages alive at time $T$ |
| $k$ | those lineages that descend from the mutation |
| $a$ | age of the mutation at time $T$, i.e. it arose $a$ generations before $T$ |
| $t_c, t_p$ | ages of the child and parent nodes of the mutation's edge |
| $w$ | fraction of that edge lying above $T$ |
| $X(T)$ | derived-allele population frequency at $T$, a random variable |
| $f_a$ | density of $X(T)$ for a mutation of age $a$ from a single copy |
| $M_m(\tau)$ | $m$-th moment of the neutral diffusion |
| $B$ | generator of the moment recursion, eq. 10 |
| $\varphi(k, n_T, a)$ | expected frequency at $T$ given the tree, eq. 6 |
| $p_i(T)$ | the same after marginalising the mutation age, eq. 8 |
| $r(T)$ | probability of observing the derived allele, eq. 2 |
| $\varepsilon$ | symmetric per-allele genotype-error probability |
| $g_i$ | observed allele of the ancient sample at site $i$ |
| $N_e$ | diploid effective population size |
| $\tau$ | diffusion time, eq. 9 |

---

## References

- Kimura, M. (1955). *Solution of a process of random genetic drift with a
  continuous model.* PNAS 41:144-150.
- Griffiths, R.C. (2003). *The frequency spectrum of a mutation, and its age, in a
  general diffusion model.* Theoretical Population Biology 64:241-251. Section
  3.3.1 gives the Dirichlet and Beta structure behind eq. 4; section 4.3 gives the
  structured coalescent within the mutant class.
- Ewens, W.J. (2009). *Mathematical Population Genetics: Introduction to the
  Stochastic Theory* (lecture notes, Guanajuato, March 2009).
  [PDF](http://www.rgwinther.com/Ewens2009MathematicalPopulationGeneticsTheGuanajuatoLectures.pdf)
- Deng, Y., Nielsen, R. & Song, Y.S. (2025). *Robust and accurate Bayesian
  inference of genome-wide genealogies for hundreds of genomes.* Nature Genetics
  57:2124-2135.
