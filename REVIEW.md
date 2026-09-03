# Review of `sample_age_dating`

## Overall assessment

The neutral moment-recursion precomputation is a promising foundation, especially
for haploid or pseudo-haploid inference. The resulting posteriors should not yet be
used for biological conclusions, however: the inference stage marginalizes ARG
draws in the wrong order, diploid likelihoods require more than the conditional
mean frequency, and several implementation bugs can produce incorrect results or
runtime failures.

## Major findings

### 1. ARG draws are marginalized in the wrong order

`posterior_sample_age_infer.py` averages allele frequencies over ARG draws
separately at each site and then multiplies the site likelihoods. An ARG draw is a
joint chromosome-wide genealogy, so the appropriate chromosome likelihood is
generally

$$
\mathcal L_c(T)=\frac{1}{M}\sum_g\prod_i \ell_{ig}(T),
$$

not

$$
\prod_i \ell_i\left(T;\frac{1}{M}\sum_g\phi_{ig}(T)\right).
$$

The current calculation destroys the across-site dependence within ARG draws and
will generally understate uncertainty. This affects both haploid and diploid
inference. Equations 11--12 of `MATH.md` encode the same ordering problem.

### 2. Diploid likelihoods cannot use only the expected frequency

Equation 3 substitutes $\bar p=E[X]$ into a binomial genotype likelihood. Diploid
genotype probabilities are nonlinear in the latent population frequency and need
second moments:

$$
P(AA)=E[q(X)^2], \qquad P(Aa)=2E[q(X)(1-q(X))].
$$

In general, $E[X^2]\ne E[X]^2$. The conditional first moment is sufficient for one
haploid Bernoulli observation, but not for `--ploidy 2`. The moment machinery could
be extended to tabulate the needed conditional second moments.

### 3. The ascertainment-independence claim is too strong

`MATH.md` states that ascertainment in a separate panel does not bias the estimate.
Separate individuals do not make ascertainment independent of the underlying
frequency trajectory. Inclusion because a SNP was observed or polymorphic in
another panel carries information about population frequency and potentially about
$T$. The model should condition on the ascertainment scheme or clearly describe
its omission as an approximation.

### 4. Age interpolation violates the mutation-existence boundary

`phi_lookup()` interpolates between precomputed mutation-age rows without
reapplying $p_T=0$ wherever $T\ge t_i$. For example, interpolation at mutation age
20 can return a positive frequency at sample age 50, although the mutation did not
yet exist. After interpolation, each result should explicitly be masked using the
actual interpolated mutation age.

### 5. Two adapter fallbacks fail because their defaults are evaluated eagerly

The nested `getattr` expressions in `_chunk_sites()` and `_resolve_rows()` evaluate
their fallback arguments even when the preferred attribute exists. Thus an object
with `chromosomes` but no `chrom`, or `rows` but no `row_indices`, can still raise
`AttributeError`. These should use explicit `hasattr` branches.

### 6. REF/ALT orientation is not checked between VCFs

Panel counts and ancient calls are joined by position alone. If the two VCFs have
swapped or otherwise different REF/ALT representations, the panel ALT count is
interpreted in the wrong orientation. The panel reader should retain REF and ALT,
and inference should verify or explicitly harmonize alleles before using a site.

## Additional code concerns

- In pseudo-haploid mode, every heterozygous call is converted to ALT by
  `alt_ct >= 1`. Heterozygous diploid-encoded calls should generally be treated as
  missing unless pseudo-haploid allele sampling happened upstream.
- Mutation-age intervals outside the table are silently clipped and renormalized,
  changing the assumed age distribution instead of reporting insufficient table
  coverage.
- Multiple interval records in one draw are collapsed to their overall minimum and
  maximum, potentially filling gaps and applying incorrect weighting.
- The demography parser assumes contiguous, non-overlapping windows but does not
  validate them. Gaps can cause negative within-window elapsed times.
- `--epsilon` should be checked, at least for $0\le\varepsilon<0.5$.
- Prior input should be checked for two columns, sorted ages, finite values, and
  nonnegative density.
- Chromosome merging checks sample order but not the grids or array dimensions.
- Both SLURM scripts source `slurm/slurm_conda_bootstrap.sh`, which is absent from
  the supplied archive, and suppress the resulting error with `|| true`.
- The archive has no automated unit or integration tests.

## Mathematics assessment

The following components appear correct under the stated neutral diffusion model:

- the neutral moment recursion in equation 6;
- the binomial expansion in equation 7;
- the bridge identity in equation 8, given the stated time orientation;
- the sampling-weighted ratio in equation 9 for the conditional first moment; and
- the neutral diffusion-time transformation for time-varying $N_e$, provided the
  demographic step function is valid.

The word "exact" should be qualified as exact within the neutral diffusion model.
The complete implementation also uses a finite table, interpolation, numerical
quadrature, a plug-in conditional expectation, and currently incorrect draw
marginalization.

## README assessment

The README is well organized and provides a useful overview, input table, example
commands, and output descriptions. It should nevertheless be revised alongside the
model and implementation:

- specify the correct order of site likelihood calculation and ARG-draw averaging;
- avoid presenting diploid support as complete until conditional second moments are
  implemented;
- qualify repeated claims that the method is exact;
- document the treatment or approximation of ascertainment;
- require matching or harmonized REF/ALT representations;
- explain how heterozygous calls behave in pseudo-haploid mode;
- document or reject mutation ages outside the precomputed grid; and
- include the referenced Conda bootstrap script or remove that dependency.

## Verification performed

The Python sources pass syntax compilation. A targeted interpolation check
reproduced the positive-frequency result across the $T\ge t_i$ boundary. No source
files were changed as part of this review.
