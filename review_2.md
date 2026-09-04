# Review 2: documentation, mathematics, and implementation consistency

## Scope and review state

This is a read-only review of the current working tree, focused on:

1. whether `README.md` functions as a concise how-to;
2. whether the README is clear and matches the code;
3. whether the derivation in `MATH.md` is correct; and
4. whether the implementation computes what `MATH.md` says it computes.

No existing file was changed for this review. The tree reviewed included the
uncommitted work for judgement call 9.6. The focused numerical-guard suite passes,
but the complete suite currently reports 36 passed and 4 failed (details below).

## Overall assessment

The first- and second-moment derivations in `MATH.md` are internally coherent under
the stated neutral Wright--Fisher diffusion, and the corresponding contractions in
`MomentEngine` match those equations. The haploid and diploid observation formulas
also agree with the implementation.

The largest remaining problem is not the moment engine. It is the treatment of ARG
draws and linked sites. The program computes a product of independently
draw-averaged site likelihoods. A SINGER draw is a joint chromosome-wide ARG, so a
fully marginalized chromosome likelihood instead averages a product of site
likelihoods within each draw. Moreover, multiplying linked-site likelihoods as if the
sites were conditionally independent is at best a composite likelihood. Therefore
the reported normalized curve should not currently be described without
qualification as a posterior with calibrated credible intervals.

The README is useful operationally, but it has accumulated derivation, validation
history, and modelling caveats that belong in `MATH.md` or `NOTES.md`. It also
contains several statements that no longer match the partial-panel implementation.

## Findings, ordered by importance

### 1. Critical: ARG draws are marginalized in the wrong order

`posterior_sample_age_infer.py:376-412` averages `phi` and `phi2` over ARG draws for
each site. It then converts those averaged moments into a site likelihood and adds
that site's log-likelihood to the chromosome total (`:413-439`). Across sites, this
is equivalent to

$$
\prod_i \left[\frac{1}{M}\sum_g \ell_{ig}(T)\right].
$$

For chromosome-wide ARG posterior draws, the marginal chromosome likelihood is
instead

$$
\frac{1}{M}\sum_g \prod_i \ell_{ig}(T).
$$

These are generally unequal. The current ordering discards the correlation of
genealogies and mutation ages across sites within the same ARG draw. It can produce
a curve that is too concentrated and can change its shape, not merely its scale.

`MATH.md:344-373` accurately describes what the code currently does (average each
site over draws, then multiply sites), but that is not the fully marginalized ARG
likelihood. This issue was already identified in `REVIEW.md` and remains unresolved.

Required decision: either restructure inference to retain a likelihood vector per
ARG draw until all sites on the chromosome have been accumulated, then use a
log-sum-exp over draws; or explicitly adopt and justify the current site-wise
posterior-marginal approximation in `NOTES.md` and stop presenting it as exact ARG
marginalization.

### 2. Critical: the across-site product is an independence/composite-likelihood assumption

`MATH.md:28-29` defines the likelihood as a product over sites, and `:362-374`
factorizes only at the chromosome level. Within a chromosome, however, sites are
linked and share a genealogy. Even after correcting the ARG-draw averaging order,
the product of per-site allele-observation probabilities is not automatically the
full joint data likelihood.

If conditional independence across sites is an intended working approximation,
the result is a composite likelihood or pseudo-posterior. That assumption and its
effect on uncertainty calibration need to be stated in `MATH.md` and justified in
`NOTES.md`. In particular, nominal 95% intervals should not be presumed to have 95%
coverage without calibration or thinning/blocking of sites.

### 3. High: inference renormalizes over surviving draws, contrary to equation 11

Equation 11 averages over all $M$ ARG draws. The implementation increments `cnt`
only for draws with usable polarity, age support, table coverage, and finite moments,
then divides by `cnt` (`posterior_sample_age_infer.py:376-412`). Consequently it
computes an average conditional on a draw surviving all filters, not the stated
$M^{-1}\sum_g$.

This distinction matters when failures correlate with mutation age, frequency, or
topology. The age cutoff and numerical guard are expressly age/frequency dependent.
Also, `sites_age_filtered` and `sites_numerical_failure` are incremented only when
no draw survives; partial loss of draws is invisible in `run.json`.

The implementation should either drop a site if any required draw is invalid, or
document and diagnose the conditional-on-retention mixture. At minimum, record
numbers/fractions of draws rejected for each reason.

### 4. High: judgement call 9.2 is only partially solved at the lower table boundary

Upper coverage is now guarded through `age_tau` and `--mutation-age-max`, but
`phi_lookup()` still silently clips mutation intervals below `age[0]` at
`posterior_sample_age_infer.py:201-204`. With the documented default
`--age-min 10`, a mutation interval lying partly or wholly below 10 generations is
renormalized at the boundary. A wholly younger interval becomes the age-10 endpoint.

This is the same change-of-distribution failure that judgement call 9.2 described,
just at the lower rather than upper endpoint. Either require the table to cover the
youngest retained mutation age, or reject/count out-of-range intervals. The README
currently says `--age-*` “should span” the store (`README.md:137-139`), but the code
does not enforce the lower side.

### 5. High: the README contradicts partial-panel support

The implementation now supports exact called-panel sizes from `--min-n` through
`--n-sample`, selecting a separate moment plane for each $n$. Nevertheless:

- `README.md:20-21` says conditioning is in “the 26-haplotype panel” without noting
  that the called count can be 20--26.
- `README.md:38-42` describes lookup only by `(t_i,d_0)`, omitting $n$.
- `README.md:267-268` says the panel must be fully called and partially called sites
  are skipped. This is false: sites with at least `--min-n` calls are retained.
- `MATH.md:256-258` still describes tables only over $(d_0,t_i,T)$, omitting the new
  $n$ axis.
- The glossary (`MATH.md:408-412`) describes $n$ as the number of panel haplotypes,
  but in the likelihood it is now the number called at that site.

The README should tell users to inspect `sites_panel_below_min_n`, and the math
should state the assumption that panel missingness is ignorable conditional on the
called subset. If missingness depends on allele state, using the $n$-specific table
does not remove bias.

### 6. High: merge validates arrays but not model settings

The new 9.6 checks correctly reject sample-order, grid, and likelihood-shape
mismatches (`posterior_sample_age_infer.py:468-485`). However, `run.json` records
only `epsilon`, `chrom`, and `mutation_age_max` (`:523-525`), omitting at least
`ploidy`, `min_n`, prior provenance, frequency-table identity/version, panel VCF,
and ancient VCF. `merge()` does not read `run.json` at all.

Thus identically shaped parts generated with different ploidy, error rates, panel
thresholds, priors, or frequency tables can still be summed silently. The prior is
applied again only after merge, which is good, but per-part `ages_table.tsv` may have
used a different prior and the provenance does not make that distinction clear.

### 7. Medium: diploid probability clipping can hide invariant failures

The moment engine now fails loudly when raw moments are unreliable, but downstream
inference independently clips `phi`, `phi2`, `E[r]`, `E[r^2]`, and each of the three
diploid genotype probabilities (`posterior_sample_age_infer.py:395-433`). If an
interpolation or polarity transformation violates

$$
E[r]^2\le E[r^2]\le E[r],
$$

the heterozygote probability can be negative and is silently changed to `1e-300`.
After independent clipping, the three values need not sum to one, whereas equation
3b says that they do.

Small roundoff can be tolerated, but material violations should be counted and
rejected rather than repaired silently. Tests should cover the invariants after age
interpolation, draw averaging, and ALT/REF transformation, not only inside
`MomentEngine`.

### 8. Medium: the table cutoff is described as an exact demographic translation, but is interpolated

Inference converts $\tau=3$ back to generations with a linear interpolation between
the log-spaced `(age_tau, age)` table nodes (`posterior_sample_age_infer.py:292-297`).
For a piecewise-constant $N_e$, $\tau(t)$ is piecewise linear, but a demographic
window boundary need not coincide with either neighboring age-grid node. The inverse
obtained this way is therefore approximate. `README.md:215-220` and
`MATH.md:275-280` read as though the actual demographic curve is used directly.

This is unlikely to dominate the inference, but the documentation should call it an
interpolated conversion, or the exact cutoff age should be stored during
precomputation.

### 9. Medium: “exact” is used too broadly

The neutral moment recursion and algebraic conditioning identities are exact within
the diffusion model. The pipeline is not numerically or statistically exact: it
uses float32 table storage, a finite log-age grid, log-age interpolation, 16-point
trapezoidal branch quadrature, a hard mutation-age truncation, rejected/renormalized
draws, and a product over linked sites.

Examples needing qualification include `README.md:20-25`, the script-table entry at
`:51`, and `precompute_freq_trajectory_moments.py:7-12`. “Exact moment recursion” is
appropriate; “the expected trajectory is computed exactly” is not.

### 10. Medium: precomputation performance claims do not match the implementation

`MATH.md:205-208` suggests one small matrix exponential, and `:288` calls the matrix
exponential cost negligible. In code, `Emoments()` recomputes three exponentials for
every `(n,d0,t_i,T)` call. `build_table()` loops over seven $n$ values, all $d_0$,
all mutation ages, and all sample ages (`precompute_freq_trajectory_moments.py:252-274`).
Many exponentials are identical across $d_0$ and could be reused.

The default grid implies millions of small `expm` calls. This is an implementation
performance issue rather than a mathematical error, but the README/SLURM claim that
precomputation is “cheap” should be verified against an actual default build.

### 11. Medium: test suite is not currently green

The full test run reports 36 passed and 4 failed:

- two failures in `test_second_moment_index_shift_would_be_caught`, because the
  deliberately shifted slice has the wrong length before it can test the intended
  scientific error;
- `test_frequency_decays_toward_the_origin`, whose finite offset from the origin is
  asserted to equal the single-copy initial frequency more tightly than the bridge
  conditioning appears to justify; and
- `test_plug_in_squared_mean_is_materially_wrong`, where the observed relative gap
  is about 24%, below an asserted 30% threshold.

These failures do not presently demonstrate errors in equations 6--9a, but a failing
default suite weakens every validation claim in the README. The tests should be
corrected or the implementation should be corrected, depending on an independent
review of each expected value. The new variable-$n$ table build also lacks an
end-to-end test that creates, saves, reloads, and uses a small multi-$n$ table.

### 12. Low/medium: some error-model wording exceeds what the code represents

`MATH.md:40-42` describes $\varepsilon$ as covering recurrent mutation and
mis-polarisation as well as VCF genotype error. The code implements a symmetric,
independent REF/ALT observation flip for each ancient allele. Recurrent mutation is
a change to the evolutionary process, and polarity uncertainty is handled through
ARG draws; neither is generally equivalent to this observation model. The narrower
description at `MATH.md:50-52` and in the README is the accurate one.

## README as a how-to

### Material that belongs in the README

- a two- or three-sentence purpose statement;
- prerequisites and the `normalizeTEs`/`PYTHONPATH` setup;
- required input schemas and coordinate/allele compatibility requirements;
- runnable local and SLURM commands;
- concise option descriptions, defaults, and units;
- outputs and the counters users must inspect; and
- a short “before trusting a run” checklist.

### Material better moved to `MATH.md` or `NOTES.md`

- the first/second-moment explanation in “What it does” and `--ploidy` beyond the
  operational fact that diploid mode needs a current table;
- the detailed numerical-cancellation discussion in “Validation provenance”;
- why the $\tau=3$ cutoff is scientifically conservative;
- the full pseudo-haploid heterozygote failure explanation (keep a one-line warning
  and link to `TODO.md`/`NOTES.md`);
- the strand-flip reasoning (keep only the required same-strand check, counter, and
  link to `NOTES.md`); and
- claims about validation experiments, which fit better in a validation section of
  `MATH.md` or a dedicated `VALIDATION.md`.

The long stale-checkout explanation at `README.md:90-117` is operationally useful,
but would scan better as a compact troubleshooting subsection after the primary run
instructions. It currently delays the first runnable command.

### Additional README/code mismatches

- `cohort_ages.png` is listed as an output (`README.md:191`) even though matplotlib
  is not in `environment.yml` and plotting failures are silently suppressed. It is
  only an optional output.
- `run.json` is described generically, but it does not record all settings and its
  site counters are global accepted-site counters, not per-sample usable-site counts.
- The SLURM wrappers do not expose `MIN_N`, `MUTATION_AGE_MAX`, or `PLOIDY` environment
  variables, even though these are central documented options. Direct CLI use works;
  the advertised SLURM route silently uses Python defaults.
- The README does not give a command for running the test suite or disclose that the
  default suite currently fails.

## Mathematics assessment

### Correct components

Subject to the neutral diffusion and observation-model assumptions, the following
are correct and reflected in code:

- the diffusion-time transformation in equation 5;
- the neutral moment recursion in equation 6;
- the binomial expansion in equation 7;
- the bridge contraction in equation 8;
- the conditional first- and second-moment ratios in equations 9 and 9a;
- the mutation-existence boundary $T\ge t_i$;
- the affine ancient-allele error transform in equations 2 and 3a;
- the haploid and diploid genotype probabilities in equations 3 and 3b;
- the ALT/REF first- and second-moment transformations in equations 10 and 10a; and
- summing log-likelihoods across genuinely independent chromosomes, provided the
  within-chromosome likelihood being summed is itself correctly constructed.

### Qualifications or corrections needed

- Equations 11--12 need to distinguish site-wise draw averaging from joint
  chromosome-wide ARG marginalization.
- The site product needs an explicit conditional-independence or composite-likelihood
  status.
- Equation 11 needs to reflect rejection and renormalization of draws, or the code
  must use all $M$ draws as written.
- The varying called-panel size $n$ and ignorable-missingness assumption need to be
  propagated consistently through equations 7--10 and the glossary.
- Numerical interpolation, quadrature, float32 storage, and truncation must qualify
  claims of exactness.
- The error model should not claim to represent recurrent mutation or generic
  mis-polarisation.

## Recommended order of work

1. Decide whether the target is a full chromosome-level ARG-marginal likelihood or
   an explicitly acknowledged composite/site-wise approximation.
2. Make draw rejection behavior consistent with equation 11 and expose partial draw
   losses in diagnostics.
3. Finish age-table coverage validation at the lower boundary.
4. Enforce merge provenance/settings compatibility.
5. Replace downstream probability clipping with invariant checks and loud rejection.
6. Make the existing test suite green and add a saved-table-to-inference integration
   test for variable $n$.
7. Then streamline the README and update `MATH.md`; documentation edits made before
   decisions 1--2 are likely to need rewriting.

## Bottom line

The mathematical core for conditional frequency moments is credible and mostly
implemented as written. The current end-to-end object is not yet demonstrably the
fully marginalized posterior claimed by the README because ARG draws and linked
sites are combined using a site-wise product-of-averages approximation. That issue
should be settled before biological interpretation or interval calibration. The
README should then be shortened into an operational guide, with derivations in
`MATH.md` and modelling judgments in `NOTES.md`.
