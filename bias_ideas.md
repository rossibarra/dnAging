# Possible Sources of Age-Estimation Bias

This document tracks possible causes of systematic bias in the sample-age
estimator and tests that can distinguish them. The current simulations use the
true modern ARG, infinite-sites mutations, no genotype error, 26 modern
haplotypes, and a 10 Mb region. Across the matched simulations, the estimator
underestimates age when all sites are retained. The bias depends on effective
population size.

## Current observations

| Effective population size | Signed bias (generations) | RMSE (generations) |
| ---: | ---: | ---: |
| 10,000 | -836 | 890 |
| 50,000 | -564 | 653 |
| 100,000 | -611 | 705 |
| 200,000 | -1,093 | 1,219 |
| 500,000 | -1,360 | 1,599 |

Removing all sites with modern derived-allele frequency below 5% removes only
derived singletons in this 26-haplotype panel. It shifts estimates upward,
eventually producing severe overestimation:

| Effective population size | Signed bias after removing singletons | RMSE |
| ---: | ---: | ---: |
| 10,000 | -542 | 623 |
| 50,000 | +495 | 644 |
| 100,000 | +1,056 | 1,132 |
| 200,000 | +1,618 | 1,879 |
| 500,000 | +5,205 | 5,451 |

Consequently, singleton sites exert strong net pressure toward younger ages.
A selective-filter test removed only low-frequency sites at which the ancient
sample carries the derived allele, while retaining singleton absences. This
shifted estimates sharply upward:

| Effective population size | Signed bias after removing carried singletons | RMSE |
| ---: | ---: | ---: |
| 10,000 | -471 | 564 |
| 50,000 | +1,138 | 1,186 |
| 100,000 | +2,320 | 2,427 |
| 200,000 | +4,099 | 4,442 |
| 500,000 | +7,283 | 7,398 |

At high Ne, many estimates reach the upper edge of the 15,000-generation grid,
so the reported upward bias is partly censored and may be an underestimate.
This experiment confirms that rare carried alleles provide most of the
downward likelihood pressure.

## 1. Disproportionate leverage from rare carried alleles

For ancient genotype indicator Y and modeled carrier probability p(T), a site
contributes

    Y log p(T) + (1 - Y) log(1 - p(T)).

A carried rare allele can have very large leverage because log p(T) changes
rapidly when p(T) is small. Although carried singletons are rare, their effect
may outweigh the much more numerous singleton absences. If the model predicts
too little ancient sharing for these sites, it will select a younger age to
increase p(T).

Predictions:

- Removing only carried singletons should move estimates substantially older.
- Removing only absent singletons should move estimates younger, or much less
  strongly older.
- Likelihood scores partitioned by derived count should show d0 = 1 pulling
  toward younger ages.

Tests:

1. Repeat while excluding only absent singletons.
2. Plot the summed likelihood and its derivative by d0 class.

## 2. Double conditioning on modern derived count

The observed ARG edge already determines the set, and therefore the number, of
modern descendants carrying a mutation. Conditional on one mutation occurring
on a known edge under a constant mutation rate, its position in generations is
uniform along that edge. The current ratio-of-integrals calculation additionally
weights candidate mutation ages by a term proportional to P(d0 | t). This may
condition on the modern count a second time and put the wrong weight on parts of
the edge.

If this weighting makes ancient sharing decline too quickly with sample age,
the estimator compensates by choosing an age that is too young.

Test:

- Replace the P(d0 | t)-weighted mutation-age integral with an integral uniform
  in generations along each observed mutation edge. Compare bias and RMSE on
  exactly the same simulations.

This is the leading structural hypothesis because it directly concerns what is
conditioned on once the ARG edge is observed.

## 3. Insufficient conditioning on the observed ARG

The likelihood groups sites primarily by modern derived count and mutation-age
interval. Two edges with the same descendant count can have different local
genealogical contexts. The exact edge, its descendant set, and the surrounding
tree may provide information about ancient sharing that is absent from an
allele-count-only diffusion calculation.

Tests:

- In simulations containing the ancient lineage, calculate the exact ancient
  sharing state implied by whether that lineage descends from the mutation.
- Compare empirical sharing probabilities among edges matched for d0 and age
  interval but separated by other tree features.
- Determine whether residuals cluster by terminal versus internal edges,
  branch length, tree height, or descendant topology.

## 4. Wrong diffusion conditioning or boundary behavior

The required quantity is a past-frequency distribution conditional on modern
ascertainment, the present descendant configuration, and mutation occurrence on
the observed edge. A diffusion initialized at 1/(2Ne) and later conditioned by
a binomial modern count may not equal that distribution.

We have evidence that the diffusion frequency estimator is biased upward at low
true allele frequencies. By itself, overestimating p(T) underestimates
1 - p(T) for ancient absences and should push age estimates older. It therefore
cannot directly explain the observed overall downward bias. However,
miscalibration may vary with T or differ between carried and absent rare alleles,
so its net likelihood effect must be measured rather than inferred from mean
frequency bias alone.

Tests:

- Compare modeled p(T) with empirical carrier probabilities from msprime,
  stratified by d0, true T, and mutation-edge age interval.
- Plot calibration residuals rather than only mean predicted versus observed
  frequency.
- Examine carried and absent sites separately.

## 5. Mutation-edge integration and numerical interpolation

Mutation ages are integrated across the child-to-parent interval of the
mutation edge. Incorrect quadrature, interpolation, or weighting within long
intervals could assign too much mass near young endpoints and pull inferred
ages downward.

Tests:

- For a subset of sites, compare the current interpolated integral with dense
  direct quadrature.
- Stratify likelihood scores by edge length and endpoint ages.
- Increase the interpolation resolution and check convergence of both MAP and
  posterior mean.

At Ne = 500,000, a handful of entire very-young moment-table rows were NaN from
numerical cancellation and had to be interpolated along the age axis. Results
at this Ne should be checked against a more stable table calculation.

## 6. Effective-population-size scaling or convention mismatch

A factor-of-two error between haploid and diploid diffusion time, or a mismatch
between msprime population size and the table's 2Ne scaling, would create a
systematic error that changes with Ne. The simulations specify diploid Ne and
the tables use time divided by 2Ne, which appears consistent, but this should be
verified analytically and with a limiting case.

Tests:

- Recompute one set using alternative Ne and time scalings without changing the
  data.
- Check whether bias collapses when expressed in diffusion units T/(2Ne).
- Verify every conversion between generations, diffusion time, mutation-age
  cutoff, and table coordinates.

## 7. Modern-polymorphism ascertainment mismatch

The simulations and analysis both retain sites polymorphic in the modern panel,
so gross ascertainment mismatch should be limited. A subtler mismatch can arise
if the theoretical conditioning treats a random mutation and then conditions on
d0, while the data consist of mutations discovered on already-observed ARG
edges.

Test:

- Derive the likelihood under the exact simulation ascertainment procedure and
  compare it term by term with the implementation. This overlaps with the
  double-conditioning hypothesis.

## 8. Composite-likelihood dependence

Linked sites are treated as independent contributions. Linkage makes posterior
intervals too narrow and gives genomic regions with many correlated mutations
too much apparent information. Dependence alone does not normally cause bias
when marginal probabilities are calibrated, but it can amplify bias from any
misspecified marginal likelihood.

Tests:

- Thin sites by genetic distance and compare point estimates.
- Use one mutation per tree or per recombination block.
- Compare ordinary posterior intervals with block-bootstrap uncertainty.

## Lower-priority explanations

- The 20-generation T grid is too fine to explain errors of hundreds or
  thousands of generations.
- The mutation-age cutoff removes very old mutations and could matter if its
  interaction with allele count is strong, but few sites are removed in these
  simulations.
- ARG inference error cannot explain the present results because the true modern
  ARG is used.
- Genotype error cannot explain them because these runs assume no error.

## Suggested order

1. Quantify likelihood pull by d0 and ancient carried/absent state.
2. Complete the carried- and absent-singleton exclusion tests.
3. Test uniform-in-generations integration along the observed ARG edge.
4. Calibrate p(T) against empirical msprime probabilities within d0 and edge-age
   strata.
5. Only then test thinning and other composite-likelihood corrections.
