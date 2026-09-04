# Open judgement calls — to discuss

Working document. Resolved items are pruned; their decisions stay in git history
(items 1-8 in `80b93f7`, later ones in the commit that closed each). Item numbers
are deliberately not renumbered.

---

## A. Likelihood composition — **OPEN, the only genuine decision left**

`review_2.md` findings 1, 2 and 3 are three faces of one question: how the
per-site, per-draw quantities are combined into a chromosome likelihood. They
should be decided together, because a fix to any one changes the others.

**1. Draws are marginalised in the wrong order.** An ARG draw is a joint,
chromosome-wide genealogy, so the chromosome likelihood is

$$
\mathcal L_c(T) = \frac{1}{M}\sum_g \prod_i \ell_{ig}(T),
$$

but inference averages $\phi$ across draws *per site* and then multiplies, i.e.
$\prod_i \ell_i\big(T; M^{-1}\sum_g \phi_{ig}(T)\big)$. For a haploid observation
$\ell$ is linear in $\phi$, so the per-site step is exact and the error is purely
the loss of across-site dependence within a draw.

**Direction of the error is NOT established.** REVIEW.md asserted it "understates
uncertainty" with no derivation. Averaging $\phi$ across draws *flattens* each
site's likelihood in $T$, which argues the opposite. This was being tested by
simulation when the run was lost to a rate limit; the result is unknown. Do not
quote a direction until it is measured.

**The correct form is not free.** $\sum_i \log \ell_{ig}$ is $O(n_{\text{sites}})$,
so a logsumexp mixture over draws collapses onto a single dominant draw — effective
sample size approaching 1. A naive fix trades bias for variance. `TODO.md`'s msprime
harness would supply real correlated draws to measure both.

**2. The across-site product is an unstated composite likelihood.** Eq. (1) takes
$\prod_i$ over ascertained sites, which assumes independence given $T$. Linked sites
are not independent. Unlike finding 1, this error has a *known* direction: a
composite likelihood over positively correlated data yields overly narrow
posteriors. This needs saying in MATH.md whatever is decided about the ordering.

**3. The average is over *surviving* draws, not $M$.** `cnt` counts only draws that
pass polarity, age-cutoff, table-coverage and numerical-stability filters, and the
result is divided by `cnt` — an expectation conditional on retention, not eq. (11)'s
$M^{-1}\sum_g$. This matters because the age cutoff and cancellation guard are
explicitly age- and frequency-dependent, so the filters correlate with the quantity
being averaged. Partial draw loss is also invisible: `sites_age_filtered` and
`sites_numerical_failure` increment only when *no* draw survives.

**The decision.** Either (a) accept all three as a documented composite/conditional
likelihood, stating the approximations and their known or unknown directions, or
(b) restructure to per-draw accumulation — which requires solving the ESS problem
and deciding whether a site with any invalid draw is dropped entirely. At minimum,
record per-reason draw-rejection fractions so the conditioning is measurable.
