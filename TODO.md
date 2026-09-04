# TODO

Deferred work, with enough context to pick up cold. Items graduate here from
JUDGEMENT_CALLS.md once a decision has been made about *what* to do but the work
itself is out of scope for the current pass.

---

## Heterozygote handling under `--ploidy 1`

**Decision made:** the ancient genotypes are *assumed* pseudo-haploid, and that
assumption is now documented in README.md. This item covers making the code robust
when the assumption does not hold.

**Current behaviour.** `posterior_sample_age_infer.py` builds the haploid
observation as

```python
a_eff = (alt_ct >= 1).astype(np.float64)   # any ALT allele -> derived
c_eff = (cl  >= 1).astype(np.float64)      # any called allele -> one observation
```

For genuinely pseudo-haploid data this is exactly right: one allele was sampled per
site and written as a homozygous diploid call, so `alt_ct` is 0 or 2 and never 1.
There are no heterozygotes to mishandle.

**The failure mode.** If true diploid calls containing heterozygotes are run through
`--ploidy 1`, every het is silently promoted to a derived observation. That inflates
derived-allele carriage at every heterozygous site, which biases $\hat T$ — and it
does so in the *default* mode, with no warning. The bias is first-order, not a
rounding concern.

**What to implement.**

1. A het policy for `--ploidy 1` instead of the unconditional collapse. At minimum
   `--het-policy missing` (treat `alt_ct == 1` as no-call, the safe default) and
   `--het-policy sample` (draw one of the two alleles, from a seeded RNG so runs are
   reproducible). The current behaviour would become `--het-policy derived`, kept
   only for backwards compatibility and documented as unsafe.
2. A `sites_het_collapsed` counter in `stats`, so `run.json` reveals when the
   assumption is being violated. Today this is completely invisible: the only
   symptom is a subtly biased posterior.
3. A loud warning (or a hard exit) when `--ploidy 1` encounters a non-trivial
   fraction of heterozygous calls, since that almost certainly means the data are
   not pseudo-haploid and the wrong mode was selected.

**Note** that `--ploidy 2` already handles true diploid genotypes correctly, using
the conditional second moment (`table2`). So the fix here is about making the
haploid path fail safely rather than about adding missing capability.

**Where:** `run_chromosome` in `posterior_sample_age_infer.py`; the `--ploidy`
documentation in README.md; MATH.md §2 if the policy changes the likelihood.

---

## Simulation validation with msprime

**Why the existing check is not enough.** `validate_moments_vs_mc.py` compares the
moment recursion against a forward simulation that takes Gaussian increments of
variance $p(1-p)/N_c$ per step — i.e. it simulates *the diffusion*. The analytics
assume the diffusion too, so the agreement confirms the **numerics** of the
recursion while leaving the **modelling approximation** untested. A discrepancy
between discrete Wright-Fisher and its diffusion limit would pass silently, and so
would an error in how we set up the conditioning.

**What msprime buys.** It can place samples at *different times*, so the quantity
the pipeline actually needs can be measured directly rather than reconstructed from
an intermediate moment. Simulate under the same step-function $N_e(t)$ with:

- $n$ modern samples at time 0 (the ARG panel), and
- one or more ancient lineages sampled at time $T$,

then place neutral mutations, and for each one read off its true age $t_i$ (the
branch it sits on), its derived count $d_0$ among the modern samples, and whether
each ancient lineage carries it. Binning by $(t_i, d_0)$ gives an empirical estimate
of exactly

$$
P(\text{ancient carries derived} \mid d_0, t_i) \;=\; E[p_T \mid d_0, t_i],
$$

which is the table's contents, measured end to end under the full coalescent.

**What it would validate that nothing currently does.**

1. The diffusion approximation itself, not merely its numerical evaluation.
2. The $T \ge t_i \Rightarrow p_T = 0$ boundary, as a simulated fact rather than an
   asserted one.
3. The conditioning setup — that conditioning on $d_0$ and $t_i$ is being done the
   way the derivation intends.
4. The second moment $E[p_T^2 \mid d_0, t_i]$, by sampling **two** ancient lineages
   at time $T$ and measuring the probability that *both* carry the derived allele.
   That is a direct, independent check on `table2` and hence on `--ploidy 2`.
5. As a by-product, the harness would give real (rather than synthetic) ARG draws,
   which is what REVIEW.md finding 1 — the ARG-draw marginalisation order — needs in
   order to be settled empirically.

**Cost.** Rare $(t_i, d_0)$ cells need many replicates, so this belongs in a
`@pytest.mark.slow` test or a standalone script rather than the default suite.
