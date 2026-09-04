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
