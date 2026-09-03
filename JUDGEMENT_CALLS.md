# Open judgement calls — to discuss

Working document. Everything here came out of the first fix pass over REVIEW.md's
"fix these first" items (2026-09-03).

Closed items have been pruned. The **full record of all resolved decisions — items
1, 4, 5, 6, 7 and 8 — is in commit `80b93f7`**; read that rather than reconstructing
them. Item numbers below are deliberately NOT renumbered, so the gaps are expected
and references from commit messages stay valid.

Resolve what remains, fold the outcomes into MATH.md / NOTES.md, and delete this
file.

---

## 2. Alternating-sum cancellation destroys the FIRST moment at large `tau_i` — **OPEN, highest priority**

Not one of the five requested fixes. Pre-existing, found by the agent, and I
confirmed it independently via the cancellation amplification
`A = sum|term| / |sum term|`, which needs no high-precision reference:

| `tau_i` | `d0=1` | `d0=8` | `d0=13` | sig. digits left |
|---|---|---|---|---|
| 0.5 | 3e7 | 2e10 | 8e9 | 5–8 |
| 3 | 7e10 | 3e13 | 7e12 | **2–3** |
| 5 | 6e11 | 3e14 | 6e13 | **1–2** |
| 10 | 9e13 | 6e15 | 8e15 | **0** |

At `tau_i = 10` the raw first moment returns `-0.908` (`d0=8`) and `-2.266`
(`d0=13`); `Efreq`'s long-standing `np.clip(..., 0, 1)` silently turns those into
`0.0`, where a high-precision reference gives ~`0.24` and ~`0.50`.

Why it matters: this is the **first** moment, i.e. the `--ploidy 1` haploid path —
the mode REVIEW.md calls the trustworthy foundation. It has never been loud about
failing. Even at `tau_i = 0.5` only 5–8 digits survive.

**Blast radius is unknown and I did not guess it.** It depends on your `Ne(t)` curve
and the mutation ages SINGER actually reports. `--age-max` defaults to `4e7`
generations, so for any plausible `Ne` the upper decades of that log-spaced axis sit
in the zero-digit regime — but whether real sites land there is an empirical
question. **Action: compute `tau_i = integral dt/(2 Ne)` for the real Ne file against
the real mutation-age distribution.** That decides urgency.

Options:
- **Guard (cheap, honest).** Compute the amplification per table entry, write `NaN`
  below ~2 surviving digits, and have `phi_lookup` propagate `NaN` to a skipped site
  instead of `nan_to_num`-ing it to `0.0`. Loses sites honestly rather than keeping
  them wrongly.
- **Real fix (redesign).** Reformulate to a positivity-preserving representation —
  evolve sample-count probabilities through the coalescent dual rather than expanding
  `Binom(d0; n, x)` into an alternating polynomial in `x`. No cancellation at all.
- **High precision at build time.** Viable in principle since precompute is one-time,
  but naive mpmath is too slow at 26 x 100 x 300 evaluations unless the
  lower-bidiagonal `expm(B*u)` closed form is used instead of a general `expm`.

My recommendation: guard now so nothing is silently wrong, then scope the
reformulation as its own piece of work.

---

## 3. `E[p^2]` clamped to `[E[p]^2, E[p]]` in `Emoments` — **OPEN**

The agent added this clamp. Both bounds are exact facts for `X` in `[0,1]`
(Cauchy–Schwarz, and `X^2 <= X`), and without them the cancellation in item 2 yields
*negative* genotype probabilities.

The cost: it converts a detectable failure into a silent one. At `tau_i = 8, d0 = 8`
it returns `E2 = E1 = 0.4605` where the truth is `0.297`.

The agent verified Cauchy–Schwarz on the **unclamped** ratios in a separate
reimplementation, so its test could not pass vacuously — that part is sound.

**Decision: keep the clamp, or replace it with the loud `NaN` guard from item 2?**
My recommendation: replace. Given the first moment is also garbage in that regime,
silence is the wrong failure mode. The clamp is only papering over item 2.

---

## 9. Deferred real work from REVIEW.md's "Additional code concerns" — **OPEN**

Promoted out of item 7, where it was wrongly filed as a closed scope record. None of
this is dismissed; it is simply out of scope for the first pass. Roughly in priority
order:

1. ~~Pseudo-haploid het-to-derived collapse~~ — **RESOLVED, moved to
   [TODO.md](./TODO.md).** Decided: the ancient genotypes are *assumed*
   pseudo-haploid, so `alt_ct` is 0 or 2 and the collapse is exact. The assumption is
   now documented in README.md (option text plus a sanity check). Making the haploid
   path fail safely when the assumption is violated — a het policy, a
   `sites_het_collapsed` counter, and a warning — is deferred to TODO.md.
2. **Mutation ages outside the table are silently clipped and renormalised**,
   changing the assumed age distribution instead of reporting insufficient table
   coverage. Interacts with item 2.
3. **Multiple interval records in one draw are collapsed to their overall min and
   max**, potentially filling gaps and misweighting the branch quadrature.
4. **`--epsilon` is unvalidated** — should require $0 \le \varepsilon < 0.5$.
5. **The prior file is unvalidated** — no check for two columns, sorted ages, finite
   values, or non-negative density.
6. **`merge()` checks sample order but not the $T$ grids or array dimensions**, so
   mismatched chromosome parts can be summed silently.
7. **`read_panel_alt` requires the panel fully called** (`tot_called == n_expected`)
   and drops any partially-called site without counting it in `stats`.
8. **Both SLURM scripts source a missing `slurm/slurm_conda_bootstrap.sh`** and
   suppress the failure with `|| true`.
9. **No automated test suite** — **IN PROGRESS.** A `tests/` suite is being built
   now (pytest, run under `/opt/anaconda3/bin/python`, with `normalize_tes` stubbed
   since it is not importable locally). Until it lands, only the second-moment check
   in `validate_moments_vs_mc.py` is permanent; the first fix pass relied on
   throwaway scratchpad scripts.
