# Open judgement calls — to discuss

Working document. Everything here came out of the first fix pass over REVIEW.md's
"fix these first" items (2026-09-03).

Closed items have been pruned. The **full record of resolved decisions 1, 4, 5, 6,
7 and 8 is in commit `80b93f7`**. Items 2 and 3 were resolved by replacing silent
moment clipping with a cancellation-based `NaN` guard and adding a default
`--mutation-age-max 3` cutoff in diffusion units. This is about 60,000 generations
at constant $N_e=10{,}000$; inference translates it through the actual $N_e(t)$
curve stored with the frequency table's age axis.

Resolve what remains, fold the outcomes into MATH.md / NOTES.md, and delete this
file.

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
