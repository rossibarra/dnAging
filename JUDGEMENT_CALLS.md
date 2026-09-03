# Open judgement calls — to discuss

Working document. Everything here came out of the first fix pass over REVIEW.md's
"fix these first" items (2026-09-03). Resolve these, then fold the outcomes into
MATH.md / NOTES.md / AGENTS.md and delete this file.

Status key: **OPEN** needs your call · **CLOSED** recorded, no decision needed

---

## 1. Untracked `NOTES.md` during the fix pass — **CLOSED**

The fix agent reported that my instructions described the tree as all-untracked when
it was actually committed and clean, and flagged an unexplained untracked `NOTES.md`
sitting in the repo.

Not a problem: that `NOTES.md` was mine, written mid-pass while the agent was
running. The agent correctly left it alone. No action.

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

## 4. One `phi_lookup` with a `key=` argument, not a separate `phi2_lookup` — **CLOSED unless you object**

The agent parameterised the existing lookup by table plane rather than duplicating
it, so the Fix-3 boundary mask and the branch quadrature cannot drift apart between
the first- and second-moment paths.

I reviewed this and think it is the right call — the duplication risk was real. No
action unless you disagree.

---

## 5. Two additions nobody asked for — **CLOSED, keep them**

Both tied to new behaviour, both in files you want kept clean:

- a README sanity-check bullet for `sites_allele_mismatch` and Ne-window tiling;
- a MATH.md section 6 sentence on REF/ALT harmonisation.

Both are correctly placed and stay. The README bullets are operational (check this
counter, tile your Ne windows), which is HOWTO. The MATH.md sentence sits directly
after eq. (10a), where $c_{\text{alt}}$ is defined, and states a precondition for
that equation to mean what it says — so it belongs with the derivation rather than in
NOTES.md, which is for approximations and their conditions. Reviewing this placement
is what surfaced item 8.

---

## 6. AGENTS.md rules 2/3 and 6 were not followed — **CLOSED, rules amended**

The fix agent edited MATH.md and README.md (which I had instructed) without the
permission request rules 2/3 required, and skipped rule 6's `.bak` copies. That
instruction was mine, so the deviation was on me, not the agent.

**Resolved by changing the rules rather than the behaviour.** AGENTS.md now says:

- markdown is documentation, not data, so modifying `*.md` needs no permission
  request (rule 3);
- but **recoverability still applies to markdown** — modify a file only when it is
  tracked-and-clean since the last commit, or after making a `.bak`. Not being a
  data file waives the permission request, not the safety net;
- `.bak` is required only when a file is not tracked-and-clean, since git otherwise
  holds the pristine copy (rule 6), and `*.bak` is gitignored.

Real data files (`*.vcf`, `*.tsv`, `*.npz`) are gitignored, hence untracked, so
"clean since the last commit" is never true for them and they still require a `.bak`
or an explicit ask. The relaxation only frees tracked text.

---

## 7. REVIEW.md findings that were wrong, or code already correct — **CLOSED**

No action needed on any of these, now or later.

- **`_trapz` was never buggy.** REVIEW.md's finding 5 named three eager-`getattr`
  sites; only two were real. `getattr(o, "a", getattr(o, "b"))` evaluates the
  fallback as a function *argument*, so it always runs — the bug. But
  `getattr(np, "trapezoid", None) or getattr(np, "trapz")` uses `or`, which
  short-circuits and never evaluates its right side when the left is truthy. Same
  idiom by eye, opposite behaviour.
- **The `1.0 - phi` REF-derived branch is exact.** When the ancestral allele is ALT,
  REF is derived, so the table is queried at `d0 = n - c_alt` and converted with
  `1 - phi`. That is exact by linearity, $E[1-X] = 1 - E[X]$, and it composes
  correctly with the existence boundary: at $T \ge t_i$ the lookup returns 0, so
  $1-0=1$, i.e. ALT frequency is 1 before the derived REF allele arose. Correct —
  everyone carried the ancestral (ALT) allele then.
- **Quadrature is correctly *not* clamped against `t_lo`/`t_hi`.** A partial nonzero
  value for a sample age inside the branch interval is the right answer, because the
  mutation age is uncertain within the branch: it may be older than $T$ (allele
  exists) or younger (it does not), and the average reflects that mixture. Only the
  per-node interpolation needed a mask, which is fix 3.
- REVIEW.md finding 1 (ARG-draw marginalisation order) — under separate
  investigation with Codex; a decision in its own right when it lands, not part of
  this pass.
- REVIEW.md finding 3 (ascertainment) — resolved as exactly ignorable under panel
  nesting; see NOTES.md.

---

## 8. Orientation check is strand-blind at palindromic sites — **CLOSED**

Found while reviewing item 5's placement. The REF/ALT harmonisation added in this
pass compares bases without complementing them, so at A/T and G/C sites a strand
flip is indistinguishable from an allele swap and $d_0$ can be silently inverted.
Non-palindromic sites fail safe (skipped, not inverted).

**Decided: documentation only, no code change.** Both VCFs are on the same strand
here, so the condition holds. Recorded in NOTES.md with the failure table and the
`sites_allele_mismatch` diagnostic; README gained a strand check pointing at it. No
`--exclude-palindromic` flag — `--include-positions` already covers the mitigation
if strand provenance ever becomes uncertain.

---

## 9. Deferred real work from REVIEW.md's "Additional code concerns" — **OPEN**

Promoted out of item 7, where it was wrongly filed as a closed scope record. None of
this is dismissed; it is simply out of scope for the first pass. Roughly in priority
order:

1. **Pseudo-haploid het-to-derived collapse — first-order bias, should be next.**
   `a_eff = (alt_ct >= 1)` maps *every* heterozygous call to derived, systematically
   inflating derived-allele carriage and biasing $\hat T$. This is in `--ploidy 1`,
   the default, and the mode REVIEW.md calls the trustworthy foundation. Hets should
   be treated as missing unless pseudo-haploid allele sampling already happened
   upstream.
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
9. **No automated test suite.** The first pass used throwaway verification scripts
   in a scratchpad; only the second-moment check in `validate_moments_vs_mc.py` is
   permanent.

---

## Also unresolved, not from the agent's list

- ~~**Discovery panel: plain polymorphism, or a MAC/MAF cutoff?**~~ **RESOLVED:**
  no filter. 1500 sequenced lines, and any SNP present in that panel was used. The
  panel also contains the ARG lines, so ascertainment is exactly ignorable and
  NOTES.md stands as written. (Had there been a cutoff, `MAC >= 5` would have shifted
  `E[X | d0=1]` by `+8.5%` and `MAF >= 1%` by `+25.5%` — nesting gives no protection
  against a frequency threshold, which is why the condition is recorded in NOTES.md
  for any future site set.)
- **Nothing is committed.** Five modified files plus `NOTES.md` and this file are
  sitting in the working tree.
- **The adapter contract is unverified.** No part of the real `normalize_tes`
  interface was exercised — `_chunk_codes`, actual chunk attribute names, real
  store/polarity shapes. The verification harness supplied objects matching the
  shapes the code *assumes*, so the logic is proven and the contract is not. No
  full-pipeline run; no real `.npz` table built.
