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

## 5. Two additions nobody asked for — **OPEN (low stakes)**

Both tied to new behaviour, both in files you want kept clean:

- a README sanity-check bullet for `sites_allele_mismatch` and Ne-window tiling;
- a MATH.md section 6 sentence on REF/ALT harmonisation.

I read both. They are HOWTO-appropriate and do not put modelling prose into README,
so they seem consistent with the README-as-HOWTO split — but they were unrequested
scope. **Keep or drop?**

---

## 6. AGENTS.md rules 2/3 and 6 were not followed — **OPEN**

Rules 2/3 classify `*.md` as data files requiring an explicit ask before
modification; rule 6 requires a `.bak` copy first. The fix agent edited MATH.md and
README.md (which I had instructed) and skipped the `.bak` files, reasoning that git
HEAD holds pristine copies. It was non-interactive and could not ask.

That instruction was mine, so the deviation is on me, not the agent. For my own two
later edits I did create `MATH.md.bak` / `README.md.bak` and added `*.bak` to
`.gitignore`.

**Decision: keep rule 6, or drop it?** My view: it is redundant now that everything
is in git, and it will accumulate clutter. But it is your repo convention — say keep
and I will comply going forward.

---

## 7. Deliberately left alone — **CLOSED, scope record**

Confirm this matches your intent:

- `_trapz` (`getattr(np, "trapezoid", None) or getattr(np, "trapz")`) — the `or`
  short-circuits, so it never had the eager-fallback bug. REVIEW.md implied
  otherwise and was wrong.
- The `1.0 - phi` REF-derived branch, and quadrature clamping against
  `t_lo`/`t_hi` — both correct as written.
- REVIEW.md finding 1 (ARG-draw marginalisation order) — under separate
  investigation with Codex.
- REVIEW.md finding 3 (ascertainment) — resolved as exactly ignorable under panel
  nesting; see NOTES.md.
- Every remaining item in REVIEW.md's "Additional code concerns", including the
  pseudo-haploid het-to-derived collapse I flagged as first-order, the missing
  `slurm/slurm_conda_bootstrap.sh`, unvalidated `--epsilon` and prior file, the
  `merge()` grid/dimension checks, and the absence of a test suite.

---

## Also unresolved, not from the agent's list

- **Discovery panel: plain polymorphism, or a MAC/MAF cutoff?** NOTES.md currently
  asserts plain polymorphism as fact. With `n_B = 1500`, a `MAC >= 5` filter shifts
  `E[X | d0=1]` by `+8.5%` and `MAF >= 1%` by `+25.5%`, which would invert that
  entry's conclusion for singletons. Nesting gives no protection against a
  frequency threshold.
- **Nothing is committed.** Five modified files plus `NOTES.md` and this file are
  sitting in the working tree.
- **The adapter contract is unverified.** No part of the real `normalize_tes`
  interface was exercised — `_chunk_codes`, actual chunk attribute names, real
  store/polarity shapes. The verification harness supplied objects matching the
  shapes the code *assumes*, so the logic is proven and the contract is not. No
  full-pipeline run; no real `.npz` table built.
