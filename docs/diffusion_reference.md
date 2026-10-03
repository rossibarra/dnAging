# Diffusion pipeline reference

Detailed option semantics and validation provenance for the present-count
diffusion pipeline (approach 1). How to run it is in [README.md](../README.md);
the model is in [MATH.md](../MATH.md).

## Options in detail

- `--ploidy` — ploidy of the **ancient** genotypes: `1` = haploid / pseudo-haploid
  (one allele per called site; homozygous calls collapsed — the right choice for
  pseudo-haploid aDNA, and the default), `2` = true diploid genotypes (ALT dosage
  0/1/2; keeps het-vs-homozygote information). Using `2` on pseudo-haploid data
  written as homozygous diploid would double-count every site. **`1` assumes the
  ancient calls contain no true heterozygotes** — see the [README sanity checks](../README.md#sanity-checks-before-trusting-results).
  `2` builds the three genotype probabilities from the first **and second**
  conditional moments — Hardy–Weinberg holds only *given* the latent frequency, and
  $E[p_T^2]\neq E[p_T]^2$ — so it requires a table built by the current precompute
  script (it carries the `table2` plane); an older table makes `--ploidy 2` exit
  with a message rather than silently substituting the squared mean. `--ploidy 1`
  needs only the first moment and works with either table.
- `--epsilon` — symmetric **effective per-allele discordance** probability,
  default `0.01` in this inference CLI; it must satisfy
  $0\le\varepsilon<0.5$. With estimated ARGs it absorbs both ancient genotype
  error and residual ARG-induced incompatibility (mutation edge, age or polarity
  errors). These contributions are confounded: a fitted epsilon is **not** a
  measured genotyping-error rate. This is a robustness approximation, not an
  explicit ARG-error process, and does not replace averaging posterior ARG draws.
  The insertion driver uses the same interpretation; simulation wrappers may
  override the default, so check their submitted `EPSILON` and output `run.json`.
  Zero error can produce impossible observations on an estimated ARG; numerical
  clipping is not a scientifically modeled error floor. See [MATH.md](../MATH.md) eq. (2).
- `--mutation-age-max` — hard cutoff on mutation age in diffusion units, defaulting
  to $\tau=3$. Mutation-age intervals wholly above the corresponding generation-age
  cutoff are discarded; intervals crossing it are truncated. The corresponding
  generation age is interpolated from the table's demographic time axis. For
  constant $N_e=10{,}000$, $\tau=3$ is about 60,000 generations; see MATH.md §5 for
  numerical details. This is a
  numerical-reliability cutoff, not a claim that every older mutation is biologically
  uninformative.
- `--include-positions` — restrict to a QC'd / approximately-neutral site set.
- `--min-n` — minimum number of called ARG-panel haplotypes required at a site,
  default `20`. Sites with **at least** that many calls are used, not only fully
  called ones: precomputation builds a separate moment plane for every called-panel
  size from `--min-n` to `--n-sample`, inference uses the plane matching the site's
  exact called count, and sites below the threshold are skipped and reported as
  `sites_panel_below_min_n`. At inference it must be **at least** the `--min-n` the
  table was built with — a larger value simply leaves the lowest planes unused,
  while a smaller one exits with the missing panel sizes listed.
- `--prior-file` — `T density` prior with exactly two columns and at least two rows,
  interpolated onto the grid; ages must be finite, unique, and strictly increasing,
  while densities must be finite, non-negative, and not all zero. Default uniform.
- `--samples-file` — run a subset of the ancient samples.

---

## Validation provenance

Both moment planes of the table were checked against a forward Wright–Fisher Monte
Carlo (agreement to MC noise, including rare present-counts —
`validate_moments_vs_mc.py`), reproduce the $T \ge t_i \Rightarrow p_T=0$ boundary,
and match Kimura's limit at constant $N_e$. See [MATH.md](../MATH.md) §5.

The alternating conditioning sums become numerically unstable at large diffusion
times. Table construction measures cancellation for each moment and writes `NaN`
when only roughly 1–2 significant digits remain. Inference propagates that failure
and skips the affected draw (or site if no reliable draws remain), rather than
silently clipping a corrupted moment into the valid probability range. Inspect
`sites_numerical_failure` and `sites_age_filtered` in `run.json`. Newly built tables
must extend beyond $\tau=3$; precomputation exits if `--age-max` is too small, and
the inference step likewise rejects a table that does not cover its requested
cutoff.

Intervals reaching below `--age-min` are retained and counted in
`sites_age_clipped_low`. This is valid for sample ages at or above `--age-min`; use a
lower `--age-min` if younger samples matter. See [NOTES.md](../NOTES.md).

If the interval store reports more than one branch interval for a mutation in any
ARG draw, that multiply mapped mutation is excluded completely. The number excluded
is reported as `sites_multiple_mapped` in `run.json`.
