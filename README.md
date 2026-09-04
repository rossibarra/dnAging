# sample_age_dating

Estimate the age of ancient samples from the derived alleles they carry, using
SINGER ARG mutation ages plus a demography-aware, age-conditioned allele-frequency
model. Produces a per-sample posterior over age and a summary table.

The statistical and population-genetic derivation — the model and how we compute
it, with numbered equations — is in **[MATH.md](./MATH.md)**. Read that for the
model; this file is how to run it. Modelling judgement calls — approximations we
make deliberately and the conditions they rely on — are in
**[NOTES.md](./NOTES.md)**.

---

## What it does (one paragraph)

For each ancient sample and each ascertained site, carrying the derived allele
has probability equal to the derived-allele **frequency at the sample's age**
$T$ (and not carrying it, $1-\text{that}$). We compute that frequency as its
**expectation conditioned on the mutation's age $t_i$ and its present count $d_0$
in the ARG panel** — counted over the $n$ panel haplotypes **called at that site**,
anywhere from `--min-n` to 26 and not necessarily all 26 — from the neutral
Wright–Fisher moment recursion under your $N_e(t)$ curve (exact in the diffusion
model; the table you run against is a numerical tabulation of it — MATH.md §5).
Genotypes enter by allele dosage (ploidy-aware; see `--ploidy`); diploid genotype
probabilities are nonlinear in the latent frequency, so they use the
conditional **second** moment $E[p_T^2]$ as well, from the same recursion. Multiply
across sites, average over ARG draws, sum across chromosomes → a posterior over $T$
per sample.

---

## Pipeline

```
  Ne(t) curve ─► [1] precompute ─► freq_table.npz ┐
                                                   ├─► [2] infer (per chrom) ─► [3] merge ─► ages_table.tsv
  store + polarity + panel VCF + ancient VCF ──────┘
```

1. **Precompute** the frequency table `E[p_T | d0, t_i, n]` (and its second-moment
   plane `E[p_T^2 | d0, t_i, n]`) once, with a separate plane for every called-panel
   size `n` (demography-specific, independent of samples/sites).
2. **Infer** per chromosome: look up the table by `(t_i, d0, n)` for every
   site/draw — `n` = the panel haplotypes called at that site — form the per-sample
   likelihood.
3. **Merge** chromosomes into genome-wide per-sample posteriors.

---

## Scripts

| file | role |
|---|---|
| `precompute_freq_trajectory_moments.py` | build the first- and second-moment planes in `freq_table.npz`, one per called-panel size $n$ |
| `posterior_sample_age_infer.py` | per-chromosome inference for all samples; also does the merge |
| `slurm/run_precompute.sbatch` | STEP 1 as a batch job |
| `slurm/run_infer.sbatch` | STEP 2 (array over chromosomes) + STEP 3 (merge) |

---

## Inputs you need

| input | what it is | where it comes from |
|---|---|---|
| `--ne` | a **piecewise-constant `N_e(t)`** as a TSV of time windows (columns `time_left`, `time_right`, `effective_population_size`; optional `series`) | any demographic inference in that form — e.g. ARGtest `coalescence_ne_plots_from_ts.py --num_bins ~50` |
| `--store` | interval store (`snp-age-interval-v1`) giving `t_i` per site/draw | your normalizeTEs build |
| `--draw-polarity` | per-draw polarity table (ancestral base per site×draw) | `build_draw_polarity` |
| `--panel-vcf` | VCF of the **26 ARG-panel haplotypes** (gives `d0` = ALT count, and the called count `n` at each site) | the panel the ARG was built on |
| `--vcf` | the **multi-sample ancient VCF** (all ancient samples, e.g. 430) | your aDNA calls |
| `--include-positions` *(optional)* | `chrom pos` site list (e.g. an approximately-neutral set) | your QC |

Notes:
- The **panel VCF** and the **ancient VCF** are different files. The panel gives
  present allele counts in the panel; the ancient VCF gives each sample's genotype.
- Chromosome labels must match across the store, both VCFs, and the ARG.
- The model is **neutral** (see MATH.md §7); restrict to a neutral site set with
  `--include-positions` if selection is a worry.
- The across-site product is a PRF-style **composite likelihood**: it retains all
  quality-controlled SNPs and does not pretend that local LD is absent. Point
  estimates use the one-site marginal model, while posterior intervals are nominal
  unless calibrated by genome-scale simulation or a linkage-aware block bootstrap.
  See MATH.md §6.

---

## Dependencies

The conda environment is declared in [environment.yml](./environment.yml):

```bash
conda env create -f environment.yml     # once
conda activate dnaging
```

That covers `numpy`, `scipy` (the precompute uses `scipy.linalg.expm`), plus
`pytest` and `mpmath` for the test suite. `matplotlib` is optional (cohort plot).

**`normalize_tes` is deliberately not in that environment**, because it is not a
package — it is the [normalizeTEs](https://github.com/rossibarra/normalizeTEs)
checkout, which has no `pyproject.toml`/`setup.py` and so cannot be pip-installed.
It has to be put on `PYTHONPATH`:

```bash
export PYTHONPATH=/path/to/normalizeTEs:${PYTHONPATH:-}
```

Only the **inference** step needs it; precompute does not, and neither does the test
suite. The adapter layer expects `normalize_tes.snp_age_store`,
`normalize_tes.build_draw_polarity`, `normalize_tes.individual_age_spectrum` and
`normalize_tes.snp_position_resolution`.

All four exist as a `normalize_tes/` package on `normalizeTEs` `main`, along with
the symbols the adapter uses (`open_snp_age_store`, `open_draw_polarity`, `NO_CALL`,
`read_vcf_chunks`, `resolve_native_position_requests`).

> **Make sure the checkout is current.** `git status` reports agreement with the last
> *fetched* `origin/main`, so a stale clone can look up to date while still having
> the old flat layout (top-level `snp_age_store.py` and no `normalize_tes` package),
> against which every adapter import fails. Run `git fetch` and confirm
> `normalize_tes/__init__.py` exists before blaming the pipeline.

**On the cluster**, activate the environment (and export `PYTHONPATH`) *before*
`sbatch`: SLURM defaults to `--export=ALL`, so the job inherits both. Each sbatch
script preflights exactly what it imports and exits 90 with an actionable message
if the environment is missing, rather than failing later inside Python.

---

## Run it

### 1. Precompute the table (once)

```bash
python precompute_freq_trajectory_moments.py \
    --ne coalescence-ne-estimates.tsv \
    --n-sample 26 \
    --min-n 20 \
    --t-min 0 --t-max 30000 --n-t 300 \
    --age-min 10 --age-max 4e7 --n-age 100 \
    --output freq_table.npz
```

or `NE=... OUT=freq_table.npz sbatch slurm/run_precompute.sbatch`.

The `--t-*` grid is the **sample-age grid** and becomes THE grid the whole
analysis uses (inference reads it back from the table). Set `--t-max` above any
plausible sample age; `--age-*` should span your store's mutation ages.

It runs once per $N_e(t)$ curve and is reused for every chromosome and sample. The
default grid takes roughly an hour on one core; refining `--n-age`, `--n-t`, or the
panel-size span increases that cost. See MATH.md §5 for the computational details.

### 2. Infer, per chromosome (all samples at once)

```bash
python posterior_sample_age_infer.py \
    --freq-table freq_table.npz \
    --store STORE --draw-polarity POLARITY \
    --panel-vcf panel26.vcf.gz \
    --vcf ancient.vcf.gz \
    --chrom chr1 \
    --ploidy 1 \
    --min-n 20 \
    --mutation-age-max 3 \
    --epsilon 0.01 \
    --output results/Tage/chr1
```

As a SLURM array over `chroms.txt`:

```bash
TABLE=freq_table.npz STORE=... POLARITY=... PANELVCF=panel26.vcf.gz \
VCF=ancient.vcf.gz OUTROOT=results/Tage \
  sbatch --array=0-9 slurm/run_infer.sbatch
```

### 3. Merge chromosomes → genome-wide posteriors

```bash
python posterior_sample_age_infer.py \
    --freq-table freq_table.npz \
    --merge results/Tage/chr* \
    --output results/Tage/genome
```

or `MERGE=1 TABLE=freq_table.npz OUTROOT=results/Tage sbatch slurm/run_infer.sbatch`.

Merge validates that every part has the same sample order and the exact sample-age
grid from the frequency table, and that every likelihood matrix has shape
`(n_samples, n_grid)`. Mismatches stop with an error rather than being broadcast or
summed across different ages.

---

## Outputs (in each `--output/`)

| file | contents |
|---|---|
| `ages_table.tsv` | **the result** — one row per sample: `map_T, mean_T, median_T, ci95_lower_T, ci95_upper_T` (ARG generations) |
| `ll_marginal.npy` | `(N_samples × grid)` per-sample log-likelihood; rows align with `samples.txt` |
| `grid.npy`, `samples.txt` | the T grid and sample order |
| `run.json` | site/sample counts and settings |
| `cohort_ages.png` | histogram of per-sample MAP ages |
| `posterior/<sample>.tsv` | full posterior curve per sample (only with `--per-sample-tsv`) |

Ages are in **ARG generations**; convert to years with your generation time.

---

## Key options

- `--ploidy` — ploidy of the **ancient** genotypes: `1` = haploid / pseudo-haploid
  (one allele per called site; homozygous calls collapsed — the right choice for
  pseudo-haploid aDNA, and the default), `2` = true diploid genotypes (ALT dosage
  0/1/2; keeps het-vs-homozygote information). Using `2` on pseudo-haploid data
  written as homozygous diploid would double-count every site. **`1` assumes the
  ancient calls contain no true heterozygotes** — see the sanity check below.
  `2` builds the three genotype probabilities from the first **and second**
  conditional moments — Hardy–Weinberg holds only *given* the latent frequency, and
  $E[p_T^2]\neq E[p_T]^2$ — so it requires a table built by the current precompute
  script (it carries the `table2` plane); an older table makes `--ploidy 2` exit
  with a message rather than silently substituting the squared mean. `--ploidy 1`
  needs only the first moment and works with either table.
- `--epsilon` — symmetric **per-allele ancient-VCF genotype-error** probability,
  default `0.01`. It must satisfy $0\le\varepsilon<0.5$ and models VCF call error,
  not ARG uncertainty.
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
and match Kimura's constant-$N_e$ limit. See MATH.md §5.

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
lower `--age-min` if younger samples matter. See [NOTES.md](./NOTES.md).

If the interval store reports more than one branch interval for a mutation in any
ARG draw, that multiply mapped mutation is excluded completely. The number excluded
is reported as `sites_multiple_mapped` in `run.json`.

## Sanity checks before trusting results

- Confirm the ancient samples are **not** in the ARG/ascertainment panels.
- **The ancient genotypes are assumed pseudo-haploid** under the default
  `--ploidy 1`: one allele sampled per site and written as a homozygous diploid
  call, so the ALT dosage is 0 or 2 and never 1. Under that assumption the haploid
  collapse (`alt_ct >= 1` → derived) is exact. If the calls are genuinely diploid
  and contain heterozygotes, `--ploidy 1` promotes **every het to a derived
  observation**, inflating derived carriage and biasing $\hat T$ — silently, with no
  counter to reveal it. Use `--ploidy 2` for true diploid genotypes. Het-aware
  handling for the haploid path is deferred; see [TODO.md](./TODO.md).
- **The panel VCF does not have to be fully called**, so check
  `sites_panel_below_min_n` in `run.json`. `d0` is formed from however many of the 26
  haplotypes are called at a site, and the moment plane for that exact count is used;
  only sites with fewer than `--min-n` calls are skipped, and they are counted there.
  A large count means many sites carry too few called panel haplotypes — either the
  panel VCF is poorly called over your site set, or `--min-n` is set too high for it.
- **Confirm the panel and ancient VCFs are on the same strand**, then check
  `sites_allele_mismatch` in `run.json`. The two files are joined on position; a site
  whose REF/ALT are exactly swapped between them is harmonised
  (`c_alt → n-c_alt`), and one whose alleles cannot be matched is **skipped** and
  counted there. That check compares bases only and never complements them, so a
  large count is the signature of a **strand** disagreement, not just a different
  reference — and if strands do disagree, the A/T and G/C sites that *did* match are
  silently mis-oriented rather than skipped. See [NOTES.md](./NOTES.md).
- The $N_e$ TSV windows must **tile** the time axis: precompute exits if consecutive
  windows leave a gap or overlap, since the diffusion-time integral assumes
  contiguity.
- Expect **broad** posteriors — array ascertainment limits the age information
  (MATH.md §2, §7). A tight interval on a single sample deserves suspicion.
