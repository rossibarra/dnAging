# sample_age_dating

Estimate the age of ancient samples from the derived alleles they carry, using
SINGER ARG mutation ages plus a demography-aware, age-conditioned allele-frequency
model. Produces a per-sample posterior over age and a summary table.

![Graphical abstract: modern panel and ARG plus ancient genotypes feed the
ancient-lineage insertion likelihood, which yields a per-sample age
posterior](docs/graphical_abstract.svg)

*Graphical abstract of the ancient-lineage insertion approach (MATH.md §8), the
preferred direction on the `insertion` branch. The validation figures are from
simulations; interval calibration under linkage is still open.*

## Status: which approach this README runs

| approach | MATH.md | status | where to look |
|---|---|---|---|
| 1. Present-count diffusion | §3–6 | **production** — the pipeline documented below | this file |
| 2. Tree-at-$T$ beta-binomial | §7 | archived methodological reference | [working_betabinom.md](./working_betabinom.md) |
| 3. Ancient-lineage insertion | §8 | preferred direction, in development; point-calibrated in simulation, intervals not yet calibrated under linkage | [working_insertion.md](./working_insertion.md) |

`posterior_sample_age_infer.py` stays the production entry point until insertion
consumes the production SNP/ARG store, handles posterior ARG draws and chromosomes,
and passes an end-to-end real-data-shaped simulation (see
[PROJECT_MAP.md](./PROJECT_MAP.md), "Branch and production policy"). Everything
from "What it does" onward describes approach 1.

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
| `pipeline/precompute_freq_trajectory_moments.py` | build the first- and second-moment planes in `freq_table.npz`, one per called-panel size $n$ |
| `pipeline/posterior_sample_age_infer.py` | per-chromosome inference for all samples; also does the merge |
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
- The model is **neutral** (see MATH.md §9); restrict to a neutral site set with
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
python pipeline/precompute_freq_trajectory_moments.py \
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
python pipeline/posterior_sample_age_infer.py \
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
python pipeline/posterior_sample_age_infer.py \
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

| option | default | short version |
|---|---|---|
| `--ploidy` | `1` | `1` = pseudo-haploid ancient calls (assumes no true hets); `2` = true diploid, needs a table with the second-moment plane |
| `--epsilon` | `0.01` | effective per-allele discordance, $0\le\varepsilon<0.5$; absorbs genotype error **and** residual ARG error, so it is not a genotyping-error rate |
| `--mutation-age-max` | $\tau=3$ | numerical-reliability cutoff on mutation age (diffusion units) |
| `--min-n` | `20` | minimum called panel haplotypes per site; must be ≥ the table's `--min-n` |
| `--include-positions` | all sites | restrict to a QC'd / approximately-neutral site set |
| `--prior-file` | uniform | two-column `T density` prior, interpolated onto the grid |
| `--samples-file` | all samples | run a subset of the ancient samples |

Full semantics, and the validation provenance (Monte Carlo checks, numerical-failure
handling, the `run.json` counters `sites_numerical_failure`, `sites_age_filtered`,
`sites_age_clipped_low`, `sites_multiple_mapped`), are in
[docs/diffusion_reference.md](./docs/diffusion_reference.md).

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
  (MATH.md §2, §9 and the ascertainment appendix). A tight interval on a single
  sample deserves suspicion.
