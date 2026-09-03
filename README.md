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
in the 26-haplotype panel** — exactly, from the neutral Wright–Fisher moment
recursion under your $N_e(t)$ curve. Genotypes enter by allele dosage (ploidy-aware;
see `--ploidy`); diploid genotype probabilities are nonlinear in the latent
frequency, so they use the conditional **second** moment $E[p_T^2]$ as well, from
the same recursion. Multiply across sites, average over ARG draws, sum across
chromosomes → a posterior over $T$ per sample.

---

## Pipeline

```
  Ne(t) curve ─► [1] precompute ─► freq_table.npz ┐
                                                   ├─► [2] infer (per chrom) ─► [3] merge ─► ages_table.tsv
  store + polarity + panel VCF + ancient VCF ──────┘
```

1. **Precompute** the frequency table `E[p_T | d0, t_i]` (and its second-moment
   plane `E[p_T^2 | d0, t_i]`) once (demography-specific, independent of
   samples/sites).
2. **Infer** per chromosome: look up the table by `(t_i, d0)` for every
   site/draw, form the per-sample likelihood.
3. **Merge** chromosomes into genome-wide per-sample posteriors.

---

## Scripts

| file | role |
|---|---|
| `precompute_freq_trajectory_moments.py` | build `freq_table.npz` = `E[p_T \| d0, t_i]` + `E[p_T^2 \| d0, t_i]` (exact moment recursion) |
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
| `--panel-vcf` | VCF of the **26 ARG-panel haplotypes** (gives `d0` = ALT count) | the panel the ARG was built on |
| `--vcf` | the **multi-sample ancient VCF** (all ancient samples, e.g. 430) | your aDNA calls |
| `--include-positions` *(optional)* | `chrom pos` site list (e.g. an approximately-neutral set) | your QC |

Notes:
- The **panel VCF** and the **ancient VCF** are different files. The panel gives
  present allele counts in the 26; the ancient VCF gives each sample's genotype.
- Chromosome labels must match across the store, both VCFs, and the ARG.
- The model is **neutral** (see MATH.md §7); restrict to a neutral site set with
  `--include-positions` if selection is a worry.

---

## Dependencies

- Python with `numpy` and `scipy` (the precompute uses `scipy.linalg.expm`);
  `matplotlib` optional (cohort plot).
- `normalize_tes` importable (the inference adapter layer calls
  `open_snp_age_store`, `open_draw_polarity`, `read_vcf_chunks`,
  `resolve_native_position_requests`). See the **ADAPTER LAYER** banner at the top
  of `posterior_sample_age_infer.py` — verify those names against your repo.

---

## Run it

### 1. Precompute the table (once)

```bash
python precompute_freq_trajectory_moments.py \
    --ne coalescence-ne-estimates.tsv \
    --n-sample 26 \
    --t-min 0 --t-max 30000 --n-t 300 \
    --age-min 10 --age-max 4e7 --n-age 100 \
    --output freq_table.npz
```

or `NE=... OUT=freq_table.npz sbatch slurm/run_precompute.sbatch`.

The `--t-*` grid is the **sample-age grid** and becomes THE grid the whole
analysis uses (inference reads it back from the table). Set `--t-max` above any
plausible sample age; `--age-*` should span your store's mutation ages.

### 2. Infer, per chromosome (all samples at once)

```bash
python posterior_sample_age_infer.py \
    --freq-table freq_table.npz \
    --store STORE --draw-polarity POLARITY \
    --panel-vcf panel26.vcf.gz \
    --vcf ancient.vcf.gz \
    --chrom chr1 \
    --ploidy 1 \
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
  cutoff are discarded; intervals crossing it are truncated. The conversion uses
  the `age_tau` axis stored in the frequency table, where
  $\tau(t)=\int_0^t ds/(2N_e(s))$. For constant $N_e=10{,}000$, $\tau=3$ is about
  60,000 generations. This is a
  numerical-reliability cutoff, not a claim that every older mutation is biologically
  uninformative.
- `--include-positions` — restrict to a QC'd / approximately-neutral site set.
- `--prior-file` — `T density` prior (two columns), interpolated onto the grid;
  default uniform.
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
- Confirm the panel VCF has the 26 haplotypes fully called at the sites you use
  (partially-called sites are skipped when forming `d0`).
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
