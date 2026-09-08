# Codex handoff: pilot ancient-sample age inference on HPC

## Objective

Run a first real-data diffusion pilot for **50 ancient samples** using the current
`main` implementation. Start with the supplied SNP-age interval store (ideally
one ARG posterior draw), run one chromosome as a smoke test, then run all
chromosomes and merge the likelihoods into genome-wide age posteriors.

Do not change scientific defaults or source code merely to make the run pass.
Diagnose mismatched inputs explicitly. Keep all generated tables, logs, and
results under a new run directory. Do not overwrite an earlier run.

## Inputs the user will provide

Ask for or set these absolute paths:

```bash
PROJECT=/path/to/dnAging
NORMALIZE_TES=/path/to/normalizeTEs
NE_FILE=/path/to/coalescence_ne.tsv
STORE=/path/to/snp_interval_store
POLARITY=/path/to/draw_polarity
PANEL_VCF=/path/to/arg_panel.vcf.gz
ANCIENT_VCF=/path/to/ancient_50.vcf.gz
SAMPLES_FILE=/path/to/ancient_50.samples.txt
CHROMS_FILE=/path/to/chroms.txt
RUN_ROOT=/path/to/new/output/dnaging_pilot_YYYYMMDD
```

The user may initially mention only the Ne file, store, and ancient samples.
Inference also requires `POLARITY` and `PANEL_VCF`; locate the matching products
next to the store or ask for them. Never substitute the ancient VCF for the panel
VCF. `--n-sample` is the maximum number of **modern ARG-panel haplotypes**, not
the number of ancient individuals.

The samples file must contain one ancient VCF sample ID per line. Confirm that it
contains exactly 50 unique IDs and that every ID occurs in the ancient VCF.

## Non-negotiable scientific settings for this pilot

- Use the diffusion implementation on `main`.
- Use `--marginalise uniform`; denominator weighting performed worse in tests.
- Use `--ploidy 1` for pseudo-haploid calls, even if they are encoded as 0/0 and
  1/1. Stop and ask if true heterozygous genotypes are intended.
- Use `--mutation-age-max 3` diffusion units.
- Use `--min-n 20`, unless the modern panel has fewer than 20 haplotypes; do not
  lower it silently.
- For real data, begin with `--epsilon 0.01` as the documented per-allele error
  assumption. This value must be reported as an assumption and later tested for
  sensitivity; do not describe it as estimated.
- Use the exact table engine. Do **not** pass `--float64`.
- Drop multiply mapped/multiple-mutation sites through the store's eligibility
  policy; do not restore them.

## 1. Checkout and environment preflight

```bash
cd "$PROJECT"
git status --short
git rev-parse --abbrev-ref HEAD

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dnaging
export PYTHONPATH="$NORMALIZE_TES:${PYTHONPATH:-}"

python -c 'import numpy, scipy, normalize_tes; print("imports OK")'
test -f "$NORMALIZE_TES/normalize_tes/__init__.py"
mkdir -p "$RUN_ROOT/logs" "$RUN_ROOT/results"
```

Record the git commit, conda environment, paths, and command lines in
`$RUN_ROOT/provenance.txt`. Do not pull, switch branches, or modify tracked files
without checking with the user.

## 2. Validate inputs before submitting expensive jobs

Perform read-only checks:

1. Confirm every path exists and compressed VCFs have readable tabix indexes.
2. Confirm chromosome labels agree among the store, panel VCF, ancient VCF, and
   `CHROMS_FILE` (`1` versus `chr1` is a common silent failure).
3. Confirm REF/ALT and coordinates agree between both VCFs at a sample of sites.
4. Inspect store metadata and report its number of posterior draws. One draw is
   appropriate for pipeline validation. If the supplied store contains more,
   using all of them is scientifically acceptable, but it is no longer a
   one-draw pilot. Do not destructively edit the store to select a draw.
5. Confirm the polarity table was built from the same draw set and has matching
   row/draw dimensions.
6. Determine the maximum modern-panel haplotype count `N_SAMPLE`. Confirm the
   table will cover every called size from `MIN_N=20` through `N_SAMPLE`.
7. Confirm the Ne file has columns `time_left`, `time_right`, and
   `effective_population_size`, non-overlapping intervals, and coverage through
   all required ages.

Stop on coordinate, allele, draw-count, or sample-identity disagreement.

## 3. Choose the sample-age grid

Convert the oldest plausible archaeological/radiocarbon age to ARG generations
using the generation time agreed with the user. Let `TMAX` exceed that value by a
reasonable margin; the posterior must not pile up at its upper boundary. For the
pilot, use 300--500 grid points:

```bash
TMAX=30000       # replace after inspecting the dates and generation time
NT=400
N_SAMPLE=26      # replace with the modern ARG-panel haplotype count
MIN_N=20
TABLE="$RUN_ROOT/freq_table.npz"
```

Choose `AGEMAX` to extend beyond the generation age corresponding to
`mutation-age-max = 3` under the supplied demographic history. The table maximum
must be larger than the inference cutoff. Avoid the historical `4e7` default if
it is far beyond the required cutoff, because very old exact rows are expensive.
Document how `AGEMAX` was obtained rather than guessing silently.

## 4. Build the frequency table once

Submit the supplied wrapper after activating the environment in the submission
shell:

```bash
cd "$PROJECT"
NE="$NE_FILE" OUT="$TABLE" NSAMPLE="$N_SAMPLE" MIN_N="$MIN_N" \
TMAX="$TMAX" NT="$NT" AGEMAX="$AGEMAX" NAGE=100 \
PROJECT="$PROJECT" sbatch slurm/run_precompute.sbatch
```

The wrapper reserves **2 days** on one core. Cost is dominated by the oldest age
rows, not the grid size: at n=26 one (panel, 300 sample-age) row takes ~0.02 s at
tau_i = 5e-4 and ~2.9 s at tau_i = 2.3, but ~42 s at tau_i = 2000 where ~930
digits are needed. `AGEMAX` is therefore the knob that matters — see the note in
the wrapper header. If the job still times out, cut `AGEMAX` before raising the
reservation.

After completion, require a zero exit status and inspect the log. Then validate:

`table`, `table2` and `log_den` are NaN for `d0 > n` by construction, so a blanket
`isfinite().all()` fails whenever `MIN_N < N_SAMPLE`. Check the **eligible** cells
only:

```bash
TABLE="$TABLE" python - <<'PY'
import json, os, numpy as np
d = np.load(os.environ["TABLE"], allow_pickle=True)
meta = json.loads(str(d["meta"]))
n_panel = np.asarray(d["n_panel"])
print({k: d[k].shape for k in ("table", "table2", "log_den")})
print("n_panel", n_panel, "T", (d["Tgrid"][0], d["Tgrid"][-1]))
print("arithmetic", meta["arithmetic"], "unreachable", meta.get("unreachable_entries"))
bad = 0
for i, n in enumerate(n_panel):
    for k in ("table", "table2", "log_den"):
        cells = d[k][i, :n]                      # d0 = 1..n only
        n_nan = int(np.isnan(cells).sum())
        if n_nan:
            bad += n_nan
            print(f"  NaN in eligible cells: panel n={n} plane {k}: {n_nan}")
    # the d0 > n region must be NaN, not zero -- a wrong zero is invisible later
    pad = d["table"][i, n:]
    assert pad.size == 0 or np.isnan(pad).all(), f"panel n={n}: d0>n cells are not NaN"
assert bad == 0, f"{bad} NaN in eligible cells"
assert meta.get("unreachable_entries", 0) == 0, "precision ceiling was hit"
print("table OK")
PY
```

Never replace NaN with zero: a wrong zero satisfies every moment identity and so
becomes invisible downstream.

## 5. One-chromosome smoke test

Use the first chromosome and all 50 selected samples. Request per-sample curves
for this small pilot:

```bash
CHROM=$(head -n 1 "$CHROMS_FILE")

python "$PROJECT/posterior_sample_age_infer.py" \
  --freq-table "$TABLE" \
  --store "$STORE" --draw-polarity "$POLARITY" \
  --panel-vcf "$PANEL_VCF" --vcf "$ANCIENT_VCF" \
  --samples-file "$SAMPLES_FILE" --chrom "$CHROM" \
  --ploidy 1 --min-n 20 --mutation-age-max 3 \
  --marginalise uniform --epsilon 0.01 --per-sample-tsv \
  --output "$RUN_ROOT/results/smoke_$CHROM"
```

Before scaling up, inspect `run.json` and require:

- exactly 50 output samples in the requested order;
- nonzero `sites_used` and plausible retained-site counts;
- no unexpected mass rejection from missing draws, polarity, multiply mapped
  mutations, table coverage, or coordinate resolution;
- finite likelihoods and posterior summaries;
- no sample posterior concentrated at `T=0` or `TMAX` without investigation.

Compare posterior ages with radiocarbon dates only after converting both to the
same units and defining whether the radiocarbon value is a point estimate or a
calendar-age distribution.

## 6. Genome-wide chromosome array

The existing wrapper runs all selected ancient samples together for one
chromosome per array task:

```bash
NCHR=$(awk 'NF && $1 !~ /^#/' "$CHROMS_FILE" | wc -l)
TABLE="$TABLE" STORE="$STORE" POLARITY="$POLARITY" \
PANELVCF="$PANEL_VCF" VCF="$ANCIENT_VCF" \
SAMPLES_FILE="$SAMPLES_FILE" CHROMS_FILE="$CHROMS_FILE" \
EPSILON=0.01 PLOIDY=1 MIN_N=20 MUT_AGE_MAX=3 MARGINALISE=uniform \
OUTROOT="$RUN_ROOT/results/genome_parts" PROJECT="$PROJECT" \
sbatch --array="0-$((NCHR-1))" slurm/run_infer.sbatch
```

Every scientific setting is passed explicitly, so the array cannot silently
diverge from the smoke test if a Python default changes. The wrapper duplicates no
defaults of its own — each of `EPSILON`, `PLOIDY`, `MIN_N`, `MUT_AGE_MAX`,
`MARGINALISE` and `PER_SAMPLE_TSV` is pass-if-set. Still confirm them in each
part's `run.json`.

If cluster policy requires a different account, partition, memory or time,
override with `sbatch` flags rather than editing the shared wrapper.

Monitor failures and rerun only failed array indices. Do not merge partial output
— the merge step now refuses to run unless every chromosome part is present.

## 7. Merge chromosomes

After every chromosome succeeds:

```bash
TABLE="$TABLE" MERGE=1 CHROMS_FILE="$CHROMS_FILE" \
OUTROOT="$RUN_ROOT/results/genome_parts" \
PROJECT="$PROJECT" sbatch slurm/run_infer.sbatch
```

`CHROMS_FILE` must be passed here too: the merge derives its parts from that file
rather than globbing, and **exits 93 listing any missing chromosome** instead of
silently merging what happens to be on disk. (It previously globbed
`"$OUTROOT"/chr*`, which matched nothing whenever chromosomes are named `1..22`
rather than `chr1..chr22`.)

The wrapper writes the merge beneath
`$RUN_ROOT/results/genome_parts/genome/`. Check that `ages_table.tsv` has exactly
50 rows and that `ll_marginal.npy` has shape `(50, NT)`.

## 8. Required report to the user

Return:

- git commit and input paths/checksums;
- number of ARG draws actually used;
- table settings and runtime;
- per-chromosome and total retained/rejected site counts;
- peak memory and wall time for the smoke test and array jobs;
- genome-wide `ages_table.tsv` and posterior curves;
- inferred mean, MAP, and 95% interval beside each radiocarbon estimate;
- residual summaries and plots against radiocarbon age, missingness, coverage,
  damage/error proxies, and chromosome;
- counts of posteriors touching either T-grid boundary;
- explicit caveat that a one-draw run validates the pipeline but does not
  propagate ARG posterior uncertainty.

Do not tune epsilon, priors, filters, or the generation time to improve agreement
with the same 50 radiocarbon dates. First report the prespecified run. Sensitivity
analyses must be separately labelled, and final validation should use held-out
dated samples.
