# Non-sharing-ratio estimator on real data: findings and status (2026-11-11, overnight session)

Everything here is aggregate counts only. Ancient samples appear as S01.. (VCF column order).
No ancient VCF or radiocarbon table was copied into the repo. Probe scripts: `probe_arg.py`,
`probe_fixed.py`; their logs are in `logs/`. How to run: `RUNNING.md`.

**Update 2026-11-12:** the first real-data check against radiocarbon ages found no signal. Ancient samples show a
constant ~31% excess of non-sharing (H = 0.27) that is not explained by read depth or simple reference bias;
population structure is the untested leading hypothesis. See [`REAL_DATA_OFFSET_REPORT.md`](REAL_DATA_OFFSET_REPORT.md).

## Start here: map of the non-sharing-ratio work (as of 2026-11-12)

Method: the report "Dating ancient genomes with inferred ARGs" (Yun Deng, 2026-11-11). Open items: `TODO_nonsharing_ratio.md`
(repo root). Estimator code: `pipeline/nonsharing_ratio*.py`; tests: `tests/test_nonsharing_ratio*.py` (run via `hpc_run`).

| What | Where | One-line result |
|---|---|---|
| Real-data inputs, conventions, fixed sites, plumbing | this file, sections 1-7 below (written 2026-11-11) | wiring works on chr1; position shift 0; ancestral/derived from the ARG |
| Held-out stand-in on real data (chr10, 6 controls, 20-haplotype panel, 10 POLEGON draws) | `batches/nonsharing_ratio_heldout_20261111/` | the test bed for everything on real data |
| Real-data age recovery and the constant offset | `REAL_DATA_OFFSET_REPORT.md` (sections 1-5) | r about 0 with radiocarbon age; ancient f is 1.31x controls for every sample |
| Depth and reference-bias tests | `REAL_DATA_OFFSET_REPORT.md` sections 3, 4 | neither explains the offset |
| Simulations, no structure (constant size; bottleneck) | `batches/nonsharing_ratio_realsims_20261112/`, `batches/nonsharing_ratio_realsims_bottleneck_20261112/`; summary in report section 6 | method recovers ages at A0 10k and 20k (r 0.90-0.94) |
| More draws / more chromosomes | `batches/nonsharing_ratio_heldout_20261111/{more_draws_check,block_check}/`; report section 7 | neither helps |
| Structure test (ancient from a population split 6,000 generations ago, same bottleneck) | `batches/nonsharing_ratio_realsims_split6k_bottleneck_20261112/` | launched 2026-11-12; read its README for the result |
| Graphical abstract | `graphical_abstract.{png,svg}` (built by `scripts/build_nonsharing_ratio_abstract.py`) | schematic only |

Conventions worth remembering: time is in generations (1 generation = 1 year assumed for real data); A0 is applied to the
child-node time of the mutation edge (`--age-mode child`; the default edge-midpoint mode keeps 88-92% of sites and is
not what was used for the later results); radiocarbon dates and sample metadata live only in the gitignored `private/` directory (`private/cmp.py`, `private/depth.py`,
`private/dates.txt`, `private/ages.txt`, the metadata table `private/meta_v10_9July2026.csv` converted from the .xlsx by `private/xlsx_to_csv.py`, and `private/check_vcf_order.py`). `radiocarbon_goodvcfs/` (ancient VCFs and PLINK sets) is gitignored too. The ancient VCF sample names are the metadata `FAM_Alias` values, and `check_vcf_order.py` confirmed on 2026-11-12 that all 50 trusted50 columns match the dates in `dates.txt` by name (0 mismatches) and are never committed;
`hpc_run` jobs cannot see the head node's scratchpad, so scripts and outputs live in the repo; use bcftools and plink by
absolute path inside `hpc_run`; `python -I` when a stray `inspect.py` could shadow the stdlib.

## 1. Where the real inputs are (paths only)

| Input | Path | Notes (verified by reading headers / one draw) |
|---|---|---|
| ARG draws (raw SINGER, whole genome concatenated) | `/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/results/combined/tsz_files/run.combined.<N>.tsz` | 75 `.tsz` files (the README there says they are the 75 argtest draws; file 100 is ~1.1 GB; `burin_tsz/` holds burn-in). 26 haploid samples, 25,030,335 sites, sequence length 2,131,846,815. No `chrom_offsets` metadata (keys: kept_intervals, mu_position, mu_rate), so `--chrom-fai` is required. |
| Chromosome offsets | `/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/ref.fa.fai` (also `/quobyte/jrigrp/beil/logan_collab/chrom_offset.txt`) | `chromosome_layout` reconstructs offsets with +1 padding per chromosome; reproduces the ARG length exactly after chr10. Keys are `1`..`10`. |
| Modern panel VCF (whole genome) | `/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/results/combined/nam_te_vcf/combined.vcf` | haploid, 26 columns (first column is named REF), variable sites only (`SC=variant`), ~8.7% missing calls at ARG sites. |
| Modern panel VCF, chr1 only | `radiocarbon_goodvcfs/filtered_trusted_20261012/panel_chr1.vcf.gz` | made by `extract_panel_chr1.sh` from the above |
| Ancient pseudo-haploid VCF (chr1 only, numeric chr label) | `radiocarbon_goodvcfs/filtered_trusted_20261012/chr1_numeric.trusted50.mappable_q37.vcf.gz` (50 samples); `chr1_numeric.all55.mappable_q37.vcf.gz` (55) | Whole-genome versions exist only as PLINK files `ARG_site_ancient_haps_mind95.{trusted50,all55}.mappable_q37.{bed,bim,fam}` in the same directory; a genome-wide numeric VCF has to be exported (see `export_chr1_vcf.sh`, which also relabels `NC_050096.1` -> `1`). The raw `chr1_accession.*.vcf` uses accession labels, which the CLI does not map. |
| Ne trajectory | `radiocarbon_goodvcfs/insertion_pilot_20261012/demography/posterior_mean.from_zero.tsv` (series `posterior_mean`); raw per-draw estimates `.../argtest/results/combined/ne_plot/coalescence-ne-estimates.tsv` | for `--ne` |
| Radiocarbon info | outside the repo (see `real_data_plan.md`); not needed by this estimator except to choose A0 | I did not read it |
| Per-draw polarity / interval store | `logan_try/polarity_draw100/`, `logan_try/store_draw100/` (draw 100 only) | per-ARG-site tables; do not cover non-ARG sites |

POLEGON caveat: the real draws above are the raw SINGER `.tsz` (`run.combined.*`). I could not find
POLEGON-re-dated real-data draws in any run record: `batches/singer_polegon_pilot_20261102/` and
`batches/*/` POLEGON runs are simulations (their trees point at `batches/nine_demographies_20260916/runs/...`).
The CLI takes any tszip/tskit draw, so POLEGON-dated draws can be dropped in unchanged (node times are the only thing used),
but who produced real POLEGON draws and where they live is a question for you.

## 2. Coordinate convention (real data)

Command: `docs/nonsharing_ratio/probe_arg.py` on `run.combined.100.tsz` (log `logs/probe_polarise_draw100.log`), chr1 only,
`nonsharing_ratio_io.polarise` with `position_shift` in {-1, 0, +1}:

| VCF | shift | ARG-site matches | allele mismatches |
|---|---|---|---|
| ancient (975,095 eligible SNP records) | -1 | 20,405 (2.1%) | 59,590 (74% of position hits) |
| ancient | **0** | **697,329 (71.5%)** | 33,054 (4.5% of position hits) |
| ancient | +1 | 20,331 (2.1%) | 59,738 (75%) |
| panel as controls (7,347,783 records) | -1 | 178,381 | 553,793 |
| panel | **0** | **4,398,087 of 4,458,589 chr1 ARG sites** | **0** |
| panel | +1 | 178,332 | 553,793 |

Shift 0 is correct (the panel VCF is the ARG's source, so 0 mismatches is the expected outcome; the ancient VCF matches ten times
more records at 0 than at +-1). CLI default `--position-shift 0` is unchanged. The 4.5% allele mismatch at shift 0 in the ancient VCF is
sites where the ancient VCF's REF/ALT pair differs from the ARG's ancestral/derived pair (provisional REF in the PLINK export, other allele, or strand);
`polarise` drops them (counted in `meta.ancient_polarise.allele_mismatch`). Cause not investigated (unverified).
Chromosome names: ARG keys from the FAI are `1`..`10`; the numeric ancient and panel VCFs use `1`. The CLI maps `1` <-> `chr1` both
ways (unit test `test_chromosome_name_styles_are_interchangeable`); `NC_050096.1` accession labels are not mapped (use the numeric export).

## 3. Panel-fixed sites (design rule 2)

What the real data show (`logs/probe_fixed_draw100.log`, `logs/probe_fixed2_draw100.log`, chr1, draw 100):

* The ARG already contains panel-fixed-derived sites: 103,335 of 4,458,589 chr1 ARG sites (2.3%) have all 26 panel haplotypes derived
  (`edge_midpoint_ages` gives them +inf; `select_sites` includes them, as rule 2 requires). Genome-wide 1,078,466 of 25,030,335 sites have age +inf.
  In the panel VCF those sites are near-fixed, not exactly fixed: median derived fraction among called haplotypes 0.92 (5-95%: 0.78-0.96)
  versus 0.08 for other sites. So they are "all derived in the ARG's genotypes", not VCF-fixed.
* The panel VCF lists variable sites only. `scan_panel_fixed` on the whole chr1 panel VCF (7,347,783 biallelic SNP records) finds **0** sites
  fixed among called haplotypes. Sites that are monomorphic across the panel are absent from every input I found (no invariant-site VCF; the header
  mentions `SC=invariant` but no such records exist). 2.9M panel-VCF records on chr1 are not ARG sites (filtered before SINGER); they are
  variable, so they are not "panel-fixed" either.
* Polarity for any such site would have to come from outside the panel. The ARG's ancestral state equals the panel VCF REF at only 69.4% of
  chr1 ARG sites (30.6% = ALT), so "REF = ancestral" is not a valid fallback here. There is no outgroup/ancestral FASTA in the run records I read.

Implemented (optional, off unless `--panel-fixed` is given): `run-draw`'s `prepare-fixed` step streams a panel VCF (that includes monomorphic sites) and keeps
sites fixed among >= `--min-called` called haplotypes; `--panel-fixed` + `--polarity-table chrom pos ancestral` (default polarity mode `table`; `ref` exists
but is opt-in and invalid for this data) adds the sites whose panel allele is derived, with age +inf, skipping positions already ARG sites.
Unit tests cover the scanner, polarity modes, merge/ordering, and an end-to-end CLI test in which fixed sites deleted from the ARG are recovered and
the counts equal the full-ARG counts exactly.

Simulation check (`validation/nonsharing_ratio_fixedsites_validation.py`; 20 replicates, 20 Mb, Ne=20,000, 26 panel, 8 controls, edge-midpoint ages, A0=12,000,
true c(t); table `results/nonsharing_ratio/fixedsites_const20mb.tsv`, summary `results/nonsharing_ratio/fixedsites_const20mb_summary.txt`).
Mean bias in generations (SE about 80-120):

| true age | no_fixed | with_fixed | wrong polarity on 10% of candidate sites |
|---|---|---|---|
| 1,000 | -36 | +80 | +176 |
| 2,000 | -212 | +10 | +213 |
| 4,000 | -559 | +3 | +365 |
| 6,000 | -984 | -90 | +389 |
| 8,000 | -1,228 | +117 | +728 |

So excluding panel-fixed sites biases old samples young (as the report says) and including them removes it (within SE) in simulation.
A wrong polarity biases the other way (old). Limitation on real data: the ARG's own near-fixed derived sites are included, but panel-monomorphic sites outside the ARG
cannot be recovered without an invariant-site panel table and an ancestral-allele source; the real-data estimates are therefore subject to the "young" direction of that
bias to an unknown extent (the ARG already carries about 2.3% panel-fixed sites on chr1, which may recover part of the effect; unmeasured).

## 4. Other additions (all in new/own files, tests in `tests/test_nonsharing_ratio_fixedsites.py`)

* Flags: `target_flags` labels `f0`, `H<=0`, `no_data`; shown in `run-draw` stdout (`FLAG=...`), stored in the NPZ (`flags`), and a `flag` column in the merged TSV (per-draw counts).
* `--ne FILE.tsv [--ne-series NAME]`: c = 1/(ploidy*N) averaged per bin from the repo's Ne TSV (via `PiecewiseConstantNe.from_tsv`, last epoch extended); demes YAML still works; giving both is an error.
* Pseudo-haploid ancient VCF: the real file codes calls as diploid `0/0`, `1/1`, `./.` (verified on the first 2,000 records: only those three genotype strings, no hets; header lists a single `GT` FORMAT, `PR` INFO flag "provisional reference").
  `--ancient-mode sample` (default) reads `0/0`->0, `1/1`->1, `./.`->missing, haploid `0`/`1`/`.` too, and any het or half-missing -> missing. Records with ALT `.` (monomorphic in the panel) are skipped as non-SNP. Test: `test_pseudohaploid_gt_conventions`.
* `--chroms` (smoke tests), `--deidentify` (S01../C01.. labels), vectorised `panel_derived_matrix`, `edge_midpoint_ages`, `arg_site_alleles` (reference versions kept as `*_slow`, equality tested). On the real draw the vectorised versions take about 200 s and 430 s genome-wide (the slow loops were not timed to completion; I stopped one after > 10 min).
* Wrappers: `slurm/nonsharing_ratio_{prepare_fixed,run_draws,merge}.sbatch`; lint clean (0 errors, 0 warnings; `--test-only` accepted). Partition `high`, account `jrigrp` (the existing real-data wrappers use `high`; `low` would requeue and `run-draw` refuses to overwrite an existing output).

## 5. Smoke test on real data (plumbing only)

Command (job 40590270, 11 min 11 s wall, 1 CPU, MaxRSS 21.4 GB; log `logs/smoke_chr1_draw100.log`):
```
python -u pipeline/nonsharing_ratio_real_data.py run-draw --tree .../tsz_files/run.combined.100.tsz --draw-id 100 \
  --controls-vcf radiocarbon_goodvcfs/filtered_trusted_20261012/panel_chr1.vcf.gz --controls-mode sample \
  --ancient-vcf radiocarbon_goodvcfs/filtered_trusted_20261012/chr1_numeric.trusted50.mappable_q37.vcf.gz \
  --chrom-fai .../argtest/ref.fa.fai --chroms 1 --a0 20000 \
  --ne radiocarbon_goodvcfs/insertion_pilot_20261012/demography/posterior_mean.from_zero.tsv --ne-series posterior_mean \
  --deidentify --output docs/nonsharing_ratio/smoke/draw_100_chr1_PANELASCONTROLS.npz
```
**The controls here are the 26 ARG panel haplotypes. This violates design rule 1 and the results are biased; they are not age estimates.**
Observed: 50 ancient targets, 4,398,087 usable chr1 sites, 3,955,395 selected, no flagged targets, T_hat 3,646-6,873 (median 5,455) generations
(meaningless under this setup; A0=20,000 is a placeholder, which selected 90% of usable sites; my guess is that most ARG site ages in this draw exceed 20,000 generations, not checked, so A0 needs a real choice), merge (200 bootstraps) wrote
`smoke/merged_PANELASCONTROLS.tsv`. Phase times from the log: tree load + both VCF reads 465 s, ARG tables + panel matrix 111 s, polarisation 69 s, counts 22 s.

Cost, extrapolated and unverified: the whole-genome tree alone loads in ~60-180 s and 18.3 GB. VCF reading is pure Python and re-done per draw
(ancient chr1 83 s, panel chr1 264 s); scaling by genome size (about 7x chr1 sites) suggests ~1-1.5 h of reading per draw and tens of GB for the genotype lists,
so I set 100 GB / 6 h per draw task. First real full-genome run should be one draw and checked with `sacct -j ID --format=Elapsed,MaxRSS` before launching the array.
Caching parsed VCFs across draws would be the obvious optimisation (not done).

## 6. Open issues / decisions for you

1. Held-out controls VCF: not yet supplied; all real results require it (and its sample names must not overlap the panel or ancient VCF; the CLI checks control vs ancient overlap, and vs the panel if `--panel-names` is given).
2. Choice of A0 (generations, above the oldest sample age, in the ARG's time units). The ARG time scale for the raw SINGER draws vs POLEGON re-dating (`real_data_plan.md` mentions a time-scale/epsilon calibration) is unresolved here.
3. Whole-genome ancient VCF in numeric labels does not exist as a VCF (only chr1); export from the PLINK files first.
4. Panel-monomorphic (truly fixed) sites outside the ARG: need an invariant-site panel table plus an ancestral-allele source (outgroup). Until then rule 2 is only met for the ARG's own all-derived sites.
5. 4.5% allele mismatches between the ancient VCF and ARG alleles at position hits: cause unverified (strand, provisional REF, or other).
6. `PROJECT_MAP.md` / `TODO_nonsharing_ratio.md` not edited (rule: do not edit other existing files); suggested entries: the three `slurm/nonsharing_ratio_*.sbatch`, `docs/nonsharing_ratio/`.

## 7. Wrapper test

`slurm/nonsharing_ratio_run_draws.sbatch` was run for real (job 40590435, `CHROMS=1 DEIDENTIFY=1`, same panel-as-controls plumbing setup as section 5, `smoke_dir/` holds a symlink to draw 100):
COMPLETED in 6 min 15 s (faster than section 5, probably warm file cache, unverified), MaxRSS 21.5 GB, output `smoke_sbatch/draws/draw_100.npz`, identical selection counts (3,955,395 of 4,398,087 sites).
The prepare-fixed and merge wrappers were only lint-checked (`sbatch --test-only` accepted); their Python entry points were run (merge on the smoke draw; prepare-fixed scanner on the real chr1 panel VCF via the probe).
