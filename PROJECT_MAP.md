# Project map

This document describes the repository as it exists now. Paths containing large simulation or real-data artifacts are outputs, not source code; inspect their `run.json` and logs before deleting them.

## Top level

| Path | Classification | Purpose and retention |
|---|---|---|
| `precompute_freq_trajectory_moments.py` | production | Exact neutral Wright–Fisher moment-table generator (`freq_table.npz`); rerunnable from its demography/grid arguments. |
| `posterior_sample_age_infer.py` | production | Main per-chromosome and genome merge inference; consumes table, interval store, polarity, panel VCF and ancient VCF. |
| `precompute_betabinom_variable_ne.py`, `betabinom_real_data.py`, `prepare_betabinom_calls.py` | production/active | Variable-Ne beta-binomial table, real-data inference, and call/site preparation. Preserve scripts and final summaries; intermediate arrays can be regenerated. |
| `direct_frequency_age_infer.py`, `diffusion_frequency_validation.py`, `msprime_exact_time_validation.py`, `validate_moments_vs_mc.py` | validation | Maintained simulation-calibration harnesses; see the code map below. Keep code and compact summaries; bulk replicate data are reproducible from recorded seeds. |
| `slim/` | validation simulation input | `single_site_neutral.slim`, the forward Wright–Fisher model used by the single-site tests. |
| `slim_single_site_direct_age/` | generated, compact validation result | Known-frequency likelihood control. Raw simulations/chunks were intentionally removed after documentation; keep the tracked summary, PNGs and `run.json`. |
| `slim_single_site_diffusion_test_A/` | generated validation result | First, over-conservative unique-lineage diffusion check. Raw chunks were removed; compact results and its constant-Ne table remain. |
| `slim_single_site_diffusion_test_A_one_per_replicate/` | generated validation result | Corrected SLiM diffusion check with at most one panel-polymorphic lineage per replicate. Raw chunks were removed; keep compact results. |
| `slim_edge_interval_validation.py`, `edge_uniformity_validation.py`, `slim_single_site_edge_validation/` | validation, completed | T9 independent-locus control. On the same 1,243 loci, exact mutation time gives +33-generation MAP bias with all seven 95% intervals covering; the true edge plus uniform marginalisation gives +376. All mutation times lie inside their edges and the discrete within-edge placement is uniform (KS p=0.346). Preserve compact results, seeds and scripts. |
| `msprime_exact_time_ne50k/` | generated completed experiment | Simulations, exact-time and edge-interval inference, tables and results for Ne=50,000. Retain final results and `run.json`; raw replicates are reproducible. |
| `logan_try/` | active real-data experiment | Cohort selection, VCF audits, epsilon calibration, and chromosome/cohort results. Treat VCF/TSV/NPZ inputs and completed results as valuable provenance; do not delete without an external copy. |
| `results/` | generated output | Main posterior/ESS figures and run JSONs. Safe to regenerate only when inputs and settings are retained. |
| `scripts/` | analysis and provenance code | Calibration checks, draw-mixture ESS, edge-width scaling, T1/T4 attribution analyses, and `scripts/provenance/` historical big-simulation recipes (some contain machine-specific absolute paths). |
| `slurm/` | execution wrappers | sbatch wrappers for every precompute, simulation, inference, merge, and preparation stage. |
| `tests/` | test suite | Unit/regression tests for moment engines, likelihood marginalisation, orientation/chunking, numerical guards, adapters, and pipeline integration. |
| `normalizeTEs/` | pinned git submodule dependency | Supplies the `normalize_tes` interval-store, polarity and VCF adapters; inference imports it via `PYTHONPATH`. Initialize with `git submodule update --init`; **currently uninitialised (the directory is empty)**, which is why `tests/conftest.py` stubs the adapter boundary. Do not treat it as generated output. |
| `working_diffusion.md`, `working_betabinom.md`, `MATH.md`, `NOTES.md` | design/provenance docs | Derivation, modelling decisions and current investigations. Read alongside [README.md](README.md). |
| `README.md`, `HPC_REAL_DATA_CODEX.md`, `ANCIENT_TEST.md`, `TODO.md` | user/run documentation | Pipeline usage, cluster notes, ancient-data checks and open work. |
| `logs/`, `slurm-*.out`, `__pycache__/`, `.pytest_cache/` | generated transient output | Job logs/caches; safe to delete after diagnosing runs (retain failed logs when provenance matters). |

## Code map

| File | Role |
|---|---|
| `precompute_freq_trajectory_moments.py` | Builds exact-arithmetic first- and second-moment diffusion tables under piecewise-constant Ne. |
| `posterior_sample_age_infer.py` | Production diffusion inference, including VCF/store adapters, edge-age lookup, ARG-draw marginalisation, epsilon handling and chromosome merging. |
| `prepare_betabinom_calls.py` | Converts the panel/ancient calls and ARG information into the site archive consumed by the beta-binomial branch. |
| `precompute_betabinom_variable_ne.py` | Builds variable-Ne beta-binomial probability tables. |
| `betabinom_real_data.py` | Runs beta-binomial inference from prepared calls and tree sequences. |
| `direct_frequency_age_infer.py` | Validates the terminal Bernoulli likelihood using exact SLiM frequency trajectories—no frequency approximation. |
| `diffusion_frequency_validation.py` | SLiM Test A: exact focal-mutation origin plus sampled modern count passed through the diffusion lookup. |
| `msprime_exact_time_validation.py` | Generates the 100 Ne=50K msprime replicates and runs the paired exact-mutation-time and true-edge-interval inference tests. |
| `slim_edge_interval_validation.py` | T9 inference/merge harness for 10,000 independent one-base SLiM replicates, supporting paired exact-time and true-edge analyses. |
| `edge_uniformity_validation.py` | T9 audit of mutation containment, discrete uniform placement within the true edge, and recurrent-replacement count mismatches. |
| `scripts/edge_width_scaling.py` | T8: paired exact-time vs true-edge inference on identical sites within edge-width strata, which separates interval *width* from interval *representation*. |
| `scripts/draw_mixture_ess.py` | T7: reconstructs per-draw log-likelihoods from saved epsilon data to measure draw-mixture ESS and the cost of per-site draw averaging. |
| `validate_moments_vs_mc.py` | Compares moment recursion against a Monte Carlo simulation of the diffusion itself; numerical validation, not a discrete-population test. |

## Test map

`tests/conftest.py` supplies fake stores, VCF chunks, arguments and tables shared by
the production-inference tests. `tests/_reference.py` is the independent
high-precision mathematical reference.

| Test file | Protects |
|---|---|
| `test_moment_engine.py` | Neutral Wright–Fisher moment equations, first/second moments, boundaries and Monte Carlo agreement. |
| `test_exact_moment_engine.py` | High-precision engine versus the independent reference, including float64 cancellation failures and adaptive precision. |
| `test_numerical_guards.py` | Failure-on-invalid numerics, CLI defaults, priors, panel sizes and merge validation. |
| `test_lookup_existence_boundary.py` | Exact `p_T=0` mutation-existence boundary and branch interpolation/integration behavior. |
| `test_branch_marginalisation_commutes.py` | That within-edge mutation-age marginalisation at the frequency level equals it at the likelihood level, which holds only while the site likelihood stays affine in the tabulated moments. Also checks the knot-split integrator against an independent quadrature. Research regression test, not basic I/O. |
| `test_chunk_orientation.py` | Record-by-sample VCF chunk orientation and allele-count correctness. |
| `test_normalizetes_adapter.py` | Compatibility with the current `normalizeTEs` chunk API and filtering of adapter events. |
| `test_draw_marginalization.py` | Correct order: accumulate sites within each ARG draw, then marginalise draws. |
| `test_draw_retention.py` | Complete-draw requirements, polarity availability and rejection accounting. |
| `test_draw_subset.py` | Explicit ARG-draw selection and validation. |
| `test_epsilon_calibration_output.py` | Preservation of site-level sufficient data used to profile epsilon. |
| `test_betabinom_real_data.py` | Variable-Ne beta-binomial construction, orientation and likelihood behavior. |
| `test_direct_frequency_age_infer.py` | Known-frequency parsing, Bernoulli accumulation and posterior normalization. |
| `test_diffusion_frequency_validation.py` | Recurrent-lineage parsing and reproducible one-lineage modern-panel selection for SLiM Test A. |
| `test_msprime_exact_time_validation.py` | Reproducible seed streams and 26-modern-plus-one-ancient VCF parsing for the paired msprime test. |

Run the default suite through SLURM rather than on the head node, e.g.
`HPC_MEM=16G ~/.claude/bin/hpc_run 'python -m pytest tests/ -q'`; the wrapper activates
`environment.yml` itself. Note that `pytest.ini` declares a `slow` marker but sets no
`addopts`, so slow tests are *not* actually deselected by default despite the comment there.

## SLURM wrapper map

| Stage | Wrappers |
|---|---|
| Production diffusion | `run_precompute.sbatch` → `run_infer.sbatch` (submit per chromosome) → merge mode through `run_infer.sbatch`. |
| Real-data beta-binomial | `prepare_betabinom_calls.sbatch` plus `precompute_betabinom_variable_ne.sbatch` → `run_betabinom_real.sbatch` → `merge_betabinom_real.sbatch`. |
| Original SLiM simulation | `run_single_site_slim.sbatch` → `merge_single_site_slim.sbatch`. The 48 GB raw output has been removed. |
| Known-frequency SLiM control | `run_direct_frequency_age_chunks.sbatch` → `merge_direct_frequency_age.sbatch`. |
| Diffusion SLiM Test A | `precompute_constant_ne_frequency_table.sbatch` → `run_diffusion_frequency_validation.sbatch` → `merge_diffusion_frequency_validation.sbatch`. |
| msprime exact-time control | `run_msprime_exact_time_simulations.sbatch` plus `precompute_msprime_ne50k_exact_time.sbatch` → `run_msprime_exact_time_inference.sbatch` → `merge_msprime_exact_time_inference.sbatch`. |
| Paired true-edge control | Reuses the msprime simulations/table, then `run_msprime_edge_interval_inference.sbatch` → `merge_msprime_edge_interval_inference.sbatch`. |
| Edge-width scaling (T8) | Reuses the msprime simulations/table, then `run_edge_width_scaling.sbatch` → `merge_edge_width_scaling.sbatch`. |
| SLiM independent-edge control (T9) | `run_slim_edge_validation_simulations.sbatch` → optional retry consolidation → `run_slim_edge_validation_inference.sbatch` → `merge_slim_edge_validation.sbatch`; `check_slim_edge_uniformity.sbatch` runs the placement audit. Completed exact and edge arms are retained separately. |

## Generator → output relationships

* `precompute_freq_trajectory_moments.py` (or `slurm/run_precompute.sbatch`) → frequency table (`.npz` plus grid/metadata), consumed by `posterior_sample_age_infer.py` and validation scripts.
* `posterior_sample_age_infer.py` via `slurm/run_infer.sbatch` → one output per chromosome (`ages_table.tsv`, likelihood/grid/sample arrays, `run.json`), then a genome-wide merge directory.
* `slim/single_site_neutral.slim` via `slurm/run_single_site_slim.sbatch` → replicate `frequencies/` and `samples/`; `merge_single_site_slim.sbatch` combines them. Direct/diffusion harnesses analyse those files; their merge wrappers produce calibration plots and JSON summaries. These raw files were deleted after the documented controls passed and must be regenerated before rerunning the SLiM analyses.
* `slim/single_site_edge_validation.slim` via `run_slim_edge_validation_simulations.sbatch` → 10,000 compact modern-only tree sequences plus historical genotype manifests. `slim_edge_interval_validation.py` produces `results_exact_time/` (+33-generation MAP bias) and `results_random_focal/` (+376); `edge_uniformity_validation.py` produces `edge_uniformity/` (all 1,243 mutations contained, uniform-placement KS p=0.346). The independent-locus result rules out LD as the edge-bias explanation.
* `msprime_exact_time_validation.py`: simulation → inference → merge, using `msprime_exact_time_ne50k/{simulations,exact_time_inference,edge_interval_uniform_inference}`. Exact-time and edge-interval branches are completed; edge result is `edge_interval_uniform_results/run.json` (MAP bias +363.56, coverage 0.44), while exact-time results are under `msprime_exact_time_ne50k/results/` — not the top-level `results/`, which is a different directory. Edge-width strata are in `edge_width_scaling_results/` (T8).
* `prepare_betabinom_calls.py` → prepared cohort/site archive; variable-Ne precompute → beta-binomial table; `betabinom_real_data.py` jobs → merged real-data posterior. Wrappers are the corresponding files in `slurm/`.
* `scripts/check_frequency_calibration.py`, `draw_mixture_ess.py`, `compare_t1_marginalisation.py`, and `t4_attribution_sweep.py` read existing outputs and write diagnostics; they do not define production inference.

## SLURM ordering

Main pipeline: **precompute → infer array (one chromosome per task) → merge**. SLiM: **simulate array → merge replicates → direct/diffusion analysis array → merge analysis**. msprime: **simulate array → precompute table → infer array → merge**, repeated for exact-time and edge-interval branches. Real beta-binomial: **prepare calls and precompute variable-Ne table → inference jobs → merge**. Verify output completeness before submitting a merge.

## Retention and reproducibility

Retain source, tests, documentation, seeds, `run.json`, final tables and figures, and real-data inputs. Raw simulation chunks and job logs are removable only after confirming seeds, parameters, and compact summaries remain. Generated directories are not source. In particular, preserve `logan_try/` unless its real-data inputs and outputs have an external archive. Retain `slim_single_site_edge_validation/results_exact_time/`, `results_random_focal/`, and `edge_uniformity/`; the 10,000 compact trees are reproducible from the fixed seeds and may be removed only after the completed T9 result is committed and archived. The 2.6 GB `msprime_exact_time_ne50k/simulations/` directory can eventually be regenerated from its per-replicate metadata, but retain it until all paired edge analyses are documented. Reproduction requires [environment.yml](environment.yml), the pinned `normalizeTEs` submodule, and external VCF/ARG/demography files listed in [README.md](README.md). Mathematical definitions are in [MATH.md](MATH.md), with implementation notes in [working_diffusion.md](working_diffusion.md) and [working_betabinom.md](working_betabinom.md).
