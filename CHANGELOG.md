# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/) in its `0.x` form: a minor bump marks a new approach or milestone, a
patch bump a fix. Nothing here is validated on real data yet, which is why the major version is 0.

Versions 0.1.0 to 0.4.1 were tagged retroactively on 2026-11-13, on the commits that ended each milestone.

## [Unreleased]

## [0.4.1] - 2026-11-13

### Fixed
- Non-sharing ratio: the coalescence rate c(t) taken from an ARG now uses local-tree lineage-pair moments
  (`lineage_moments`) instead of lineage counts alone. In practice this rescales T_hat. On the held-out chr10 real-data
  check (10 draws, child-node A0 20,000) the median T_hat goes from about 10,200 to 15,700 generations, while f and the
  correlation with calibrated age stay the same.
- Panel-fixed sites keep the actual number of called panel haplotypes as their reference weight, so uncalled panel
  samples no longer contribute (site, reference) pairs.
- `merge` rejects draws whose A0, age mode, c(t) source, chromosomes, position shift, panel-fixed setting, panel size,
  hazard grid or count-array shape differ from the first draw's.

## [0.4.0] - 2026-11-13

### Added
- Non-sharing-ratio estimator: age from the excess of non-shared derived alleles in an ancient genome relative to
  held-out modern controls, inverted through the integral of c(t) (`pipeline/nonsharing_ratio*.py`).
- c(t) from the ARG's own lineage counts, from an Ne trajectory, or from a demes model.
- Real-data CLI (`run-draw`, `merge`, `prepare-fixed`), SLURM wrappers in `slurm/`, and the guide
  `docs/nonsharing_ratio/RUNNING.md`.
- Simulation validation with SINGER + POLEGON ARGs (r 0.90 to 0.95 with true age), including bottleneck and
  split-population structure scenarios.
- Real-data offset report (`docs/nonsharing_ratio/REAL_DATA_OFFSET_REPORT.md`) based on a held-out-controls stand-in
  for chromosome 10.

### Changed
- The front-page README now describes the non-sharing ratio. Earlier approaches are summarised in
  `docs/PREVIOUS_APPROACHES.md`, and the old front page is kept as `docs/PRESENT_COUNT_DIFFUSION_README.md`.
- Open work moved to GitHub issues labelled `non-sharing-ratio`. Issues #1 to #5, filed against earlier approaches,
  were closed as not planned.

## [0.3.0] - 2026-11-09

### Added
- POLEGON insertion calibration: diagnostics on POLEGON-dated true-topology ARGs matched to the production setup,
  per-chromosome calibration, controlled recovery tests and SLURM workflows.
- Methods and results guide `POLEGON_INSERTION_VALIDATION.md`.

## [0.2.0] - 2026-10-28

### Added
- Radiocarbon epsilon calibration for the insertion model (`pipeline/radiocarbon_calibration.py`): profiles the
  effective discordance epsilon against radiocarbon densities, with cross-validation and a bootstrap over dated
  individuals.
- `--epsilon-grid` for `insertion_real_data.py`, which gives a joint age and epsilon likelihood per ARG draw, with
  checkpoints for long runs.
- Joint calibration of ARG time scale and epsilon.

## [0.1.0] - 2026-10-12

The state of `main` before the move to the non-sharing ratio (branch `insertion` points at the same commit).

### Added
- Present-count diffusion approach: an age posterior from derived-allele frequency, conditioned on mutation age and
  present panel count, under a Wright-Fisher moment recursion. The frequency table is computed in exact arithmetic.
- Ancient-lineage insertion approach: the ancient lineage is inserted into the modern ARG at age T. Validated on
  linked msprime ARGs, across varied simulations and under variable demography.
- Mathematical derivation in `MATH.md` and working records (`working_diffusion.md`, `working_insertion.md`).
- Graphical abstract for the README and slides.

### Changed
- Scripts reorganised into `pipeline/`, `validation/` and `scripts/`.
- epsilon is now interpreted as an effective discordance rather than a genotyping error rate.

### Removed
- The tree-at-T beta-binomial approach was archived to branch `betabinom` (tag `betabinom-archive-2026-09-12`).

[Unreleased]: https://github.com/rossibarra/dnAging/compare/v0.4.1...HEAD
[0.4.1]: https://github.com/rossibarra/dnAging/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/rossibarra/dnAging/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/rossibarra/dnAging/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/rossibarra/dnAging/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/rossibarra/dnAging/releases/tag/v0.1.0
