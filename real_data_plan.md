# Real-data calibration and validation plan

## Objective

Use ancient individuals with high-confidence radiocarbon dates to:

1. estimate the effective discordance parameter, epsilon, for the insertion
   likelihood on real data; and
2. measure out-of-sample genetic age-estimation accuracy.

Calibration and validation must be separated. A sample used to estimate epsilon
cannot also provide an unbiased assessment of dating accuracy. The primary
analysis will therefore use cross-validation, with radiocarbon information used
only for the training samples in each fold.

This repository must not store sample-level radiocarbon metadata or ancient VCFs.
The analysis should read those inputs from their protected external locations and
write only aggregate summaries or explicitly de-identified validation results.

## Interpretation of epsilon

For estimated ARGs, epsilon is an **effective per-allele discordance parameter**,
not a pure sequencing-error rate. It may absorb:

- ancient-genotype error;
- residual allele-orientation or polarity error;
- SINGER topology error;
- POLEGON node-time error;
- Ne misspecification; and
- other model mismatch that makes an observed ancient allele incompatible with
  the inferred ARG.

The fitted value must be reported as effective discordance. It must not be
labelled or interpreted as a biochemical damage or genotype-call error rate.

The earlier real-data estimate near 6% was obtained with the diffusion pipeline
and ages fixed near radiocarbon values. It is useful context, but it is not an
insertion estimate and should not be transferred without recalibration.

## Inputs

Required external inputs are:

- the high-confidence radiocarbon sample set;
- sample identifiers matching the ancient VCF columns;
- calibrated radiocarbon-age information;
- the ancient pseudo-haploid VCF;
- the modern-panel VCF used for ARG inference;
- POLEGON-dated SINGER posterior ARG draws for each chromosome;
- the estimated Ne trajectory; and
- the existing insertion call preparation and common-draw masks.

For radiocarbon ages, the preferred input is the full calibrated probability
density for each sample. If only a `cal_BP_95` interval is available, the first
implementation may use a uniform density over that interval. A midpoint should
not be treated as an exact age.

All ages will use the same time direction and units as insertion inference:
generations before present. The generation time must be recorded in each run.

## Statistical model

Let \(L_s(T,\varepsilon)\) be the insertion genetic likelihood for sample \(s\)
at age \(T\) and effective discordance \(\varepsilon\), after the required
within-chromosome ARG-draw marginalisation and across-chromosome combination.
Let \(p_{\mathrm{RC},s}(T)\) be the calibrated radiocarbon-age density.

The radiocarbon-integrated likelihood contribution of training sample \(s\) is

$$
L_s(\varepsilon)
=
\int L_s(T,\varepsilon)\,p_{\mathrm{RC},s}(T)\,dT.
$$

For training set \(\mathcal S_{\mathrm{train}}\), the calibration likelihood is

$$
L_{\mathrm{train}}(\varepsilon)
=
\prod_{s\in\mathcal S_{\mathrm{train}}} L_s(\varepsilon).
$$

The first implementation will estimate one global epsilon on a prespecified
grid. It will retain the complete profile rather than only the maximum. This is
important because a flat or multimodal profile indicates that the data do not
identify a transferable epsilon.

## Primary cross-validation analysis

Use leave-one-out cross-validation when the high-confidence set is small enough.
Use stratified five-fold cross-validation when leave-one-out computation is
prohibitive. Folds should preserve age-range coverage and, where possible, avoid
placing all members of one library or sequencing batch in both training and
validation sets.

For each held-out fold:

1. Fit epsilon using only the training samples and their radiocarbon densities.
2. Retain the training epsilon profile and its uncertainty.
3. Infer each held-out sample's age from genetic data alone. Do not include its
   radiocarbon density in the age posterior.
4. Propagate training uncertainty in epsilon into the held-out age posterior,
   rather than using only a plug-in maximum-likelihood value.
5. Save de-identified fold-level diagnostics and aggregate results.

Every sample must be evaluated only in a fold where it was excluded from epsilon
calibration. An additional all-sample fit may be produced as the final operational
epsilon estimate, but it is not a validation result.

## Epsilon grid and sensitivity analysis

The initial grid should cover at least 0 to 0.15, with enough resolution near the
profile maximum to resolve differences of approximately 0.001. The final grid
must be declared before viewing held-out age errors.

Include the following prespecified comparisons:

- epsilon fixed at 0.01;
- epsilon fixed at 0.05;
- cross-validated profiled epsilon; and
- the earlier diffusion estimate near 0.06 as a labelled external reference,
  not as a fitted insertion result.

These comparisons determine whether profiling improves transfer to unseen
samples and whether the real-data result resembles the simulation behavior.

## Uncertainty propagation

Three sources of uncertainty should be kept distinct:

1. **ARG uncertainty:** marginalise posterior ARG draws within each chromosome,
   then combine independently inferred chromosomes.
2. **Genomic sampling uncertainty:** resample physical blocks within chromosome.
   Begin with the implemented 5 Mb blocks and perform a block-length sensitivity
   analysis if conclusions depend on the interval width.
3. **Calibration-sample uncertainty:** bootstrap dated individuals or use the
   variation among cross-validation training folds.

For each genomic/bootstrap replicate, refit epsilon on the training samples and
recompute the held-out age estimate. Resampling only the final age estimates
would fail to propagate calibration uncertainty.

If the fully nested bootstrap is too expensive initially, proceed in stages:

1. cross-validation with the full epsilon profile but no nested bootstrap;
2. individual-level bootstrap of epsilon;
3. genomic block bootstrap for selected folds; and
4. a nested analysis after computational cost and stability are known.

## Accuracy metrics

No single radiocarbon midpoint is exact truth. Report complementary metrics.

Point-estimate metrics:

- signed error relative to the radiocarbon posterior mean or median;
- mean bias;
- median bias;
- MAE;
- RMSE; and
- regression slope and intercept of genetic estimates against radiocarbon
  central estimates.

Interval/distribution metrics:

- fraction of genetic MAP estimates inside `cal_BP_95`;
- overlap between genetic and radiocarbon 95% intervals;
- whether the radiocarbon central estimate lies in the genetic interval;
- calibration by true-age/radiocarbon-age bin; and
- a distributional score comparing the complete genetic and radiocarbon age
  densities when full calibrated densities are available.

Report ordinary genetic posterior intervals and block-bootstrap intervals
separately. They answer different questions and must not be combined under one
generic "95% CI" label.

## Heterogeneity diagnostics

After fitting a global epsilon, test whether calibration or held-out age error
varies with:

- radiocarbon age;
- radiocarbon interval width;
- library and sequencing batch;
- coverage and missingness;
- damage or terminal-misincorporation metrics;
- population or geographic group;
- chromosome;
- estimated Ne over the relevant age interval; and
- ARG/site incompatibility counts.

These analyses should be prespecified as diagnostics, not used to repeatedly
retune epsilon on the validation samples.

If heterogeneity is reproducible, the next model should be hierarchical: a
shared epsilon distribution with sample- or library-level deviations. Avoid
fitting an unconstrained epsilon independently for every individual, which could
absorb age signal and overfit sparse samples.

## Leakage and confounding safeguards

- Ancient validation samples must not be included in the modern ARG panel.
- The held-out sample's radiocarbon age must not influence epsilon, filtering,
  site selection, or model choice for its fold.
- Folds and the epsilon grid must be fixed before inspecting held-out errors.
- Samples from the same individual or replicate libraries should remain in the
  same fold.
- Site masks must not be selected using held-out age agreement.
- Do not choose epsilon by minimizing pooled error across the validation folds;
  that would reuse the validation outcome for calibration.
- Report exclusions and failed samples before comparing methods.

Reservoir corrections, laboratory offsets, or disputed radiocarbon dates should
be resolved before defining the high-confidence set. Otherwise the calibration
will treat radiocarbon-model error as ARG/genotype discordance.

## Implementation stages

### Stage 1: input audit

- Confirm VCF/metadata identifier matches without writing protected identifiers
  into the repository.
- Confirm pseudo-haploid genotype encoding.
- Verify calibrated-age units and interval ordering.
- Summarise the number of eligible samples across the age range.
- Record library/batch groupings needed for fold construction.

### Stage 2: insertion epsilon profiling

- Add a mode that evaluates the insertion likelihood over a two-dimensional
  \((T,\varepsilon)\) grid without recomputing ARG/tree quantities unnecessarily.
- Integrate each training sample over its radiocarbon density.
- Combine training samples to produce the fold-specific epsilon profile.
- Save sufficient diagnostics to count incompatibilities and identify whether a
  small number of sites dominates the profile.

### Stage 3: cross-validated age inference

- Generate deterministic, age-stratified folds.
- Fit epsilon on each training set.
- Infer held-out genetic ages without radiocarbon information.
- Produce aggregate accuracy and calibration summaries.

### Stage 4: uncertainty and heterogeneity

- Add individual and genomic-block bootstraps.
- Propagate epsilon uncertainty into held-out age estimates.
- Test prespecified sample/library covariates.
- Decide whether a global or hierarchical epsilon model is warranted.

### Stage 5: operational fit

After the validation design and model are frozen, fit epsilon using all eligible
high-confidence radiocarbon samples. Record the full profile, uncertainty, ARG
draws, Ne input, block scheme, sample inclusion criteria, and software revision.
Use that fit for samples lacking reliable dates.

## Planned outputs

Repository-safe outputs should include:

- aggregate epsilon profile by fold;
- cross-validated accuracy table without protected sample identifiers;
- genetic-versus-radiocarbon plot with radiocarbon uncertainty;
- residuals against age, coverage, damage, batch, and missingness;
- posterior and block-bootstrap coverage summaries;
- sensitivity comparison for fixed 0.01, fixed 0.05, and profiled epsilon;
- machine-readable run metadata containing code revision and model settings; and
- a concise list of exclusions and their reasons.

Protected sample-level inputs and identifiable result tables remain outside the
repository.

## Decision criteria

A single global insertion epsilon is acceptable only if:

- fold-specific profiles have overlapping, interior maxima;
- held-out bias does not change direction strongly across age or major library
  groups;
- profiled epsilon improves or at least preserves held-out MAE/RMSE relative to
  prespecified fixed values;
- improvement is not merely cancellation among subgroups; and
- uncertainty intervals have useful, reproducible calibration.

If those conditions fail, do not tune a global epsilon further. Determine
whether the failure tracks genotype quality, POLEGON/SINGER incompatibility,
estimated demography, or radiocarbon uncertainty, then introduce the smallest
model extension supported by the diagnostics.
