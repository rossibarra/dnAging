# Joint ARG Time-Scale and Epsilon Calibration Plan

## Objective

Extend the ancient-lineage insertion analysis so that radiocarbon-dated samples
jointly estimate:

1. a multiplicative correction, \(\varsigma\), for systematic compression or expansion
   of inferred ARG times; and
2. an effective discordance probability, `epsilon`, that prevents isolated
   ARG/genotype contradictions from assigning an age zero likelihood.

The calibrated model will then estimate ages for samples without trusted dates.
Both parameters must be estimated using only the training samples within each
cross-validation fold.

The immediate motivation is the inferred-ARG simulation result. Profiling only
`epsilon` selected `epsilon = 0.01`, but produced a mean age bias of about -970
generations and a regression slope of approximately 0.81. A fixed
`epsilon = 0.05` reduced the MAE substantially, but `epsilon` should not be
forced to compensate for a systematic ARG time-scale error.

## Model

Let:

- `u` be time in the inferred ARG;
- \(t = \varsigma u\) be corrected time in generations;
- `T_i` be the true sampling age of ancient sample `i`;
- `G_i` be its observed genotype calls;
- `R_i(T)` be its calibrated radiocarbon-age density; and
- `epsilon` be the effective per-allele discordance probability.

The initial model will use a multiplicative correction:

\[
t = \varsigma u, \qquad \varsigma > 0.
\]

This preserves time zero for modern samples. An affine transformation should
not be the default because a nonzero intercept would move the present. A more
flexible, monotone time warp can be considered later only if a multiplicative
scale leaves clear age-dependent residual structure.

For a dated training set `D`, estimate the joint calibration surface with:

\[
L(\varsigma,\epsilon)
=
\prod_{i\in D}
\int L_i(T\mid \varsigma,\epsilon)R_i(T)\,dT.
\]

Here, \(L_i(T\mid\varsigma,\epsilon)\) is the insertion likelihood for the genotype data
of sample `i`. The radiocarbon density is integrated rather than replaced with
a point estimate. On a discrete age grid, the log likelihood is:

\[
\ell(\varsigma,\epsilon)
=
\sum_{i\in D}
\log\sum_k
\exp\{\ell_i(T_k\mid \varsigma,\epsilon)\}
R_i(T_k)\,\Delta T_k.
\]

This expression asks which combination of ARG time scale and discordance rate
makes the genetic likelihood agree with the known-age distributions.

For an unknown or held-out sample `j`, propagate calibration uncertainty rather
than using only a point estimate:

\[
p(T_j\mid G_j,D)
\propto
p(T_j)
\sum_{\varsigma,\epsilon}
L_j(T_j\mid \varsigma,\epsilon)
p(\varsigma,\epsilon\mid D).
\]

The weights \(p(\varsigma,\epsilon\mid D)\) can initially be grid-based likelihood weights
with explicit priors. A bootstrap approximation will be used for reported
uncertainty.

## Interpretation and identifiability

\(\epsilon\) and \(\varsigma\) address different failure modes:

- `epsilon` softens isolated conflicts, such as an ancient derived allele whose
  inferred mutation is younger than the sample. It is an effective discordance
  rate, not necessarily a sequencing-error rate.
- \(\varsigma\) corrects a systematic distortion of inferred node and mutation times. A
  value above one stretches an ARG whose inferred times are too young.

The parameters will be correlated: increasing `epsilon` can make older ages
possible by weakening too-young mutation conflicts. A broad set of known
samples spanning the age range is therefore important. The joint likelihood
surface and its boundary behavior must always be reported.

Without known-age samples or another external time anchor, \(\varsigma\) is generally not
identifiable from ancient genotypes alone because ARG time and sample age can
move together. `epsilon` may be weakly informed by discordance, but it remains
confounded with sample age, genotype error, polarization error, topology error,
and mutation-time error. Joint estimation of both parameters from only unknown
samples should not be treated as reliable.

Modern samples anchor time zero but do not adequately calibrate the deeper time
scale. A trusted mutation rate, recombination rate, or demographic history may
provide partial scale information, but the empirical correction considered here
still requires dated samples for validation.

## Correct application of the time scale

The final implementation must apply \(\varsigma\) inside the insertion likelihood:

- multiply inferred ARG node times and mutation-edge endpoints by \(\varsigma\);
- leave the candidate ancient age grid in real generations;
- evaluate the coalescent hazard using the real-generation Ne history; and
- apply `epsilon` only in the genotype-emission probability.

Merely multiplying the final MAP ages by a fitted constant is not the final
model. That shortcut does not correctly update coalescent exposure, likelihood
curvature, or age intervals. A post-hoc axis transformation may be used as a
fast diagnostic, but its output must be labeled approximate.

## Proposed code changes

### 1. Scale-aware insertion kernel

Modify `pipeline/insertion_likelihood.py`:

- Add a `time_scale=1.0` argument to the insertion-probability functions.
- Apply the scale consistently to all tree node times and mutation-edge
  endpoints used by the kernel.
- Continue evaluating `PiecewiseConstantNe` in real generations.
- Validate that every scale is finite and positive.
- Preserve exact existing behavior when `time_scale == 1`.

Avoid mutating the tree sequence or writing one rescaled tree per grid point.
The scale should be applied numerically while evaluating the kernel.

### 2. Joint likelihood grid

Modify `pipeline/insertion_real_data.py`:

- Add `--time-scale-grid`, accepting `START:STOP:STEP` or comma-separated
  values, analogous to `--epsilon-grid`.
- Accumulate a likelihood tensor with logical dimensions:

  ```text
  time_scale x epsilon x chromosome x sample x age
  ```

- Store the time-scale grid in `time_scale_grid.npy`.
- Store summed likelihoods per draw in
  `ll_by_chrom_scale_epsilon.npy`.
- Include both grids in checkpoint compatibility metadata.
- Extend checkpoint state to include the scale/epsilon likelihood tensor.
- Retain the current \(\varsigma=1\) outputs for backward compatibility when no scale
  grid is requested.

The first diagnostic grid should be deliberately coarse, for example:

```text
time scale: 0.80 to 1.60 by 0.05
epsilon:    0.00 to 0.12 by 0.01
```

After locating the joint maximum, refine the scale grid locally. The simulation
slope of approximately 0.81 suggests \(\varsigma\) near
\(1 / 0.81 = 1.24\), but this is an
initial diagnostic rather than a fixed correction.

### 3. Runtime and memory control

The scale dimension increases computation substantially. Implement the
following before a full real-data grid:

- vectorize probability calculations across scale values where possible;
- process scales in configurable chunks;
- stream likelihood results to checkpoint files instead of retaining all
  chromosomes and draws in memory;
- retain atomic checkpoints and record the next chromosome, target row, and
  scale chunk;
- log target-row and scale-chunk progress for ETA estimation; and
- profile whether tree seeking, probability calculation, or likelihood
  accumulation dominates runtime.

Do not launch all ten real-data ARG draws on a dense scale grid until a single
chr1 draw has established runtime and memory requirements.

### 4. Merge across inferred ARG draws and chromosomes

Extend the merge path in `pipeline/insertion_real_data.py`:

- marginalize ARG draws independently within each chromosome;
- then sum independent chromosome likelihoods;
- preserve both the scale and epsilon axes; and
- validate identical grids, samples, retained sites, and chromosome layouts
  across draws.

The merged output should contain:

```text
grid.npy
time_scale_grid.npy
epsilon_grid.npy
ll_marginal_scale_epsilon.npy
samples.txt
run.json
```

### 5. Joint radiocarbon calibration

Extend `pipeline/radiocarbon_calibration.py` or add a narrowly scoped
`pipeline/joint_scale_epsilon_calibration.py`:

- Read the merged scale/epsilon/age likelihood tensor.
- Integrate each training sample over its full radiocarbon density.
- Add explicit priors for \(\varsigma\) and \(\epsilon\); use uniform grid priors initially
  and record them in metadata.
- Normalize and save the two-dimensional joint calibration surface.
- Marginalize over the joint calibration surface when dating held-out and
  unknown samples.
- Report marginal and joint MAP estimates, credible intervals, boundary mass,
  and parameter correlation.

Suggested outputs:

```text
joint_scale_epsilon_profile.tsv
scale_marginal.tsv
epsilon_marginal.tsv
cross_validated_ages.tsv
unknown_ages.tsv
calibration_bootstrap.tsv
run.json
```

### 6. Leakage-safe cross-validation

Use age-stratified K-fold cross-validation:

1. Divide dated samples into folds spanning the radiocarbon-age range.
2. Estimate the complete \(p(\varsigma,\epsilon\mid\text{training dates})\) surface using only
   training samples.
3. Date the held-out samples without using their dates.
4. Repeat for every fold.

Report:

- mean bias;
- MAE and RMSE;
- regression slope and intercept of estimated versus radiocarbon age;
- genetic interval coverage of radiocarbon means;
- overlap with radiocarbon 95% intervals;
- metrics by age bin; and
- the fitted \(\varsigma\) and \(\epsilon\) for every fold.

The goal is not merely to maximize the training-date likelihood. The chosen
model must improve held-out predictive accuracy and remove the age-dependent
residual trend.

### 7. Bootstrap uncertainty

Use two complementary bootstraps:

1. **Dated-sample bootstrap:** resample known individuals and re-estimate the
   joint \((\varsigma,\epsilon)\) surface. This measures sensitivity to the calibration
   panel.
2. **Genomic block bootstrap:** resample linked genomic blocks and propagate
   the result through calibration and unknown-age estimation. This measures
   genomic sampling variability and linkage sensitivity.

The final interval should propagate both sources. A practical first
implementation can use nested resampling with fewer inner replicates, followed
by a larger run once runtime is established. ARG-draw uncertainty must remain
separate from block and calibration-panel uncertainty.

## Tests

Add focused unit and integration tests:

1. `time_scale = 1` exactly reproduces the current insertion likelihood.
2. Scaling node and mutation times by \(\varsigma\) matches an explicitly rescaled tiny
   tree sequence.
3. Simulated data with a known imposed scale recovers that scale.
4. Simulated isolated contradictions increase estimated `epsilon` without
   inducing a scale shift.
5. Simulated pure time compression changes \(\varsigma\) without requiring inflated
   `epsilon`.
6. Jointly perturbed simulations recover both parameters within uncertainty.
7. Cross-validation never includes held-out dates in the calibration surface.
8. Marginalization over a one-element scale grid matches the existing
   epsilon-only implementation.
9. Checkpoint/restart produces byte-equivalent numerical arrays to an
   uninterrupted run.
10. Merge order does not affect the result, and incompatible grids fail loudly.
11. Bootstrap output is deterministic under a fixed seed.
12. Boundary diagnostics trigger when the scale or epsilon maximum is at the
    edge of its grid.

## Validation sequence

### Phase A: Fast simulation diagnostic

- Apply an approximate post-hoc likelihood-axis scaling to the completed
  inferred-ARG simulation output.
- Determine whether a single multiplicative scale removes most of the observed
  slope and age-dependent bias.
- Use this only to choose a sensible exact scale grid.

### Phase B: Exact simulation validation

- Recompute the inferred-ARG insertion likelihood over the joint grid.
- Use the same 20 known simulated samples and the same held-out design.
- Require materially improved held-out MAE, bias, regression slope, and
  interval coverage compared with epsilon-only calibration.
- Compare inferred-ARG results with true-ARG controls to separate inference
  error from insertion-model error.

### Phase C: chr1 real-data pilot

- Run one ARG draw on chromosome 1 with checkpoints.
- Estimate \((\varsigma,\epsilon)\) from the 50 trusted radiocarbon samples.
- Produce explicitly provisional ages for the five uncertain samples.
- Inspect the joint surface for ridges and boundary maxima.

### Phase D: ARG-draw and genome expansion

- Run all selected ARG draws for chromosome 1.
- Add chromosomes incrementally after verifying retained-site consistency.
- Repeat cross-validation and all uncertainty calculations after each
  expansion.

## Acceptance criteria

Before using calibrated real-data ages as primary results:

- the simulation joint maximum is not pinned to a grid boundary;
- the time-scale estimate is stable across cross-validation folds;
- held-out regression slope is close to one and intercept close to zero;
- age-dependent residual bias is substantially reduced;
- cross-validated MAE and RMSE improve over epsilon-only calibration;
- bootstrap intervals have reasonable simulated coverage;
- estimates are not dominated by one chromosome, one ARG draw, or a small
  number of dated samples; and
- results remain qualitatively stable under reasonable priors and grid
  refinements.

## Important data requirements

- Use full calibrated radiocarbon densities whenever available, on the same
  generation-before-present scale as the genetic age grid.
- Record the generation-time conversion explicitly if calendar years and
  generations differ.
- Treat the five uncertain samples as completely unavailable during parameter
  calibration.
- Preserve sample IDs exactly as they appear in the likelihood arrays.
- Record which dates, calibration curves, ARG draws, Ne history, mutation rate,
  chromosome masks, and genotype filters were used in every run.

## Deliverables

1. Scale-aware insertion likelihood code with checkpoint support.
2. Joint scale/epsilon merge and calibration code.
3. Unit and simulation-recovery tests.
4. Cross-validation and bootstrap summaries.
5. A simulation comparison of epsilon-only versus joint calibration.
6. Provisional chr1 ages for the five uncertain samples.
7. Final multi-draw, multi-chromosome calibrated age estimates after the
   simulation acceptance criteria are met.
