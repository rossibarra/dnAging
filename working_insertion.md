# Working record: direct ancient-lineage insertion

This document is the experimental and implementation record for the direct
ancient-lineage insertion likelihood (MATH.md, approach 3). Diffusion-specific
diagnostics remain in `working_diffusion.md`; beta-binomial work remains in
`working_betabinom.md`.

## Established controls

The insertion likelihood is point-calibrated when the modern ARG and demography
are known. On the 300-replicate constant-size benchmark it had MAP bias +12.9
generations and slope 1.001. On 1,000 estimates under independently varying,
piecewise-constant demographies with correctly specified 1% observation error,
MAP bias was +7.7 generations. These controls show that the core insertion
likelihood, variable-Ne integration, and nonzero-error transformation do not by
themselves produce sustained age bias.

Linked-site posterior intervals remain overconfident. The implemented 5 Mb block
bootstrap improves coverage but does not make it nominal in all settings.

### Completed three-epoch factorial

The earlier three-chromosome, three-epoch simulation already compared true ARG +
true Ne, POLEGON-dated SINGER ARG + true Ne, and POLEGON-dated SINGER ARG +
estimated Ne on the same ten ancient samples. This established that inferred-ARG
error is the primary source of the large zero-epsilon failure, with estimated Ne
modifying rather than creating it.

The true-ARG/true-Ne, epsilon-zero condition had mean bias +12.6 generations,
RMSE 57.1, and 10/10 interval coverage. With the estimated ARG, epsilon zero was
catastrophic: mean bias was -2,232 with true Ne and -2,217 with estimated Ne.
At epsilon 0.01, mean bias was -237 and RMSE 289 with true Ne, compared with
-362 and RMSE 476 with estimated Ne. Thus ARG inference/dating generated the
main discrepancy in that dataset, while Ne estimation amplified it.

This factorial does not need to be rediscovered. Repeating its missing arms on
the nine-demography suite would test whether that attribution remains stable
across constant, bottleneck/recovery, and recent-growth histories.

## Inferred-ARG, inferred-Ne benchmark

### Design

The nine-demography suite is under
`batches/nine_demographies_20260916/`. It contains three independent runs from
each of three demographic models: constant size, bottleneck/recovery, and recent
growth. Every run has three 50 Mb chromosomes and 21 ancient samples at true ages
0, 500, ..., 10,000 generations.

The modern data were analysed with `singer-snakemake`. This means the inferred
ARGs combine **SINGER topologies with POLEGON node/branch dating**; they must not be
described as SINGER-only ARGs. There are 100 posterior draws per chromosome, and
age inference uses draws 50--99. ARGtest supplies the estimated Ne trajectory.
Inference marginalises ARG draws within chromosome, combines chromosomes after
that marginalisation, and uses a 5 Mb block bootstrap.

The simulated ancient genotypes contain **no observation error**:
`epsilon_applied = 0`, zero flipped alleles. Consequently, positive epsilon in
the analyses below is deliberately misspecified as genotyping error. It acts as
an effective discordance or robustness parameter for error in the POLEGON-dated
SINGER ARG, polarity, and other model mismatch.

### Results

All nine analyses completed with 50/50 retained ARG draws and successful merges
at both epsilon values.

| Epsilon | Mean bias | MAE | RMSE | 5 Mb bootstrap coverage | Posterior coverage |
|---:|---:|---:|---:|---:|---:|
| 0.01 | -592.9 | 616.7 | 769.1 | 55/189 (29.1%) | 11/189 (5.8%) |
| 0.05 | +111.9 | 363.2 | 452.8 | 82/189 (43.4%) | 28/189 (14.8%) |

The epsilon=0.05 improvement in the pooled mean is not uniform calibration. It
is cancellation among large, demography-specific errors:

| Demography/run | Mean bias | MAE | RMSE | 5 Mb bootstrap coverage |
|---|---:|---:|---:|---:|
| constant r1 | +4.8 | 85.7 | 125.4 | 19/21 |
| constant r2 | +40.5 | 154.8 | 180.9 | 20/21 |
| constant r3 | +59.5 | 173.8 | 211.9 | 16/21 |
| bottleneck/recovery r1 | -245.2 | 345.2 | 393.5 | 7/21 |
| bottleneck/recovery r2 | -269.0 | 369.0 | 415.6 | 7/21 |
| bottleneck/recovery r3 | -338.1 | 385.7 | 454.3 | 7/21 |
| recent growth r1 | +540.5 | 540.5 | 589.9 | 1/21 |
| recent growth r2 | +633.3 | 633.3 | 690.8 | 3/21 |
| recent growth r3 | +581.0 | 581.0 | 627.0 | 2/21 |

At epsilon=0.01, all three model classes were biased young: constant runs by
-771 to -802 generations, bottleneck/recovery by -410 to -500, and recent growth
by -505 to -588. Raising epsilon to 0.05 nearly removes constant-size bias,
reduces but does not remove young bias after bottleneck/recovery, and reverses the
recent-growth error into a large old bias.

### Interpretation

Five percent effective discordance is a better pooled robustness setting than
one percent for this benchmark, but it does **not** solve calibration. Because
the data were generated with epsilon=0, the apparent optimum cannot be
interpreted as a genotyping-error estimate. A single global epsilon is absorbing
demography-dependent errors in inferred ARG topology/dating and estimated Ne.
The near-zero pooled bias at epsilon=0.05 therefore conceals misspecification.

The current inferred-ARG condition cannot attribute the problem specifically to
SINGER: `singer-snakemake` ran POLEGON, so topology error and node-time error are
confounded. Nor can it attribute all remaining variation to the ARG, because
ARGtest's Ne trajectory is used in the same run.

## Next diagnostic sequence

1. **Separate SINGER topology from POLEGON dating.** Retain the inferred SINGER
   topology while changing its dating treatment, and compare POLEGON node/edge
   time distributions with truth by genomic region and true sample age. Establish
   first whether `use-polegon: false` yields branch lengths on a usable generation
   scale; do not treat that switch as a valid control without verifying units.
2. **Localise the likelihood shift.** At each true sample age, measure the
   log-likelihood difference between the true age and MAP by chromosome, ARG
   draw, derived-carrier sites, ancestral-carrier sites, and sites incompatible
   at the true age. Relate these contributions to edge child/parent times and
   inferred-versus-true mutation-bearing clades.
3. **Profile effective epsilon instead of selecting it by pooled age RMSE.** Fit
   epsilon within each simulated demography/run and inspect whether its profile
   is stable across age and chromosome. Since true observation epsilon is zero,
   any positive fitted value directly measures aggregate model discordance. A
   demography-specific optimum would confirm that a universal epsilon is not a
   transportable error parameter.
4. **Optionally extend the completed factorial across demographies.** Add true
   ARG + true Ne and estimated ARG + true Ne arms to the nine-demography suite if
   the SINGER-versus-POLEGON diagnostics do not already identify the mechanism.
   This is replication and interaction testing, not an unperformed foundational
   control.
5. **Keep point calibration and interval calibration separate.** Continue the
   5 Mb block bootstrap, but do not tune epsilon to improve coverage. After the
   point-bias source is isolated, evaluate block-length sensitivity or a
   replicate-clustered/Godambe interval calculation.

The immediate priority is steps 1 and 2. Additional epsilon values alone cannot
identify the cause; they only trace the response to increasingly downweighted
ARG/genotype incompatibilities.
