# Previous approaches (before the non-sharing ratio)

The repository's front page now describes the non-sharing-ratio estimator. Earlier approaches estimated age from
allele-frequency models rather than from non-sharing. They are not the current direction and are kept here for reference.
Last updated 2026-11-13.

| Approach | What it did | Where |
|---|---|---|
| 1. Present-count diffusion (the former `main`) | Age posterior from derived-allele frequency at the sample's age, conditioned on mutation age and present panel count, under a Wright-Fisher moment recursion | The old front page is preserved as [PRESENT_COUNT_DIFFUSION_README.md](PRESENT_COUNT_DIFFUSION_README.md). Code: `pipeline/posterior_sample_age_infer.py`, `pipeline/precompute_freq_trajectory_moments.py`. Derivation: MATH.md sections 3-6. Record: [../working_diffusion.md](../working_diffusion.md). Last state of the old `main`: commit `bf27303` (branch `insertion` points at it) |
| 2. Tree-at-$T$ beta-binomial | Beta-binomial likelihood from the tree at the sample's age | Branch `betabinom` (archived); MATH.md section 7; [../working_betabinom.md](../working_betabinom.md) |
| 3. Ancient-lineage insertion | Insert the ancient lineage into the modern ARG at age $T$; the product over sites gives a posterior for $T$. Point-calibrated in simulation with a known ARG; intervals overconfident under linkage | Branch `insertion`; [../working_insertion.md](../working_insertion.md); MATH.md section 8; graphical abstract [graphical_abstract_readme.png](graphical_abstract_readme.png) |
| 3a. Insertion validation with POLEGON | POLEGON-dated ARG insertion calibration and its validation guide | Branch `feature/polegon-insertion-validation`; [../POLEGON_INSERTION_VALIDATION.md](../POLEGON_INSERTION_VALIDATION.md) |
| 3b. Insertion with radiocarbon epsilon calibration | Joint ARG time-scale and epsilon calibration against radiocarbon ages | Branch `feature/radiocarbon-epsilon-calibration`; [../TIME_SCALE_EPSILON_CALIBRATION_PLAN.md](../TIME_SCALE_EPSILON_CALIBRATION_PLAN.md) |

Status in the old README (2026-11): approach 1 was the production pipeline; approach 3 was the preferred direction, in
development, until it consumed the production SNP/ARG store, handled posterior ARG draws and chromosomes, and passed an
end-to-end real-data-shaped simulation (see [../PROJECT_MAP.md](../PROJECT_MAP.md), "Branch and production policy").

## Issues closed when the project moved to the non-sharing ratio (2026-11-13)

These were open against approach 1 or the insertion approach. They were closed as not planned, not deleted, and can be reopened.

| Issue | Title |
|---|---|
| #1 | Insertion likelihood: pseudo-haploid only, but README doesn't say so and code doesn't enforce it |
| #2 | Approach 1: heterozygotes silently become derived under --ploidy 1 (add --het-policy, counter, warning) |
| #3 | Approach 1: check moment inequalities and genotype probabilities under --ploidy 2 instead of clipping |
| #4 | End-to-end msprime validation of E[p_T \| d0, t_i] and the second moment |
| #5 | Insertion: localise age bias with inferred SINGER+POLEGON ARGs and inferred Ne |

Open work on the current approach is tracked in GitHub issues with the label `non-sharing-ratio`.
