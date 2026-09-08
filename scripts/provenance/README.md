# Provenance scripts for the Ne-sweep results

These four scripts produced the Ne sweep that `working_diffusion.md` discusses,
including the figures that were in the removed `bias_ideas.md`. They lived only
in `/tmp` and are committed here for reconstruction, **not** as maintained tools.

They are recorded as-found. Do not run them unmodified: each hardcodes absolute
local paths (`/Users/jeffreyross-ibarra/src/dnAging`, `~/Projects/mutrates/...`,
`/tmp/...`) and several read a `haploid_simulations_100mb_batch 2` directory for
ages and seeds that is not in this repository. On a cluster those paths must be
supplied; most are already overridable by environment variable.

## What each did

**`simulate-ne10k-infinite-batch.py`** — generated the simulations under
`/tmp/dnAging-ne{Ne}-infinite-10mb/simulation_{01..10}`. msprime, 26 modern
haploid samples at time 0 plus one ancient haploid at the replicate's true age;
`population_size=NE`, `ploidy=2`, Hudson model; mutation and recombination rates
both 1e-8; `InfiniteSites(NUCLEOTIDES)` with `discrete_genome=False`. Sites are
restricted to those polymorphic in the 26 modern samples, then modern and ancient
tree sequences are simplified separately so **the published ARG contains only the
modern panel**. Provenance and metadata are stripped from the modern trees. Ages
and seeds are reused from an earlier 100 Mb batch so replicates are comparable
across sweeps; the simulations themselves were regenerated. Writes `run.json`
with the true ages.

**`build-ne10k-moment-table.py`** — built the frequency tables for that sweep.
Note it uses the **legacy float64 `MomentEngine`** and its `max_cancellation`
guard, writing `NaN` where the guard fires, and emits **format 4**. So the sweep
tables predate the exact-arithmetic work in d75b9b9 / 54b9113 / 33b23d6. T4
measured what that costs on this data: -0.2 generations, so the old tables are
not the source of the bias.

**`dnAging-bigsims-estimate.py`** — the inference harness, and the important one
to understand. **It marginalises the mutation age by DENOMINATOR WEIGHTING**, not
uniformly: `denominator_rows` recomputes the alternating-sum denominators from
`MomentEngine._moms` and `eng.coeff`, then integrates `den * phi` over `den`
by dense trapezoid. Its own output records `"conditioning": "ratio of
integrals"`. It imports production `load_table` but **never calls production
`phi_lookup`, normalizeTEs, or the ARG-draw pipeline**, so its numbers never
described the shipping uniform estimator. This is why `bias_ideas.md`'s table
could not be reproduced from any production code path, and why it was removed —
see `working_diffusion.md`, "State of the bias". Reproducing its weighted result
independently (format5 `log_den` plus analytic per-knot integration) agrees to
within one generation, which validates both implementations.

Also supports options used in the singleton experiments: `SIM_MIN_DERIVED_FREQ`,
`SIM_DROP_RARE_CARRIED`, `SIM_EPSILONS` (the sweep used `0`), and a 100-block
bootstrap.

**`summarize-batch2-posteriors.py`** — turned the saved posteriors into
bias/RMSE/coverage. Reports the **posterior mean**, with HPD95 by shortest
interval. Worth knowing when comparing against MAP-based numbers elsewhere.

## Reproducing the current baseline instead

For the uniform estimator that main actually ships, use
`scripts/t4_attribution_sweep.py`, which evaluates production `phi_lookup`
directly on the true tree sequences and can select either marginalisation.
