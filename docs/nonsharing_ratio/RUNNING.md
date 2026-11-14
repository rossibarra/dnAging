# Running the non-sharing-ratio estimator on real data

Inputs you supply:

| Variable | What |
|---|---|
| `CONTROLS_VCF` | new VCF of held-out modern genomes, NOT in the ARG panel (phased diploid -> `CONTROLS_MODE=phased`, the default; haploid or unphased -> `CONTROLS_MODE=sample`) |
| `ANCIENT_VCF` | ancient pseudo-haploid VCF (`0/0`,`1/1`,`./.` or haploid calls), same chromosome labels as the ARG (`1`..`10` or `chr1`..; accession labels such as `NC_050096.1` are not mapped) |
| `DRAW_DIR` | directory with the ARG draws (`*.tsz` / `*.trees` directly inside; draw id = trailing number in the name). Real raw SINGER draws: `/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/results/combined/tsz_files` |
| `CHROM_FAI` | `/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/ref.fa.fai` (the draws carry no `chrom_offsets` metadata) |
| `A0` | cutoff in generations, above the oldest sample age, in the ARG's time units (you must choose; see the A0 issue and report section 6) |
| `OUTROOT` | NEW directory for outputs (nothing is overwritten; a second run into the same files fails) |
| optional `NE_TSV` (+ `NE_SERIES`, default `posterior_mean`) or `DEMES_YAML` | c(t) source; default is the ARG's own lineage counts |

Checks before submitting: control sample names are not in the ARG panel nor the ancient VCF (pass `PANEL_NAMES=file` with the panel's
sample names to have the CLI verify it); controls and ancient VCFs use `REF`/`ALT` consistently with the panel's genome (polarity is by allele string, so flipped REF/ALT is fine).

```bash
cd /quobyte/jrigrp/beil/logan_collab/dnAging
export OUTROOT=$PWD/results/nonsharing_ratio/real_YYYYMMDD     # new directory
export CONTROLS_VCF=/path/to/new_controls.vcf.gz
export ANCIENT_VCF=/path/to/ancient_pseudohaploid.vcf.gz
export DRAW_DIR=/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/results/combined/tsz_files
export CHROM_FAI=/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/ref.fa.fai
export A0=...                                                  # generations
export NE_TSV=$PWD/radiocarbon_goodvcfs/insertion_pilot_20261012/demography/posterior_mean.from_zero.tsv   # optional

# 1. ONE draw first; check time and memory
sbatch --export=ALL --array=0-0 slurm/nonsharing_ratio_run_draws.sbatch
sacct -j <JOBID> --format=JobID,Elapsed,MaxRSS,State    # when finished; then read nsrDraw-<JOBID>_0.out

# 2. all draws (N = number of *.tsz in DRAW_DIR; 75 for the real directory); skip index 0 if already done
sbatch --export=ALL --array=1-74%10 slurm/nonsharing_ratio_run_draws.sbatch

# 3. merge after the array finishes (writes $OUTROOT/ages_merged.tsv)
sbatch --export=ALL --dependency=afterok:<ARRAY_JOBID> slurm/nonsharing_ratio_merge.sbatch
```

Notes
* Per task: 1 CPU, 100 GB, 6 h, account `jrigrp`, partition `high` (edit the `#SBATCH` lines, or pass `-p low` on the sbatch line; `low` can preempt and requeue, and a requeued task fails because it refuses to overwrite its half-written output: delete `draws/draw_<id>.npz` first).
  These resources are an extrapolation from a chr1 run (21 GB, 11 min): confirm with the single-draw run.
* Each task re-reads both VCFs (pure Python): expect about an hour per draw for a genome-wide ancient VCF plus controls.
* Output: `$OUTROOT/draws/draw_<id>.npz` (per-block counts) and `$OUTROOT/ages_merged.tsv` with columns
  `sample T_hat_mean T_hat_median draw_sd ci_lo ci_hi n_draws flag`. `flag` lists degenerate targets per draw count:
  `f0` (no non-sharing at all: T_hat = 0 by convention), `H<=0` (no excess non-sharing: T_hat truncated to 0), `no_data`. Do not interpret flagged rows as ages.
* Panel-fixed sites outside the ARG (optional; real data currently lack an invariant-site panel table and an ancestral-allele source):
  `python pipeline/nonsharing_ratio_real_data.py prepare-fixed --panel-vcf PANEL_WITH_INVARIANT_SITES.vcf --output $OUTROOT/panel_fixed.npz`
  (or `sbatch slurm/nonsharing_ratio_prepare_fixed.sbatch`), then add `PANEL_FIXED=$OUTROOT/panel_fixed.npz POLARITY_TABLE=anc_alleles.txt`
  (`chrom pos ancestral_allele` per line) to the draw job. Without an ancestral-allele table this option is refused (REF is not ancestral in the maize ARG).
* Smoke test on a chromosome subset: add `CHROMS=1` (restricts the ARG, the VCF reads and c(t) from the ARG to that chromosome).
* Without held-out controls the estimator is invalid (design rule 1). Using the panel as controls is a plumbing test only.
