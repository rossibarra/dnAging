# Non-sharing-ratio estimator on SINGER+POLEGON ARGs of msprime simulations (realsims, 2026-11-12)

## Design (as specified)
* 5 independent replicates, constant Ne = 20,000 (msprime `ploidy=1`, `population_size=2*Ne`, i.e. pair coalescence rate 1/(2Ne)),
  mu = 1.25e-8, rho = 1e-8, 20 Mb, all haploid, BinaryMutationModel, discrete genome (same conventions as `validation/nonsharing_ratio_validation.py`).
* Per replicate: 26 modern panel haplotypes (the only samples the ARG is built on) + 10 held-out modern controls (not in the ARG) +
  10 ancient haplotypes (46 haplotypes: 36 modern + 10 ancient, as the user specified). 50 ancient samples in total.
* Seeds: ancestry 20261112+100*rep, mutations +1, ancient ages RNG (numpy default_rng) +2. Ancient ages: stratified uniform, one draw in each of
  ten 600-generation bins over [0,6000] (exact ages in `sim/rep*/truth.tsv`).
* VCFs (CHROM `1`, POS = msprime position + 1, REF = ancestral `A`, ALT = derived `C`): `panel_full.vcf.gz` (all sites, pre-dropping),
  `panel_segregating.vcf.gz` (sites segregating in the panel; SINGER input), `controls.vcf.gz` (all sites), `ancient.vcf.gz` (haploid GT; unique-to-ancient sites removed).
* Ancient ascertainment: sites whose derived allele is in no modern haplotype were removed from the ancient VCF. Removed per replicate (rep1..5):
  10,645 / 10,233 / 10,321 / 10,482 / 10,220 sites (of ~93-95k simulated; >99.8% of them carry the derived allele in at least one ancient sample, the rest are in nobody).
* SINGER 0.1.8-beta custom + POLEGON 0.1.3-alpha via the singer-snakemake workflow, same settings as `batches/nonsharing_ratio_heldout_20261111`
  (150 MCMC samples, thin 10, POLEGON 150/10/5 reps, 1 Mb chunks), `polarised: False` (sites randomly polarised and re-polarised by SINGER as in the real data, although REF is truly ancestral here).
  mu/rho set to the simulation values; no masks. 5 draws per replicate, evenly spaced post-burn-in: 75, 93, 112, 131, 149 (of 0..149; burn-in 75).
  This is the real-data recipe, not the older direct-run recipe of `singer_polegon_pilot_20261102`.
* Estimator: `run-draw --age-mode child --a0 {10000,20000,50000}`, held-out controls (`--controls-mode sample`), `--panel-names`, `--panel-fixed` from
  `prepare-fixed --min-called 21` on `panel_full.vcf.gz`, `--fixed-polarity ref`, c(t) from the ARG (default; grid max 20000). Merge: mean of T_hat over the 5 draws (also ran `merge` -> `merged/rep*_a0*.tsv`, which includes bootstrap CIs; the table/figure use the draw mean computed by `05_analyze.py` from the draw NPZs, equal to merge's `T_hat_mean`).

## Deviations / choices to know
* None from the 46-haplotype design. Choices not fixed by the brief: sampling ages stratified; `polarised: False`; 12 CPUs/job; draw ids above; hapmap with constant 1 cM/Mb added so the sequence length is 20,000,000 (+1 in the ARG) and `--chrom-fai` (`chr1.fai`) reconstructs it.
* Failed attempt: job 40594403 (simulation rep 1) crashed building the true-tree export (mutation time older than parent after simplify); fixed by setting mutation times unknown in that export only; no other effect.

## Verified checks
* ARG time scale (observed): SINGER's own Ne estimate 19,990 (config 20,000); sequence length 20,000,001; 26 samples; tree sequence diversity 9.99e-4 (=4 Ne mu); root time per tree median 82.5k (draw 75) / 77.5k (draw 149) generations vs expectation ~77k; max node time ~5e5. So times are in generations consistent with the simulation.
* POLEGON-dated, not raw SINGER (observed): the workflow's merge rule reads `chunks/*_polegon_*.txt` (100 files for rep 1 = 5 draws x 20 chunks) and `use-polegon` is not disabled; node times in `_polegon_` differ from raw `_nodes_` files.
* (c) No ancient/control in the ARG: SINGER input VCF has exactly P01..P26 (checked for all 5 reps); every ARG has 26 samples; the CLI's disjoint-name check against `panel_names.txt` passed.
* (d) Position/allele convention: checked on draw 75 (A0 10000) of each replicate (NPZ meta): rep1 76,298 of 76,307 ARG sites matched, 0 allele mismatches; reps 2-5: 76,901/76,911, 78,650/78,671, 77,910/77,934, 76,905/76,924, 0 mismatches each (controls and ancient identical). The unmatched ~9-24 sites are ARG sites absent from the VCFs; I did not investigate them (unverified; probably multi-hit sites that the CLI skips). Position shift 0 is right.
* Panel-fixed: `prepare-fixed` found ~17-18k fixed sites per rep, of which only the fixed-derived ones are kept (973/775/893/866/820 for reps 1..5), added with age +inf.

## Results (50 ancient samples, mean over 5 draws; bias = inferred - true, generations)
`summary_by_a0.tsv`, per-sample `per_sample_results.tsv`, site selection `site_selection.tsv`, figure `inferred_vs_true.png/.svg`.

| source | A0 | n | mean bias | median bias | RMSE | r | slope | 0-1500 | 1500-3000 | 3000-4500 | 4500-6000 | targets flagged (any draw) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ARG | 10,000 | 50 | -714 | -806 | 911 | 0.948 | 0.806 | -207 | -616 | -991 | -1077 | 4 |
| ARG | 20,000 | 50 | -611 | -675 | 949 | 0.904 | 0.830 | -55 | -686 | -903 | -807 | 7 |
| ARG | 50,000 | 50 | +4025 | +3147 | 5192 | 0.917 | 2.576 | +1037 | +2351 | +5137 | +8253 | 8 |
| true tree (check a) | 10,000 | 50 | +198 | +101 | 652 | 0.958 | 1.139 | -31 | +79 | +182 | +662 | - |
| true tree | 20,000 | 50 | +371 | +264 | 863 | 0.952 | 1.233 | +6 | +138 | +359 | +1140 | - |
| true tree | 50,000 | 50 | +5876 | +5168 | 7424 | 0.905 | 3.212 | +1640 | +3542 | +7343 | +11985 | - |

Sites selected (mean over draws, fraction of usable ARG sites): A0 10,000: ~28.2k (36%); 20,000: ~16.7k (21%); 50,000: ~3.1k (4%) (per rep in `site_selection.tsv`).
Flagged targets are mostly `H<=0` at the youngest ages; their T_hat (truncated 0) is INCLUDED in the statistics above (flag columns in the TSV).
Check (b), controls as pseudo-targets (leave-one-out against the other 9 controls; `controls_pseudotargets.tsv`; 250 control x draw values per A0 for the ARG, 50 for the true tree): ARG median T_hat 44/0/20 and mean 215/313/912 at A0 10k/20k/50k; true tree median 30/0/0, mean 302/432/1185. Near zero at 10k/20k (positive mean is expected from truncation at 0), but ~1000 at 50k.

## Interpretation (hypotheses, not tested)
* A0 = 10k-20k: with SINGER+POLEGON ARGs the estimator tracks true age (r 0.90-0.95) with a modest underestimate (~-600 to -700, slope ~0.8); the true tree gives a smaller, positive bias (+200 to +370), so I attribute the negative sign to ARG inference error (unverified: the A0/age definition interacts with inferred ages). Substantial across-draw spread for single samples.
* A0 = 50k: strongly biased high even with the TRUE tree, so this is the estimator at only ~3k selected sites (near the root, TMRCA ~77k), not ARG error. Why (few sites, breakdown of exp(integral c) assumption for mutations near the root) is not tested.
* Caveat: this is a constant-Ne simulation with a clean ARG panel; the real-data null result is not explained by this. Differences vs the real data: no variable Ne, no masks/missing data, only ~76k ARG sites and 20 Mb, REF known ancestral, samples haploid.
* No bugs found in pipeline/nonsharing_ratio*.py during this run. Observation only: the CLI does not print the allele-match counts to stdout (they are in the NPZ `meta`).

## Resources (sacct, observed)
| job | what | elapsed | MaxRSS |
|---|---|---|---|
| 40594477 / 40594556 | simulation rep 1 / reps 2-5 | 0:18 / 1:11 | 192 MB / 231 MB |
| 40594808-12 | SINGER+POLEGON per rep (12 CPU, 24G req.) | 15:06, 18:18, 17:28, 15:22, 14:58 | 5.7-6.2 GB (TotalCPU 1.7-2.0 h each, ~9.1 CPU-h total) |
| 40594830 | prepare-fixed + true-tree runs | 1:53 | 175 MB |
| 40595004/10/11/05/06 | run-draw x15 + merge x3 per rep (1 CPU) | 2:18, 2:32, 2:17, 2:19, 2:18 | 175-184 MB (jobs 40595004-06; 10/11 MaxRSS not recorded in my query) |
| 40595019 | analysis | 0:21 | 99 MB |
All COMPLETED; outputs were opened and checked.

## Layout
`scripts/` (01-06), `sim/rep*/` (VCFs, truth, true_panel.trees, sim_info.json), `singer_in/`, `singer_out/` (large; only `sim*/sim*/trees/*.tsz` needed), `draws/`, `merged/`, `logs/`, `singer_config_rep*.yaml`, `chr1.fai`.
