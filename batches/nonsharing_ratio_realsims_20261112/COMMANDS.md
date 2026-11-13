# Commands (run from the repo root unless noted; `hpc_run` = ~/.claude/bin/hpc_run, activates ./environment.yml)

    # 1. simulate (job 40594403 failed: tskit mutation-time error in the true-tree export; fixed, rerun 40594477 for rep 1)
    HPC_MEM=8G HPC_TIME=00:20:00 hpc_run 'python batches/nonsharing_ratio_realsims_20261112/scripts/01_simulate.py 1 batches/nonsharing_ratio_realsims_20261112/sim'   # job 40594477
    HPC_MEM=8G HPC_TIME=00:40:00 hpc_run 'for r in 2 3 4 5; do python .../01_simulate.py $r .../sim; done'                                                          # job 40594556
    # 2. SINGER inputs (done interactively, light): per rep, singer_in/sim$r/{sim$r.vcf.gz (=sim/rep$r/panel_segregating.vcf.gz), sim$r.hapmap (1 cM/Mb, 20 Mb)},
    #    singer_config_rep$r.yaml = heldout_20261111 config with mu 1.25e-8, rho 1e-8, seed 56789+r, max-masked-samples 0, chromosomes [sim$r]
    # 3. SINGER + POLEGON, one job per replicate (low, 12 CPU, 24G, 1 day)
    for r in 1 2 3 4 5; do sbatch --export=ALL,REP=$r batches/nonsharing_ratio_realsims_20261112/scripts/02_singer_polegon.sbatch; done   # jobs 40594808..40594812 (rep1..5)
    # 4. prepare-fixed (pre-dropping panel VCF) + true-tree sanity check (a)
    HPC_MEM=16G HPC_TIME=01:00:00 hpc_run 'bash batches/nonsharing_ratio_realsims_20261112/scripts/03_prepare_fixed_and_true.sh'   # job 40594830
    # 5. run-draw (5 draws x 3 A0) + merge per rep/A0, one job per rep (low, 1 CPU, 24G)
    sbatch --export=ALL,REP=$r batches/nonsharing_ratio_realsims_20261112/scripts/04_run_draws.sbatch   # 40595004 (rep1) 40595010 (rep2) 40595011 (rep3) 40595005 (rep4) 40595006 (rep5)
    # 6. analysis / figure / pseudo-target controls
    HPC_MEM=8G HPC_TIME=00:20:00 hpc_run 'python batches/nonsharing_ratio_realsims_20261112/scripts/05_analyze.py'   # job 40595019
    # inspection of a draw: scripts/06_inspect_tree.py   # job 40595007

run-draw expands to (per rep r, draw d in 75 93 112 131 149, A0 in 10000 20000 50000):

    python pipeline/nonsharing_ratio_real_data.py run-draw --tree singer_out/sim$r/sim$r/trees/sim$r.$d.tsz --draw-id $d \
      --controls-vcf sim/rep$r/controls.vcf.gz --controls-mode sample --ancient-vcf sim/rep$r/ancient.vcf.gz \
      --panel-names sim/rep$r/panel_names.txt --chrom-fai chr1.fai --a0 $A0 --age-mode child \
      --panel-fixed draws/rep$r/panel_fixed.npz --fixed-polarity ref --output draws/rep$r/draw${d}_a0$A0.npz
