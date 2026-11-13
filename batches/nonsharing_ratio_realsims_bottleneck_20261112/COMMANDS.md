# Commands (repo root /quobyte/jrigrp/beil/logan_collab/dnAging; B=batches/nonsharing_ratio_realsims_bottleneck_20261112; hpc_run=~/.claude/bin/hpc_run)

    # 1. simulate 5 reps with the bottleneck (demography_bottleneck.yaml)
    HPC_MEM=8G HPC_TIME=00:40:00 hpc_run "for r in 1 2 3 4 5; do python $B/scripts/01_simulate.py \$r $B/sim; done"
    # 2. SINGER inputs: singer_in/sim$r/{sim$r.vcf.gz (= sim/rep$r/panel_segregating.vcf.gz), sim$r.hapmap (copied from template)}; singer_config_rep$r.yaml = template configs with paths changed
    # 3. SINGER + POLEGON, one job per rep (low, 12 CPU, 24G)
    for r in 1 2 3 4 5; do sbatch --export=ALL,REP=$r $B/scripts/02_singer_polegon.sbatch; done      # 40598362..66
    # 4. prepare-fixed + true-tree runs (c(t) from tree and from demes)
    HPC_MEM=16G HPC_TIME=01:30:00 hpc_run "bash $B/scripts/03_prepare_fixed_and_true.sh"
    # 5. run-draw x5 draws x3 A0, ARG c(t) and demes c(t), + merge, one job per rep
    sbatch --export=ALL,REP=$r $B/scripts/04_run_draws.sbatch     # 40598778 (1) 40598781 (2) 40598782 (3) 40598798 (4) 40598779 (5)
    # 6. analysis, tables, figures
    HPC_MEM=8G HPC_TIME=00:30:00 hpc_run "python $B/scripts/05_analyze.py"                            # 40598805
    # 7. cleanup
    rm -rf singer_out/sim$r/sim$r/{chunks,plots,sim$r.vcf}   (r=1..5, after the 25 .tsz were present)

Extra flags vs the template for the demes runs: `--demography $B/demography_bottleneck.yaml --ploidy 1` (file in haploid counts).
