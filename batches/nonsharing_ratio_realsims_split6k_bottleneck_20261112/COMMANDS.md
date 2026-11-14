# Commands (repo root /quobyte/jrigrp/beil/logan_collab/dnAging; B=batches/nonsharing_ratio_realsims_split6k_bottleneck_20261112; hpc_run=~/.claude/bin/hpc_run)

    # 1. simulate 5 reps (split 6000 + bottleneck)                      job 40619415
    HPC_MEM=8G HPC_TIME=00:40:00 hpc_run "for r in 1 2 3 4 5; do python $B/scripts/01_simulate.py \$r $B/sim; done"
    # 2. SINGER inputs: singer_in/sim$r/{sim$r.vcf.gz, sim$r.hapmap}; singer_config_rep$r.yaml
    # 3. SINGER + POLEGON, one job per rep                              jobs 40619425-29
    for r in 1 2 3 4 5; do sbatch --export=ALL,REP=$r $B/scripts/02_singer_polegon.sbatch; done
    # 4. prepare-fixed + true-tree runs (c(t) from tree and from demes) job 40619430
    HPC_MEM=16G HPC_TIME=01:30:00 hpc_run "bash $B/scripts/03_prepare_fixed_and_true.sh"
    # 5. run-draw x5 draws x3 A0 (ARG c(t) and demes c(t)) + merge, one job per rep
    for r in 1 2 3 4 5; do sbatch --parsable --export=ALL,REP=$r $B/scripts/04_run_draws.sbatch; done   # 40627391-95
    # 6. analysis + offset ratio
    HPC_MEM=8G HPC_TIME=00:30:00 hpc_run "python -I $B/scripts/05_analyze.py && python -I $B/scripts/06_offset.py"   # 40627525
