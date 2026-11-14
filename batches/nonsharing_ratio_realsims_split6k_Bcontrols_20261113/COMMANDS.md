# Commands (repo root /quobyte/jrigrp/beil/logan_collab/dnAging; B=batches/nonsharing_ratio_realsims_split6k_Bcontrols_20261113; hpc_run=~/.claude/bin/hpc_run)

    # 1. simulate 5 reps: 26 panel (A), 10 controls A, 10 controls B, 10 ancient (B)
    HPC_MEM=8G HPC_TIME=00:40:00 hpc_run "for r in 1 2 3 4 5; do python $B/scripts/01_simulate.py \$r $B/sim; done"
    # 2. SINGER inputs: singer_in/sim$r/{sim$r.vcf.gz (= sim/rep$r/panel_segregating.vcf.gz), sim$r.hapmap (copied from the earlier split batch)}; singer_config_rep$r.yaml = earlier configs with paths changed
    # 3. SINGER + POLEGON, one job per rep (low, 12 CPU, 24G)                     jobs 40627402-06
    for r in 1 2 3 4 5; do sbatch --export=ALL,REP=$r $B/scripts/02_singer_polegon.sbatch; done
    # 4. prepare-fixed + true-tree runs, both control arms (true_*, trueB_*, truedemes_*, truedemesB_*)
    HPC_MEM=16G HPC_TIME=01:30:00 hpc_run "bash $B/scripts/03_prepare_fixed_and_true.sh"
    # 5. run-draw x5 draws x3 A0 x2 arms (+ demes c(t)) + merge, one job per rep, after the SINGER job of that rep
    sbatch --parsable --dependency=afterok:$SINGER_JOB --export=ALL,REP=$r $B/scripts/04_run_draws.sbatch     # 40627541-45 (reps 1-5)
    # 6. analysis (both arms, offset ratio, cross pseudo-targets, figures)
    HPC_MEM=8G HPC_TIME=00:30:00 hpc_run "python -I $B/scripts/05_analyze.py"                                  # 40627663
    # 7. NOT done: SINGER intermediates singer_out/sim*/sim*/{chunks,plots,sim*.vcf} (~19 GB per batch) are still present; they are gitignored and can be deleted once the .tsz files are no longer needed.
Arm A = controls.vcf.gz, arm B = controlsB.vcf.gz; the arms differ only in --controls-vcf.
