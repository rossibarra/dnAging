#!/bin/bash
# from repo root; one hpc_run (srun) per draw, each with its own log file
M=batches/nonsharing_ratio_heldout_20261111/more_draws_check
for d in 82 89 96 103 110 117 124 131 138; do
  A0=20000 EXTRA="--age-mode child" HPC_CPUS=2 HPC_MEM=32G HPC_TIME=01:30:00 \
    ~/.claude/bin/hpc_run "bash $M/scripts/run_draw_deid.sh $d" > $M/logs/draw_$d.log 2>&1 &
done
wait
