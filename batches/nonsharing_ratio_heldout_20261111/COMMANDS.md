# Exact commands (run from this directory unless noted)

    bash -lc 'module load bcftools/1.19; python3 scripts/01_select_heldout.py inputs'     # seed 20261111, k=6
    sbatch scripts/02_prep_inputs.sbatch          # job 40590168
    sbatch --dependency=afterok:40590168 scripts/03_singer_polegon.sbatch   # job 40590169 (low, 48 CPU, 64G, 24h, --requeue; resumable)
    # verification, from the repo root (hpc_run activates ./environment.yml):
    HPC_MEM=24G HPC_TIME=01:00:00 ~/.claude/bin/hpc_run 'bash batches/nonsharing_ratio_heldout_20261111/scripts/04_verify_run_draw.sh 75'
    HPC_MEM=24G HPC_TIME=00:30:00 ~/.claude/bin/hpc_run 'python batches/nonsharing_ratio_heldout_20261111/scripts/06_inspect_and_deidentify.py 75'

04_verify_run_draw.sh expands to (ancient chr10 export written only to node-local scratch, deleted on exit):

    python pipeline/nonsharing_ratio_real_data.py run-draw \
      --tree singer_out/chr.10.combined.snp.te.sorted/trees/chr.10.combined.snp.te.sorted.75.tsz --draw-id 75 \
      --controls-vcf controls/heldout_controls_chr10.vcf.gz --ancient-vcf <node-local chr10 export of the real ancient VCF> \
      --panel-names panel_names.txt --chrom-fai chr10.fai --a0 12000 --output verify/run_draw_75.npz
