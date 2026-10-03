#!/bin/bash -l
#SBATCH --job-name=nine_argtest
#SBATCH --account=jrigrp
#SBATCH --partition=low
#SBATCH --qos=jrigrp-low-qos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=nine-argtest-%j.out

set -euo pipefail
export PROJECT=/quobyte/jrigrp/beil/logan_collab/dnAging
export BATCH_ROOT=/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nine_demographies_20260916
export PYTHONPATH="$BATCH_ROOT/code:$PROJECT:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load conda
conda activate dnaging
RUN_INDEX="${RUN_INDEX:?}"
RUN_LINE=$(awk -F '\t' -v n="$((RUN_INDEX+2))" 'NR==n' "$BATCH_ROOT/runs.tsv")
IFS=$'\t' read -r INDEX RUN_NAME MODEL REPLICATE RUN_ROOT BASE_SEED SINGER_SEED <<< "$RUN_LINE"
export RUN_ROOT
test "$INDEX" = "$RUN_INDEX"
cd "$RUN_ROOT"
ROOT="$RUN_ROOT"
W="$ROOT/argtest"
SRC="$ROOT/singer/output"
mkdir -p "$W"
exec 9>"$W/.resume.guard"
flock -n 9 || { echo "Another resume process owns this workflow" >&2; exit 1; }
# Preserve existing biological inputs and configuration when resuming.
if [[ ! -s "$W/config.yaml" ]]; then
mkdir -p "$W/input"
ln -sfn /home/jri/src/argtest/scripts "$W/scripts"

# Flat per-chromosome layout: argtest's infer_mu_path searches for
# <parent_dir_name>.mut_rate.p beside the trees, so chrN/ must hold BOTH the
# .tsz files and chrN.mut_rate.p. A tree_subdir layout would hide the SINGER
# rate map and silently demote step 3 to the scalar mutation_rate fallback.
: > "$W/input/all.hapmap.tsv"
printf 'Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n' > "$W/input/all.hapmap.tsv"
: > "$W/input/sim.fai"
for c in chr1 chr2 chr3; do
    mkdir -p "$W/input/$c"
    ln -sf "$SRC/$c"/trees/*.tsz "$W/input/$c/"
    ln -sf "$SRC/$c/$c.mut_rate.p" "$W/input/$c/$c.mut_rate.p"
    LEN=$(python3 -c "import json,sys;print(int(json.load(open('$ROOT/singer/input/$c.prepare_summary.json'))['sequence_length']))")
    printf '%s\t%s\t0\t60\t61\n' "$c" "$LEN" >> "$W/input/sim.fai"
    # r = 1e-8 per bp per generation = 1 cM/Mb, uniform, past the last variant.
    printf '%s\t1\t1.0\t0.0\n' "$c" >> "$W/input/all.hapmap.tsv"
    printf '%s\t%s\t1.0\t%s\n' "$c" "$LEN" \
        "$(python3 -c "print(f'{$LEN/1e6:.6f}')")" >> "$W/input/all.hapmap.tsv"
done

python - <<PY
import yaml
from pathlib import Path
base = yaml.safe_load(open("/home/jri/src/argtest/config/snakemake.yaml"))
w = Path("$W")
cfg = dict(base)                                  # every shipped default kept
cfg["root_dir"] = str(w / "input")
cfg["hapmap"] = str(w / "input/all.hapmap.tsv")
cfg["fai"] = str(w / "input/sim.fai")
cfg["tree_pattern"] = "*.tsz"
cfg["out_dir"] = str(w / "out")
cfg["rec_fraction"] = 0.0                         # requested deviation 1
cfg["validation_first_chrom_only"] = False        # requested deviation 2
(w / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
print("config written")
PY
fi

conda activate argtest
if [[ -d "$W/.snakemake/locks" ]]; then
 cp -an "$W/.snakemake/locks" "$W/.snakemake/locks.pre_unlock_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}.bak"
fi
snakemake --snakefile /home/jri/src/argtest/Snakefile \
    --directory "$W" --configfile "$W/config.yaml" --unlock
snakemake \
    --snakefile /home/jri/src/argtest/Snakefile \
    --directory "$W" \
    --configfile "$W/config.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-32}" \
    --rerun-incomplete
