#!/bin/bash -l
#SBATCH --job-name=3epArgtest
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=4-00:00:00

set -euo pipefail
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ROOT="$PROJECT/msprime_three_epoch"
W="$ROOT/argtest"
SRC="$ROOT/singer/output"
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

module load conda
conda activate dnaging
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

conda activate argtest
snakemake \
    --snakefile /home/jri/src/argtest/Snakefile \
    --directory "$W" \
    --configfile "$W/config.yaml" \
    --cores "${SLURM_CPUS_PER_TASK:-32}" \
    --rerun-incomplete
