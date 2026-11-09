#!/bin/bash -l
# Build a single shareable archive of constant-size simulations and POLEGON artifacts.
#SBATCH --job-name=packPolegonConst
#SBATCH --account=jrigrp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=packPolegonConst-%j.out
#SBATCH --error=packPolegonConst-%j.err
set -euo pipefail

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
ARCHIVE_NAME="polegon_constant_size_bundle_20261109"
FINAL="$PROJECT/sharing/$ARCHIVE_NAME.tgz"
STAGE="$TMPDIR/$ARCHIVE_NAME"
TEMP_ARCHIVE="$TMPDIR/$ARCHIVE_NAME.tgz"

test ! -e "$FINAL"
mkdir -p "$STAGE"/{data,scripts,analysis_outputs}

for run in constant_r1 constant_r2 constant_r3; do
  SOURCE="$PROJECT/batches/nine_demographies_20260916/runs/$run"
  DEST="$STAGE/data/$run"
  mkdir -p "$DEST"/{simulations,inferred_trees,combined_arg,draw50_chunks}
  cp -p "$SOURCE/demography.tsv" "$DEST/"
  for replicate in 001 002 003; do
    mkdir -p "$DEST/simulations/replicate_$replicate"
    for file in \
      truth_with_ancient.trees known_modern_arg.trees all_samples.vcf.gz \
      mutation_truth.tsv.gz metadata.json constant_ne_epochs.tsv; do
      cp -p "$SOURCE/simulations/replicate_$replicate/$file" \
        "$DEST/simulations/replicate_$replicate/"
    done
  done
  for chromosome in 1 2 3; do
    CHUNK_SOURCE="$SOURCE/singer/output/chr$chromosome/chunks"
    CHUNK_DEST="$DEST/draw50_chunks/chr$chromosome"
    mkdir -p "$CHUNK_DEST"
    mkdir -p "$DEST/inferred_trees/chr$chromosome"
    cp -p "$SOURCE/singer/output/chr$chromosome/trees/"*.tsz \
      "$DEST/inferred_trees/chr$chromosome/"
    cp -p "$CHUNK_SOURCE"/*.yaml \
      "$CHUNK_SOURCE"/*_nodes_50.txt \
      "$CHUNK_SOURCE"/*_branches_50.txt \
      "$CHUNK_SOURCE"/*_muts_50.txt \
      "$CHUNK_SOURCE"/*_polegon_50.txt \
      "$CHUNK_SOURCE"/*.50.log \
      "$CHUNK_DEST/"
  done
  cp -p "$SOURCE/argtest/out/combined/"*.tsz \
    "$DEST/combined_arg/"

  mkdir -p "$STAGE/analysis_outputs/$run"
  cp -a \
    "$PROJECT/batches/joint_recovery_matrix_20261029/oracle/baseline/$run" \
    "$STAGE/analysis_outputs/$run/oracle"
  cp -a \
    "$PROJECT/batches/joint_recovery_matrix_20261029/inferred/$run" \
    "$STAGE/analysis_outputs/$run/inferred"
done

for file in \
  environment.yml \
  validation/simulate_variable_ne_error.py \
  validation/chunked_true_topology_polegon.py \
  validation/true_topology_polegon.py \
  validation/per_chromosome_calibration.py \
  validation/joint_scale_epsilon_recovery.py \
  validation/radiocarbon_oracle_stage1.py \
  pipeline/insertion_likelihood.py \
  pipeline/insertion_real_data.py \
  pipeline/radiocarbon_calibration.py \
  slurm/run_chunked_true_topology_polegon.sbatch \
  slurm/run_chunked_true_topology_calibration.sbatch \
  slurm/package_polegon_constant_bundle.sbatch; do
  mkdir -p "$STAGE/scripts/$(dirname "$file")"
  cp -p "$PROJECT/$file" "$STAGE/scripts/$file"
done

cat > "$STAGE/README.md" <<'EOF'
# Constant-size ARG and POLEGON diagnostic bundle

This archive contains only constant-population-size simulations and their
associated SINGER/POLEGON posterior tree outputs. It was assembled on
2026-11-09 from the dnAging project to support investigation of ARG node-time
calibration.

## Simulation design

- Population size: constant diploid Ne = 20,000.
- Three run sets: `constant_r1`, `constant_r2`, and `constant_r3`.
- Three independent 50 Mb chromosome replicates per run set.
- Mutation and recombination rates: 1e-8 per base per generation.
- Each truth ARG has 100 modern haploid samples and 21 serial samples at ages
  0, 500, ..., 10,000 generations. The age-zero serial sample is a control, so
  20 samples are genuinely ancient.
- Mutations were simulated with msprime's infinite-sites nucleotide model and
  retained only when polymorphic in the modern panel.

## Layout

`data/<run>/simulations/replicate_00N/`

- `truth_with_ancient.trees`: complete simulated ARG with modern and serial
  samples. This is the ground truth.
- `known_modern_arg.trees`: the same simulation simplified to the 100 modern
  samples, retaining true node and mutation times. This is the clean topology
  to re-date with POLEGON.
- `all_samples.vcf.gz`: modern and serial-sample genotypes.
- `mutation_truth.tsv.gz`: true mutation times, nodes, frequencies, and ancient
  alleles.
- `metadata.json`: sample-node IDs, true ages, seeds, versions, and rates.
- `constant_ne_epochs.tsv`: demographic table used by the age likelihood.

`data/<run>/draw50_chunks/chrN/`

- Exact production SINGER chunk YAML files.
- Raw SINGER draw-50 node, branch, and mutation text tables.
- POLEGON draw-50 output node times (`*_polegon_50.txt`).
- POLEGON draw-50 logs, including the exact command and parameters.

These chunk files describe inferred SINGER topologies. For a clean fixed-
topology POLEGON experiment, use `known_modern_arg.trees` with
`scripts/validation/chunked_true_topology_polegon.py` and the corresponding
chunk YAML directory.

`data/<run>/inferred_trees/`

- All final per-chromosome, POLEGON-dated SINGER posterior tree sequences.
- Files are grouped by chromosome and include every saved posterior draw.

`data/<run>/combined_arg/`

- All saved three-chromosome posterior tree sequences produced by ARGtest.

`analysis_outputs/<run>/`

- `oracle/`: likelihood and cross-validation results using true modern ARGs.
- `inferred/`: likelihood and cross-validation results using inferred,
  POLEGON-dated ARGs.

`scripts/`

- Simulation, POLEGON adapter, likelihood, calibration, and Slurm scripts.
- `environment.yml` records the project software environment.

## Important interpretation

POLEGON should be run on `known_modern_arg.trees`, not on
`truth_with_ancient.trees`. The latter retains ancient nodes only to provide
complete simulation truth. The node times in `known_modern_arg.trees` are known
truth and can be compared directly with POLEGON's estimates while topology is
held fixed.

The bundled POLEGON executable is not redistributed. The original runs used
POLEGON 0.1.3-alpha; exact invocations are present in each `*.50.log` file.

## Integrity

`SHA256SUMS` contains checksums for every other file in this archive.
EOF

(
  cd "$STAGE"
  find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
)
tar -C "$TMPDIR" -czf "$TEMP_ARCHIVE" "$ARCHIVE_NAME"
gzip -t "$TEMP_ARCHIVE"

mkdir -p "$PROJECT/sharing"
COPY_TARGET="$FINAL.tmp.$SLURM_JOB_ID"
cp -p "$TEMP_ARCHIVE" "$COPY_TARGET"
mv "$COPY_TARGET" "$FINAL"
sha256sum "$FINAL" > "$FINAL.sha256"
echo "Created $FINAL"
ls -lh "$FINAL" "$FINAL.sha256"
