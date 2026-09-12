#!/usr/bin/env python3
"""Variable-Ne msprime replicates with genotyping error on the ancient calls.

Every simulation control so far has been error-free and constant-Ne, because the
working documents park real-data concerns until the estimator is unbiased on
perfect data.  Both of those conditions are now the thing being tested: a
ten-epoch piecewise-constant demography exercises the diffusion-time change of
clock (MATH.md eq. 5), and a symmetric per-allele error on the ancient calls
exercises the epsilon term (MATH.md eq. 2, r = eps + (1-2 eps) X) that has never
had a nonzero value in a simulation.

Design choices, none of which the request fixed:

* **Epoch boundaries are log-spaced.**  Recent demography is resolved finely and
  ancient demography coarsely, which is how inferred Ne(t) is actually reported
  (this project's real analysis uses ~50 log-spaced windows from ARGtest).  The
  tenth epoch is ancestral and extends to infinity.
* **Sizes are log-uniform on [5e3, 1.5e5] and drawn independently per epoch.**
  Log-uniform because Ne is a scale parameter spanning 30x here; independent
  because the resulting sawtooth is a *harder* test of the tau integral than a
  smooth trajectory would be.  `--ne-prior uniform` and `--smooth` are available.
* **Error is applied only to the ancient calls**, which is what eps means in
  MATH.md section 2: it describes error in each ancient-VCF allele call and
  nothing else.  The modern panel is the ARG panel and stays exact.
* **Truth is retained alongside the corrupted calls.**  `mutation_truth.tsv.gz`
  carries both the true and the observed ancient genotype per site, so the
  effect of the error is separable after the fact rather than baked in.

One ARG serves every ancient sample, so asking for ten ancient haploids at ten
random ages costs essentially nothing beyond one ancient haploid, and mirrors
the real cohort where many samples share one panel.

Writes `constant_ne_epochs.tsv` per replicate in the series/time_left/time_right
format `precompute_freq_trajectory_moments.py` consumes, because each replicate
has its own demography and therefore needs its own frequency table.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import msprime
import numpy as np
import tskit

SEED_STREAMS = ("demography", "age", "ancestry", "mutation", "error")


def seeds_for_replicate(replicate: int, base_seed: int) -> dict:
    """One independent, reproducible seed per stream.

    Separate streams matter here: re-drawing the error without re-drawing the
    ancestry is what makes an error-free control on identical genealogies
    possible.
    """
    return {name: base_seed + (i + 1) * 100_000 + replicate
            for i, name in enumerate(SEED_STREAMS)}


def draw_demography(rng, args):
    """Ten piecewise-constant epochs; returns (boundaries, sizes)."""
    boundaries = np.geomspace(args.epoch_min_time, args.epoch_max_time,
                              args.n_epochs - 1)
    if args.ne_prior == "loguniform":
        sizes = np.exp(rng.uniform(np.log(args.ne_min), np.log(args.ne_max),
                                   args.n_epochs))
    else:
        sizes = rng.uniform(args.ne_min, args.ne_max, args.n_epochs)
    if args.smooth:
        # A random walk in log-Ne instead of independent draws, for a trajectory
        # that looks like an inferred one rather than a sawtooth.
        step = np.log(args.ne_max / args.ne_min) / 4.0
        walk = np.cumsum(rng.normal(0.0, step, args.n_epochs))
        sizes = np.clip(np.exp(np.log(np.sqrt(args.ne_min * args.ne_max)) + walk),
                        args.ne_min, args.ne_max)
    return boundaries, sizes


def build_demography(boundaries, sizes):
    demography = msprime.Demography()
    demography.add_population(name="pop", initial_size=float(sizes[0]))
    for time, size in zip(boundaries, sizes[1:]):
        demography.add_population_parameters_change(
            time=float(time), initial_size=float(size), population="pop")
    return demography


def write_ne_table(path: Path, boundaries, sizes, horizon):
    """The step function in the form precompute_freq_trajectory_moments reads."""
    edges = [0.0, *[float(b) for b in boundaries], float(horizon)]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["series", "time_left", "time_right",
                         "effective_population_size"])
        for left, right, size in zip(edges[:-1], edges[1:], sizes):
            writer.writerow(["variable", f"{left:.6f}", f"{right:.6f}",
                             f"{float(size):.6f}"])


def simulate(args):
    seeds = seeds_for_replicate(args.replicate, args.base_seed)
    demography_rng = np.random.default_rng(seeds["demography"])
    age_rng = np.random.default_rng(seeds["age"])
    error_rng = np.random.default_rng(seeds["error"])

    boundaries, sizes = draw_demography(demography_rng, args)
    ancient_ages = np.sort(age_rng.uniform(args.age_min, args.age_max,
                                           args.n_ancient))

    out = args.output_root / f"replicate_{args.replicate:03d}"
    if out.exists():
        raise SystemExit(f"refusing to overwrite existing replicate: {out}")
    tmp = out.with_name(out.name + f".tmp.{args.replicate}")
    tmp.mkdir(parents=True)

    # Each ancient sample is its own haploid individual at its own time, so a
    # single genealogy serves every age.
    sample_sets = [msprime.SampleSet(args.n_modern, time=0, ploidy=1)]
    sample_sets += [msprime.SampleSet(1, time=float(t), ploidy=1)
                    for t in ancient_ages]
    full = msprime.sim_ancestry(
        samples=sample_sets,
        demography=build_demography(boundaries, sizes),
        sequence_length=args.length,
        recombination_rate=args.recombination_rate,
        ploidy=2,
        model="hudson",
        random_seed=seeds["ancestry"],
    )
    full = msprime.sim_mutations(
        full, rate=args.mutation_rate,
        model=msprime.InfiniteSites(msprime.NUCLEOTIDES),
        discrete_genome=False, random_seed=seeds["mutation"],
    )
    if full.num_sites != full.num_mutations:
        raise ValueError("infinite-sites simulation produced a recurrent site")

    times = np.array([full.node(u).time for u in full.samples()])
    modern_nodes = np.asarray(full.samples())[times == 0].astype(np.int32)
    if len(modern_nodes) != args.n_modern:
        raise ValueError("unexpected modern sample-node count")
    # Match each ancient age to its node by time; ages are distinct with
    # probability one under a continuous draw, and duplicates are rejected.
    ancient_nodes = []
    for age in ancient_ages:
        hits = np.asarray(full.samples())[np.isclose(times, age)]
        if len(hits) != 1:
            raise ValueError(f"expected exactly one ancient node at time {age}")
        ancient_nodes.append(int(hits[0]))
    ancient_nodes = np.asarray(ancient_nodes, dtype=np.int32)

    # Ascertain on modern polymorphism, as the real pipeline does.
    keep = [v.site.id for v in full.variants(samples=modern_nodes)
            if 0 < int(np.count_nonzero(v.genotypes)) < args.n_modern]
    truth = full.delete_sites(np.setdiff1d(
        np.arange(full.num_sites, dtype=np.int64),
        np.asarray(keep, dtype=np.int64), assume_unique=True))
    modern_arg = truth.simplify(modern_nodes, filter_sites=False)
    if modern_arg.num_sites != truth.num_sites:
        raise ValueError("modern-only simplification changed retained site count")
    if np.any(tskit.is_unknown_time(modern_arg.tables.mutations.time)):
        raise ValueError("msprime returned unknown mutation times")

    modern_matrix = modern_arg.genotype_matrix()                 # (sites, modern)
    ancient_true = truth.genotype_matrix(samples=ancient_nodes)  # (sites, ancient)
    ancient_true = (ancient_true != 0).astype(np.int8)

    # Symmetric per-allele genotyping error: each observed ancient allele is
    # flipped independently with probability eps.  This is exactly the process
    # MATH.md eq. (2) models, so a run at eps=0.01 tests the term rather than the
    # estimator's tolerance of an unmodelled corruption.
    flips = error_rng.random(ancient_true.shape) < args.epsilon
    ancient_observed = np.where(flips, 1 - ancient_true, ancient_true)

    names = ([f"modern_{i:02d}" for i in range(1, args.n_modern + 1)]
             + [f"ancient_{i:02d}" for i in range(1, args.n_ancient + 1)])
    # Written by hand rather than through ts.write_vcf, because the ancient
    # columns must carry the CORRUPTED calls while the tree sequence keeps truth.
    with gzip.open(tmp / "all_samples.vcf.gz", "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write(f"##contig=<ID=1,length={int(args.length)}>\n")
        handle.write("#" + "\t".join(
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
             "FORMAT", *names]) + "\n")
        for k, site in enumerate(modern_arg.sites()):
            pos = max(1, int(site.position))
            gts = [str(int(g != 0)) for g in modern_matrix[k]]
            gts += [str(int(g)) for g in ancient_observed[k]]
            handle.write("\t".join(
                ["1", str(pos), str(site.id), site.ancestral_state,
                 site.mutations[0].derived_state, ".", "PASS", ".", "GT", *gts])
                + "\n")

    with gzip.open(tmp / "mutation_truth.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["site_id", "position", "mutation_id", "mutation_time",
                         "mutation_node", "modern_derived_count",
                         *[f"ancient_{i:02d}_true" for i in range(1, args.n_ancient + 1)],
                         *[f"ancient_{i:02d}_observed" for i in range(1, args.n_ancient + 1)]])
        for k, site in enumerate(modern_arg.sites()):
            mutation = site.mutations[0]
            writer.writerow([
                site.id, f"{site.position:.17g}", mutation.id,
                f"{mutation.time:.17g}", mutation.node,
                int(np.count_nonzero(modern_matrix[k])),
                *[int(v) for v in ancient_true[k]],
                *[int(v) for v in ancient_observed[k]]])

    truth.dump(tmp / "truth_with_ancient.trees")
    modern_arg.dump(tmp / "known_modern_arg.trees")
    write_ne_table(tmp / "constant_ne_epochs.tsv", boundaries, sizes,
                   args.ne_table_horizon)

    edges = [0.0, *[float(b) for b in boundaries], float("inf")]
    metadata = {
        "replicate": args.replicate,
        "seeds": seeds,
        "epochs": [{"index": i, "time_left": edges[i], "time_right": edges[i + 1],
                    "effective_population_size": float(sizes[i])}
                   for i in range(args.n_epochs)],
        "ne_prior": args.ne_prior, "ne_min": args.ne_min, "ne_max": args.ne_max,
        "smooth_ne": bool(args.smooth),
        "true_ancient_ages": [float(t) for t in ancient_ages],
        "n_modern_haploid": args.n_modern,
        "n_ancient_haploid": args.n_ancient,
        "epsilon_applied": args.epsilon,
        "ancient_alleles_flipped": int(flips.sum()),
        "ancient_alleles_total": int(flips.size),
        "realised_flip_rate": float(flips.mean()) if flips.size else 0.0,
        "sequence_length": args.length,
        "mutation_rate": args.mutation_rate,
        "recombination_rate": args.recombination_rate,
        "ancestry_model": "hudson",
        "mutation_model": "InfiniteSites(NUCLEOTIDES)",
        "discrete_genome": False,
        "sites_before_modern_ascertainment": full.num_sites,
        "modern_polymorphic_sites": int(modern_arg.num_sites),
        "num_trees": int(modern_arg.num_trees),
        "msprime_version": msprime.__version__,
        "tskit_version": tskit.__version__,
        "analysis_arg": "known_modern_arg.trees (ancient haplotypes removed)",
        "truth_tree_sequence": "truth_with_ancient.trees",
        "combined_vcf": "all_samples.vcf.gz (ancient columns carry ERROR)",
        "ne_table": "constant_ne_epochs.tsv",
    }
    (tmp / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    tmp.replace(out)
    print(json.dumps(metadata, indent=2), flush=True)


def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output-root", type=Path, required=True)
    p.add_argument("--replicate", type=int, required=True)
    p.add_argument("--base-seed", type=int, default=1730000000)
    p.add_argument("--n-epochs", type=int, default=10)
    p.add_argument("--ne-min", type=float, default=5_000.0)
    p.add_argument("--ne-max", type=float, default=150_000.0)
    p.add_argument("--ne-prior", choices=("loguniform", "uniform"),
                   default="loguniform")
    p.add_argument("--smooth", action="store_true",
                   help="random walk in log-Ne instead of independent epochs")
    p.add_argument("--epoch-min-time", type=float, default=100.0)
    p.add_argument("--epoch-max-time", type=float, default=100_000.0)
    p.add_argument("--ne-table-horizon", type=float, default=3_000_000.0,
                   help="right edge written for the ancestral epoch; must exceed "
                        "the table's --age-max")
    p.add_argument("--n-modern", type=int, default=26)
    p.add_argument("--n-ancient", type=int, default=10)
    p.add_argument("--age-min", type=float, default=0.0)
    p.add_argument("--age-max", type=float, default=10_000.0)
    p.add_argument("--length", type=float, default=10_000_000.0)
    p.add_argument("--mutation-rate", type=float, default=1e-8)
    p.add_argument("--recombination-rate", type=float, default=1e-8)
    p.add_argument("--epsilon", type=float, default=0.01,
                   help="symmetric per-allele error on ANCIENT calls only")
    return p


if __name__ == "__main__":
    simulate(build_parser().parse_args())
