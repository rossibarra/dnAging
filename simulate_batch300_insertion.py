#!/usr/bin/env python3
"""Regenerate the deterministic 300-simulation benchmark for insertion inference."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import msprime
import numpy as np
import tskit

MASTER = 20260907


def parameters(i: int) -> dict:
    rng = np.random.default_rng(MASTER + i)
    mu = float(10 ** rng.uniform(-9, -8))
    return {
        "replicate": i,
        "seed": int(rng.integers(1, 2**31 - 1)),
        "Ne": int(rng.integers(10_000, 100_001)),
        "mutation_rate": mu,
        "recombination_rate": mu * float(rng.choice([0.5, 1.0, 2.0])),
        "n_modern_haploid": int(rng.integers(10, 41)),
        "true_ancient_age": float(rng.uniform(100, 10_000)),
        "sequence_length": 10_000_000.0,
    }


def main(args):
    p = parameters(args.replicate)
    out = args.output_root / f"replicate_{args.replicate:03d}"
    if out.exists():
        raise SystemExit(f"refusing to overwrite {out}")
    tmp = out.with_name(out.name + ".tmp")
    if tmp.exists():
        raise SystemExit(f"refusing to overwrite {tmp}")
    tmp.mkdir(parents=True)

    full = msprime.sim_ancestry(
        samples=[msprime.SampleSet(p["n_modern_haploid"], time=0, ploidy=1),
                 msprime.SampleSet(1, time=p["true_ancient_age"], ploidy=1)],
        population_size=p["Ne"], sequence_length=p["sequence_length"],
        recombination_rate=p["recombination_rate"], ploidy=2,
        random_seed=p["seed"])
    full = msprime.sim_mutations(full, rate=p["mutation_rate"],
                                 random_seed=p["seed"] + 1)
    modern = np.array([u for u in full.samples() if full.node(u).time == 0], dtype=np.int32)
    ancient = [u for u in full.samples() if full.node(u).time == p["true_ancient_age"]][0]
    modern_genotypes = full.genotype_matrix(samples=modern)
    keep = [s.id for s in full.sites()
            if len(s.mutations) == 1 and s.position >= 1
            and 0 < np.count_nonzero(modern_genotypes[s.id])
            < p["n_modern_haploid"]]
    truth = full.delete_sites(np.setdiff1d(np.arange(full.num_sites), keep))
    modern_arg = truth.simplify(modern, filter_sites=False)
    modern_arg.dump(tmp / "known_modern_arg.trees")
    truth.dump(tmp / "truth_with_ancient.trees")
    individuals = [truth.node(u).individual for u in truth.samples()]
    names = [f"modern_{j:02d}" for j in range(1, p["n_modern_haploid"] + 1)] + ["ancient_01"]
    with gzip.open(tmp / "all_samples.vcf.gz", "wt") as handle:
        truth.write_vcf(handle, contig_id="1", individuals=individuals,
                        individual_names=names,
                        position_transform=lambda x: np.fmax(1, x))
    p.update({"ancestry_seed": p["seed"], "mutation_seed": p["seed"] + 1,
              "modern_polymorphic_sites": modern_arg.num_sites,
              "num_trees": modern_arg.num_trees,
              "benchmark_master_seed": MASTER})
    (tmp / "metadata.json").write_text(json.dumps(p, indent=2) + "\n")
    tmp.replace(out)
    print(json.dumps(p), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    main(parser.parse_args())
