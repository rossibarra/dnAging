#!/usr/bin/env python3
"""Simulate and analyse exact-mutation-time age calibration with msprime.

Each replicate contains 26 modern haploid samples and one ancient haploid at a
hidden uniform age. Sites are ascertained as polymorphic in the modern panel.
Inference uses the recorded msprime Mutation.time as a point mutation age and
calls the production phi_lookup and summarize functions.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

import msprime
import numpy as np
import tskit

from posterior_sample_age_infer import load_table, phi_lookup, summarize


def replicate_dir(root: Path, replicate: int) -> Path:
    return root / f"replicate_{replicate:03d}"


def seeds_for_replicate(replicate: int, base_seed: int) -> tuple[int, int, int]:
    # Separate, explicit streams. All values remain in msprime's valid uint32 range.
    return base_seed + replicate, base_seed + 100_000 + replicate, base_seed + 200_000 + replicate


def simulate(args: argparse.Namespace) -> None:
    replicate = args.replicate
    age_seed, ancestry_seed, mutation_seed = seeds_for_replicate(
        replicate, args.base_seed
    )
    age_rng = np.random.default_rng(age_seed)
    true_age = float(age_rng.uniform(args.age_min, args.age_max))

    out = replicate_dir(args.output_root, replicate)
    if out.exists():
        raise SystemExit(f"refusing to overwrite existing replicate directory: {out}")
    tmp = out.with_name(out.name + f".tmp.{replicate}")
    if tmp.exists():
        raise SystemExit(f"refusing to overwrite existing temporary directory: {tmp}")
    tmp.mkdir(parents=True)

    full = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(args.n_modern, time=0, ploidy=1),
            msprime.SampleSet(1, time=true_age, ploidy=1),
        ],
        population_size=args.ne,
        sequence_length=args.length,
        recombination_rate=args.recombination_rate,
        ploidy=2,
        model="hudson",
        random_seed=ancestry_seed,
    )
    full = msprime.sim_mutations(
        full,
        rate=args.mutation_rate,
        model=msprime.InfiniteSites(msprime.NUCLEOTIDES),
        discrete_genome=False,
        random_seed=mutation_seed,
    )
    if full.num_sites != full.num_mutations:
        raise ValueError("infinite-sites simulation produced a recurrent site")

    modern_nodes = np.asarray(
        [u for u in full.samples() if full.node(u).time == 0], dtype=np.int32
    )
    ancient_nodes = np.asarray(
        [u for u in full.samples() if full.node(u).time == true_age], dtype=np.int32
    )
    if len(modern_nodes) != args.n_modern or len(ancient_nodes) != 1:
        raise ValueError("unexpected modern/ancient sample-node counts")

    modern_polymorphic = []
    for variant in full.variants(samples=modern_nodes):
        derived_count = int(np.count_nonzero(variant.genotypes))
        if 0 < derived_count < args.n_modern:
            modern_polymorphic.append(variant.site.id)
    discard = np.setdiff1d(
        np.arange(full.num_sites, dtype=np.int64),
        np.asarray(modern_polymorphic, dtype=np.int64),
        assume_unique=True,
    )
    truth = full.delete_sites(discard)
    if truth.num_sites != truth.num_mutations:
        raise ValueError("ascertained truth is not infinite-sites")

    modern_arg = truth.simplify(modern_nodes, filter_sites=False)
    if modern_arg.num_samples != args.n_modern or modern_arg.num_sites != truth.num_sites:
        raise ValueError("modern-only simplification changed retained site count")
    if np.any(tskit.is_unknown_time(modern_arg.tables.mutations.time)):
        raise ValueError("msprime returned unknown mutation times")

    truth.dump(tmp / "truth_with_ancient.trees")
    modern_arg.dump(tmp / "known_modern_arg.trees")
    sample_individuals = [truth.node(u).individual for u in truth.samples()]
    if len(set(sample_individuals)) != args.n_modern + 1 or min(sample_individuals) < 0:
        raise ValueError("expected one haploid individual per sample node")
    names = [f"modern_{i:02d}" for i in range(1, args.n_modern + 1)] + ["ancient_01"]
    with gzip.open(tmp / "all_samples.vcf.gz", "wt") as handle:
        truth.write_vcf(
            handle, contig_id="1", individuals=sample_individuals,
            individual_names=names,
            # Continuous-genome sites can round to VCF POS=0. Keep coordinates
            # one-based-compliant while retaining site.id as the exact join key.
            position_transform=lambda positions: np.fmax(1, positions),
        )

    modern_matrix = modern_arg.genotype_matrix()
    ancient_matrix = truth.genotype_matrix(samples=ancient_nodes)[:, 0]
    with gzip.open(tmp / "mutation_truth.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "site_id", "position", "mutation_id", "mutation_time",
            "mutation_node", "modern_derived_count", "ancient_genotype",
        ])
        for site, modern_genotypes, ancient_genotype in zip(
            modern_arg.sites(), modern_matrix, ancient_matrix
        ):
            mutation = site.mutations[0]
            writer.writerow([
                site.id, f"{site.position:.17g}", mutation.id,
                f"{mutation.time:.17g}", mutation.node,
                int(np.count_nonzero(modern_genotypes)), int(ancient_genotype != 0),
            ])

    metadata = {
        "replicate": replicate,
        "true_ancient_age": true_age,
        "age_seed": age_seed,
        "ancestry_seed": ancestry_seed,
        "mutation_seed": mutation_seed,
        "Ne": args.ne,
        "n_modern_haploid": args.n_modern,
        "n_ancient_haploid": 1,
        "sequence_length": args.length,
        "mutation_rate": args.mutation_rate,
        "recombination_rate": args.recombination_rate,
        "ancestry_model": "hudson",
        "mutation_model": "InfiniteSites(NUCLEOTIDES)",
        "discrete_genome": False,
        "sites_before_modern_ascertainment": full.num_sites,
        "modern_polymorphic_sites": modern_arg.num_sites,
        "num_trees": modern_arg.num_trees,
        "msprime_version": msprime.__version__,
        "tskit_version": tskit.__version__,
        "analysis_arg": "known_modern_arg.trees (ancient haplotype removed)",
        "truth_tree_sequence": "truth_with_ancient.trees",
        "combined_vcf": "all_samples.vcf.gz",
    }
    (tmp / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    tmp.replace(out)
    print(json.dumps(metadata, indent=2), flush=True)


def vcf_calls(path: Path, n_modern: int):
    calls = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                expected = [f"modern_{i:02d}" for i in range(1, n_modern + 1)] + ["ancient_01"]
                if header[9:] != expected:
                    raise ValueError(f"unexpected VCF sample order in {path}")
                continue
            fields = line.rstrip("\n").split("\t")
            fmt = fields[8].split(":")
            gt_i = fmt.index("GT")
            gt = [sample.split(":")[gt_i] for sample in fields[9:]]
            if any(x not in {"0", "1"} for x in gt):
                raise ValueError(f"non-binary or missing call at site ID {fields[2]}")
            calls[int(fields[2])] = (sum(x == "1" for x in gt[:n_modern]), int(gt[-1]))
    return calls


def infer(args: argparse.Namespace) -> None:
    directory = replicate_dir(args.input_root, args.replicate)
    metadata = json.loads((directory / "metadata.json").read_text())
    ts = tskit.load(directory / "known_modern_arg.trees")
    if ts.num_samples != args.n_modern:
        raise ValueError(f"modern ARG has {ts.num_samples} samples")
    calls = vcf_calls(directory / "all_samples.vcf.gz", args.n_modern)
    if len(calls) != ts.num_sites:
        raise ValueError(f"VCF has {len(calls)} sites but ARG has {ts.num_sites}")
    tab = load_table(args.frequency_table)
    if tab["n_sample"] != args.n_modern:
        raise ValueError("frequency-table panel size mismatch")
    grid = np.asarray(tab["Tgrid"], dtype=np.float64)
    ll = np.zeros(len(grid), dtype=np.float64)
    counts = {"arg_sites": ts.num_sites, "used": 0, "bad_site_mutation_count": 0,
              "vcf_arg_count_mismatch": 0, "age_filtered": 0,
              "root_mutation": 0, "nonfinite_frequency": 0}
    mutation_ages = []
    cutoff = args.mutation_age_max_tau * 2.0 * float(metadata["Ne"])
    for site in ts.sites():
        if len(site.mutations) != 1:
            counts["bad_site_mutation_count"] += 1
            continue
        mutation = site.mutations[0]
        if tskit.is_unknown_time(mutation.time):
            raise ValueError(f"unknown time for mutation {mutation.id}")
        d0, ancient_genotype = calls[site.id]
        tree = ts.at(site.position)
        tree_count = int(tree.num_samples(mutation.node))
        if d0 != tree_count:
            counts["vcf_arg_count_mismatch"] += 1
            continue
        if args.mutation_age_source == "exact":
            t_lo = t_hi = float(mutation.time)
            if t_lo >= cutoff or t_lo < float(tab["age"][0]):
                counts["age_filtered"] += 1
                continue
        else:
            parent = tree.parent(mutation.node)
            if parent == tskit.NULL:
                counts["root_mutation"] += 1
                continue
            t_lo = float(ts.node(mutation.node).time)
            t_hi = min(float(ts.node(parent).time), cutoff)
            if t_lo >= cutoff or t_hi <= t_lo:
                counts["age_filtered"] += 1
                continue
        phi = phi_lookup(
            tab, d0, t_lo, t_hi,
            n_called=args.n_modern, marginalise=args.marginalise,
        )
        if phi is None or not np.all(np.isfinite(phi)):
            counts["nonfinite_frequency"] += 1
            continue
        q = np.clip(
            (1.0 - args.epsilon) * phi + args.epsilon * (1.0 - phi),
            1e-300, 1.0,
        )
        ll += np.log(q) if ancient_genotype else np.log(
            np.clip(1.0 - q, 1e-300, 1.0)
        )
        mutation_ages.append(float(mutation.time))
        counts["used"] += 1
    if counts["used"] == 0:
        raise ValueError("no sites survived exact-time inference")
    summary, density = summarize(grid, ll)
    result = {
        "replicate": args.replicate,
        "true_T": float(metadata["true_ancient_age"]),
        **summary,
        "contains_true_95": int(summary["ci95_lower_T"] <= metadata["true_ancient_age"]
                                <= summary["ci95_upper_T"]),
        "counts": counts,
        "epsilon": args.epsilon,
        "mutation_age_source": (
            "msprime Mutation.time (exact point age)"
            if args.mutation_age_source == "exact"
            else "known ARG child-parent edge interval"
        ),
        "mutation_age_marginalisation": args.marginalise,
        "mutation_age_cutoff_generations": cutoff,
        "mutation_age_min": float(np.min(mutation_ages)),
        "mutation_age_max": float(np.max(mutation_ages)),
        "frequency_lookup": "posterior_sample_age_infer.phi_lookup",
        "posterior_summary": "posterior_sample_age_infer.summarize",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_dir / f"replicate_{args.replicate:03d}.npz",
        Tgrid=grid, log_likelihood=ll, density=density,
    )
    (args.output_dir / f"replicate_{args.replicate:03d}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2), flush=True)


def merge(args: argparse.Namespace) -> None:
    rows = []
    for replicate in range(1, args.n_replicates + 1):
        path = args.input_dir / f"replicate_{replicate:03d}.json"
        if not path.is_file():
            raise ValueError(f"missing inference result {path}")
        row = json.loads(path.read_text())
        rows.append({k: v for k, v in row.items() if k != "counts"})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with (args.output_dir / "age_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    truth = np.asarray([r["true_T"] for r in rows])
    maps = np.asarray([r["map_T"] for r in rows])
    means = np.asarray([r["mean_T"] for r in rows])
    lower = np.asarray([r["ci95_lower_T"] for r in rows])
    upper = np.asarray([r["ci95_upper_T"] for r in rows])
    residual = maps - truth
    summary = {
        "n_replicates": len(rows),
        "model": rows[0]["mutation_age_source"] + " -> production diffusion frequency lookup",
        "map_bias_estimated_minus_true": float(residual.mean()),
        "map_mae": float(np.abs(residual).mean()),
        "map_rmse": float(np.sqrt(np.mean(residual ** 2))),
        "mean_bias_estimated_minus_true": float(np.mean(means - truth)),
        "map_regression_slope": float(np.polyfit(truth, maps, 1)[0]),
        "map_regression_intercept": float(np.polyfit(truth, maps, 1)[1]),
        "ci95_coverage_count": int(np.sum((lower <= truth) & (truth <= upper))),
        "ci95_coverage_fraction": float(np.mean((lower <= truth) & (truth <= upper))),
        "total_sites_used": int(sum(json.loads(
            (args.input_dir / f"replicate_{i:03d}.json").read_text()
        )["counts"]["used"] for i in range(1, args.n_replicates + 1))),
    }
    (args.output_dir / "run.json").write_text(json.dumps(summary, indent=2) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.plot([0, 10000], [0, 10000], color="0.5", ls="--", lw=1)
    ax.vlines(truth, lower, upper, color="0.65", lw=0.7)
    ax.scatter(truth, maps, s=18, alpha=0.8)
    ax.set(xlabel="True ancient-sample age (generations)",
           ylabel="Posterior MAP age (generations)", xlim=(0, 10000), ylim=(0, 10000))
    fig.tight_layout(); fig.savefig(args.output_dir / "map_vs_true_age.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2), flush=True)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sim = sub.add_parser("simulate")
    sim.add_argument("--output-root", type=Path, required=True)
    sim.add_argument("--replicate", type=int, required=True)
    sim.add_argument("--base-seed", type=int, default=1_700_000_000)
    sim.add_argument("--ne", type=float, default=50_000)
    sim.add_argument("--n-modern", type=int, default=26)
    sim.add_argument("--length", type=float, default=10_000_000)
    sim.add_argument("--mutation-rate", type=float, default=1e-8)
    sim.add_argument("--recombination-rate", type=float, default=1e-8)
    sim.add_argument("--age-min", type=float, default=0.0)
    sim.add_argument("--age-max", type=float, default=10_000.0)
    sim.set_defaults(func=simulate)
    inf = sub.add_parser("infer")
    inf.add_argument("--input-root", type=Path, required=True)
    inf.add_argument("--frequency-table", type=Path, required=True)
    inf.add_argument("--output-dir", type=Path, required=True)
    inf.add_argument("--replicate", type=int, required=True)
    inf.add_argument("--n-modern", type=int, default=26)
    inf.add_argument("--epsilon", type=float, default=0.0)
    inf.add_argument("--mutation-age-max-tau", type=float, default=3.0)
    inf.add_argument("--mutation-age-source", choices=("exact", "edge"), default="exact")
    inf.add_argument("--marginalise", choices=("uniform", "weighted"), default="uniform")
    inf.set_defaults(func=infer)
    mer = sub.add_parser("merge")
    mer.add_argument("--input-dir", type=Path, required=True)
    mer.add_argument("--output-dir", type=Path, required=True)
    mer.add_argument("--n-replicates", type=int, default=100)
    mer.set_defaults(func=merge)
    return parser


def main():
    args = build_parser().parse_args()
    if hasattr(args, "replicate") and args.replicate < 1:
        raise SystemExit("--replicate must be positive")
    args.func(args)


if __name__ == "__main__":
    main()
