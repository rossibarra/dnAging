#!/usr/bin/env python3
"""Apply direct ancient-lineage insertion to the 100-replicate msprime suite."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import tskit

from insertion_likelihood import PiecewiseConstantNe, derived_probability_uniform_edge_grid
from posterior_sample_age_infer import summarize


def replicate_dir(root: Path, replicate: int) -> Path:
    return root / f"replicate_{replicate:03d}"


def vcf_call_matrix(path: Path, n_modern: int):
    calls = {}
    ancient_names = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                ancient_names = header[9 + n_modern:]
                if not ancient_names:
                    raise ValueError("VCF has no ancient sample columns")
                continue
            fields = line.rstrip("\n").split("\t")
            gt_index = fields[8].split(":").index("GT")
            genotypes = [sample.split(":")[gt_index] for sample in fields[9:]]
            if any(gt not in {"0", "1"} for gt in genotypes):
                raise ValueError(f"non-binary call at site {fields[2]}")
            calls[int(fields[2])] = (
                sum(gt == "1" for gt in genotypes[:n_modern]),
                np.asarray([int(gt) for gt in genotypes[n_modern:]], dtype=np.int8),
            )
    return ancient_names, calls


def vcf_calls(path: Path, n_modern: int):
    """Backward-compatible reader for suites containing one ancient sample."""
    names, matrix = vcf_call_matrix(path, n_modern)
    if len(names) != 1:
        raise ValueError(f"expected one ancient sample, found {len(names)}")
    return {site: (d0, int(genotypes[0]))
            for site, (d0, genotypes) in matrix.items()}


def infer(args):
    directory = replicate_dir(args.input_root, args.replicate)
    metadata = json.loads((directory / "metadata.json").read_text())
    ts = tskit.load(directory / "known_modern_arg.trees")
    ancient_names, matrix = vcf_call_matrix(
        directory / "all_samples.vcf.gz", args.n_modern)
    if args.ancient_sample not in ancient_names:
        raise ValueError(f"ancient sample {args.ancient_sample!r} not in VCF")
    ancient_index = ancient_names.index(args.ancient_sample)
    calls = {site: (d0, int(genotypes[ancient_index]))
             for site, (d0, genotypes) in matrix.items()}
    ne = (PiecewiseConstantNe.from_tsv(
              args.ne_table, args.ne_series, args.ne_extrapolate_last)
          if args.ne_table else float(metadata["Ne"]))
    grid = np.arange(args.age_min, args.age_max + 0.5 * args.age_step,
                     args.age_step, dtype=float)
    log_likelihood = np.zeros(len(grid), dtype=float)
    counts = {"arg_sites": ts.num_sites, "used": 0, "root_mutation": 0,
              "bad_site_mutation_count": 0, "count_mismatch": 0,
              "nonfinite_probability": 0}
    for site_index, site in enumerate(ts.sites(), start=1):
        if len(site.mutations) != 1:
            counts["bad_site_mutation_count"] += 1
            continue
        mutation = site.mutations[0]
        d0, ancient_genotype = calls[site.id]
        tree = ts.at(site.position)
        if int(tree.num_samples(mutation.node)) != d0:
            counts["count_mismatch"] += 1
            continue
        parent = tree.parent(mutation.node)
        if parent == tskit.NULL:
            counts["root_mutation"] += 1
            continue
        lower = float(ts.node(mutation.node).time)
        upper = float(ts.node(parent).time)
        probability = derived_probability_uniform_edge_grid(
            tree, mutation.node, grid, lower, upper, ne
        )
        if not np.all(np.isfinite(probability)):
            counts["nonfinite_probability"] += 1
            continue
        probability = np.clip(probability, 1e-300, 1 - 1e-15)
        observed_probability = np.clip(
            args.epsilon + (1 - 2 * args.epsilon) * probability,
            1e-300, 1 - 1e-15)
        log_likelihood += (np.log(observed_probability) if ancient_genotype
                           else np.log1p(-observed_probability))
        counts["used"] += 1
        if site_index % 10000 == 0:
            print(f"replicate {args.replicate}: {site_index}/{ts.num_sites}",
                  flush=True)
    summary, density = summarize(grid, log_likelihood)
    if "true_ancient_ages" in metadata:
        truth = float(metadata["true_ancient_ages"][ancient_index])
    else:
        truth = float(metadata["true_ancient_age"])
    result = {
        "replicate": args.replicate,
        "true_T": truth,
        **summary,
        "contains_true_95": int(summary["ci95_lower_T"] <= truth
                                <= summary["ci95_upper_T"]),
        "counts": counts,
        "model": "direct ancient-lineage insertion conditional on modern tree",
        "Ne": "piecewise" if args.ne_table else float(metadata["Ne"]),
        "epsilon": args.epsilon,
        "ancient_sample": args.ancient_sample,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output_dir / f"replicate_{args.replicate:03d}.npz",
                        Tgrid=grid, log_likelihood=log_likelihood, density=density)
    (args.output_dir / f"replicate_{args.replicate:03d}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2), flush=True)


def merge(args):
    rows = []
    for replicate in range(1, args.n_replicates + 1):
        path = args.input_dir / f"replicate_{replicate:03d}.json"
        if not path.is_file():
            raise ValueError(f"missing {path}")
        rows.append(json.loads(path.read_text()))
    truth = np.asarray([row["true_T"] for row in rows])
    maps = np.asarray([row["map_T"] for row in rows])
    means = np.asarray([row["mean_T"] for row in rows])
    lower = np.asarray([row["ci95_lower_T"] for row in rows])
    upper = np.asarray([row["ci95_upper_T"] for row in rows])
    residual = maps - truth
    result = {
        "n_replicates": len(rows),
        "model": rows[0]["model"],
        "map_bias_estimated_minus_true": float(residual.mean()),
        "map_mae": float(np.abs(residual).mean()),
        "map_rmse": float(np.sqrt(np.mean(residual ** 2))),
        "mean_bias_estimated_minus_true": float(np.mean(means - truth)),
        "map_regression_slope": float(np.polyfit(truth, maps, 1)[0]),
        "map_regression_intercept": float(np.polyfit(truth, maps, 1)[1]),
        "ci95_coverage_count": int(np.sum((lower <= truth) & (truth <= upper))),
        "ci95_coverage_fraction": float(np.mean((lower <= truth) & (truth <= upper))),
        "total_sites_used": int(sum(row["counts"]["used"] for row in rows)),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    fields = [key for key in rows[0] if key != "counts"]
    with (args.output_dir / "age_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.plot([0, 10000], [0, 10000], color="0.5", ls="--")
    ax.vlines(truth, lower, upper, color="0.7", lw=0.7)
    ax.scatter(truth, maps, s=18)
    ax.set(xlabel="True ancient-sample age (generations)",
           ylabel="Insertion posterior MAP age (generations)",
           xlim=(0, 10000), ylim=(0, 10000))
    fig.tight_layout()
    fig.savefig(args.output_dir / "map_vs_true_age.png", dpi=180)
    plt.close(fig)
    print(json.dumps(result, indent=2))


def parser():
    top = argparse.ArgumentParser(description=__doc__)
    sub = top.add_subparsers(dest="command", required=True)
    run = sub.add_parser("infer")
    run.add_argument("--input-root", type=Path, required=True)
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--replicate", type=int, required=True)
    run.add_argument("--n-modern", type=int, default=26)
    run.add_argument("--ancient-sample", default="ancient_01")
    run.add_argument("--ne-table", type=Path)
    run.add_argument("--ne-series")
    run.add_argument("--ne-extrapolate-last", action="store_true",
                     help="treat the final Ne epoch as ancestral and unbounded")
    run.add_argument("--epsilon", type=float, default=0.0)
    run.add_argument("--age-min", type=float, default=0)
    run.add_argument("--age-max", type=float, default=10000)
    run.add_argument("--age-step", type=float, default=20)
    run.set_defaults(func=infer)
    combine = sub.add_parser("merge")
    combine.add_argument("--input-dir", type=Path, required=True)
    combine.add_argument("--output-dir", type=Path, required=True)
    combine.add_argument("--n-replicates", type=int, default=100)
    combine.set_defaults(func=merge)
    return top


if __name__ == "__main__":
    parsed = parser().parse_args()
    parsed.func(parsed)
