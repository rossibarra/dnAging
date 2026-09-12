#!/usr/bin/env python3
"""Calibrate direct ancient-lineage insertion on the independent T9 loci."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import tskit

from direct_frequency_age_infer import discrete_quantile, normalize_log_likelihood
from insertion_likelihood import derived_probability_uniform_edge_grid
from slim_edge_interval_validation import mutation_slim_ids


def run_chunk(args):
    diagnostic_path = args.diagnostic_dir / f"chunk_{args.chunk_id:04d}.npz"
    with np.load(diagnostic_path, allow_pickle=True) as archive:
        diagnostics = archive["diagnostics"].tolist()
        grid = np.asarray(archive["candidate_ages"], dtype=float)
    likelihood = np.zeros((7, len(grid)), dtype=float)
    counts = {"examined": 0, "edge_count_mismatch": 0, "used": 0}
    for row in diagnostics:
        replicate, focal_id, d0, edge_lower, edge_upper, *genotypes = row
        replicate = int(replicate)
        seed = args.base_seed + replicate
        stem = f"replicate_{replicate:06d}.seed_{seed}"
        ts = tskit.load(args.input_dir / "trees" / f"{stem}.trees")
        matches = [m for m in ts.mutations()
                   if str(focal_id) in mutation_slim_ids(m)]
        if len(matches) != 1:
            raise ValueError(
                f"replicate {replicate}: focal mutation has {len(matches)} matches"
            )
        mutation = matches[0]
        tree = ts.at(ts.site(mutation.site).position)
        counts["examined"] += 1
        if tree.num_samples(mutation.node) != int(d0):
            counts["edge_count_mismatch"] += 1
            continue
        probability = derived_probability_uniform_edge_grid(
            tree, mutation.node, grid, float(edge_lower), float(edge_upper), args.ne
        )
        probability = np.clip(probability, 1e-300, 1 - 1e-15)
        for sample_index, genotype in enumerate(genotypes):
            likelihood[sample_index] += (
                np.log(probability) if int(genotype)
                else np.log1p(-probability)
            )
        counts["used"] += 1
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_dir / f"chunk_{args.chunk_id:04d}.npz",
        log_likelihood=likelihood,
        candidate_ages=grid,
        counts=json.dumps(counts),
    )
    print(json.dumps(counts))


def merge(args):
    likelihood = None
    grid = None
    counts = {}
    for chunk_id in range(args.n_chunks):
        path = args.chunk_dir / f"chunk_{chunk_id:04d}.npz"
        if not path.is_file():
            raise ValueError(f"missing {path}")
        with np.load(path) as archive:
            if likelihood is None:
                likelihood = np.asarray(archive["log_likelihood"], dtype=float)
                grid = np.asarray(archive["candidate_ages"], dtype=float)
            else:
                likelihood += archive["log_likelihood"]
            for key, value in json.loads(str(archive["counts"])).items():
                counts[key] = counts.get(key, 0) + value
    posterior = normalize_log_likelihood(likelihood)
    true_ages = np.asarray([6000, 5000, 4000, 3000, 2000, 1000, 500], dtype=float)
    rows = []
    for index, true_age in enumerate(true_ages):
        density = posterior[index]
        rows.append({
            "true_age": true_age,
            "map_age": float(grid[np.argmax(density)]),
            "posterior_mean_age": float(np.sum(grid * density)),
            "ci95_lower": discrete_quantile(grid, density, 0.025),
            "ci95_upper": discrete_quantile(grid, density, 0.975),
        })
    result = {
        "model": "direct ancient-lineage insertion conditional on modern tree",
        "ne": args.ne,
        "epsilon": 0.0,
        "counts": counts,
        "map_bias_estimated_minus_true": float(np.mean(
            [row["map_age"] - row["true_age"] for row in rows]
        )),
        "all_true_ages_in_ci95": bool(all(
            row["ci95_lower"] <= row["true_age"] <= row["ci95_upper"]
            for row in rows
        )),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    with (args.output_dir / "age_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(args.output_dir / "age_posteriors.npz",
                        candidate_ages=grid, posterior=posterior,
                        log_likelihood=likelihood, true_ages=true_ages)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 6))
    maps = np.asarray([row["map_age"] for row in rows])
    ax.plot([0, 7000], [0, 7000], color="0.5", ls="--")
    ax.scatter(true_ages, maps)
    ax.set(xlabel="True sample age (generations)",
           ylabel="Insertion posterior MAP age (generations)",
           xlim=(0, 7000), ylim=(0, 7000))
    fig.tight_layout()
    fig.savefig(args.output_dir / "map_vs_true_age.png", dpi=180)
    plt.close(fig)
    print(json.dumps(result, indent=2))


def parser():
    top = argparse.ArgumentParser(description=__doc__)
    sub = top.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run-chunk")
    run.add_argument("--input-dir", type=Path, required=True)
    run.add_argument("--diagnostic-dir", type=Path, required=True)
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--chunk-id", type=int, required=True)
    run.add_argument("--base-seed", type=int, default=202609090000)
    run.add_argument("--ne", type=float, default=10000)
    run.set_defaults(func=run_chunk)
    combine = sub.add_parser("merge")
    combine.add_argument("--chunk-dir", type=Path, required=True)
    combine.add_argument("--output-dir", type=Path, required=True)
    combine.add_argument("--n-chunks", type=int, default=100)
    combine.add_argument("--ne", type=float, default=10000)
    combine.set_defaults(func=merge)
    return top


if __name__ == "__main__":
    parsed = parser().parse_args()
    parsed.func(parsed)
