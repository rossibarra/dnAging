#!/usr/bin/env python3
"""Validate focal SLiM mutation placement within its modern-only ARG edge."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import tskit
from scipy import stats

from slim_edge_interval_validation import mutation_slim_ids, read_samples

FOCAL_XOR = np.uint64(0x5EED_A6E)
PIT_XOR = np.uint64(0x0ED6E71F)


def focal_for_replicate(rows, seed):
    modern = [row for row in rows if row["role"] == "modern"]
    ids, counts = np.unique(
        [row["mutation_id"] for row in modern if row["mutation_id"] != "NA"],
        return_counts=True,
    )
    candidates = [(mid, int(n)) for mid, n in zip(ids, counts) if 0 < n < 26]
    if not candidates:
        return None
    rng = np.random.default_rng(np.uint64(seed) ^ FOCAL_XOR)
    return candidates[int(rng.integers(len(candidates)))]


def uniform_summary(values):
    values = np.asarray(values, dtype=float)
    test = stats.kstest(values, "uniform")
    return {
        "n": int(len(values)),
        "mean": float(values.mean()),
        "quantiles": {str(q): float(np.quantile(values, q))
                      for q in (0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)},
        "ks_statistic": float(test.statistic),
        "ks_pvalue": float(test.pvalue),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-dir", type=Path, required=True)
    parser.add_argument("--base-seed", type=int, default=202609090000)
    args = parser.parse_args()

    focal_rows = []
    for chunk_path in sorted(args.chunk_dir.glob("chunk_*.npz")):
        with np.load(chunk_path, allow_pickle=True) as chunk:
            focal_rows.extend(chunk["diagnostics"].tolist())
    if len(focal_rows) == 0:
        raise ValueError(f"no inference diagnostics found in {args.chunk_dir}")

    records = []
    for focal_row in focal_rows:
        replicate, focal_id, d0 = int(focal_row[0]), str(focal_row[1]), int(focal_row[2])
        seed = args.base_seed + replicate
        stem = f"replicate_{replicate:06d}.seed_{seed}"
        ts = tskit.load(args.input_dir / "trees" / f"{stem}.trees")
        matches = [m for m in ts.mutations() if focal_id in mutation_slim_ids(m)]
        if len(matches) != 1:
            raise ValueError(
                f"replicate {replicate}: focal ID {focal_id} has {len(matches)} rows"
            )
        mutation = matches[0]
        tree = ts.at(ts.site(mutation.site).position)
        parent = tree.parent(mutation.node)
        if parent == tskit.NULL:
            raise ValueError(f"replicate {replicate}: mutation above root")
        lo = float(ts.node(mutation.node).time)
        hi = float(ts.node(parent).time)
        age = float(mutation.time)
        width = hi - lo
        edge_descendants = int(tree.num_samples(mutation.node))
        first_tick = int(np.ceil(lo))
        last_tick = int(np.ceil(hi)) - 1
        n_ticks = last_tick - first_tick + 1
        if n_ticks < 1 or age != round(age):
            pit = np.nan
        else:
            jitter = np.random.default_rng(np.uint64(seed) ^ PIT_XOR).random()
            pit = (int(round(age)) - first_tick + jitter) / n_ticks
        records.append((replicate, seed, focal_id, d0, edge_descendants,
                        lo, hi, width, age, (age - lo) / width,
                        first_tick, last_tick, n_ticks, pit))

    fields = ("replicate", "seed", "mutation_id", "d0", "edge_descendants",
              "edge_lower", "edge_upper", "edge_width", "mutation_age",
              "normalized_position", "first_tick", "last_tick", "n_ticks",
              "randomized_pit")
    dtype = np.dtype([
        ("replicate", "i8"), ("seed", "i8"), ("mutation_id", "U32"),
        ("d0", "i4"), ("edge_descendants", "i4"), ("edge_lower", "f8"),
        ("edge_upper", "f8"), ("edge_width", "f8"), ("mutation_age", "f8"),
        ("normalized_position", "f8"), ("first_tick", "i8"),
        ("last_tick", "i8"), ("n_ticks", "i8"), ("randomized_pit", "f8"),
    ])
    data = np.asarray(records, dtype=dtype)
    contained = ((data["mutation_age"] >= data["edge_lower"])
                 & (data["mutation_age"] < data["edge_upper"]))
    count_matches = data["d0"] == data["edge_descendants"]
    finite = np.isfinite(data["randomized_pit"])
    result = {
        "n_retained": int(len(data)),
        "n_contained_lower_inclusive_upper_exclusive": int(contained.sum()),
        "n_outside": int((~contained).sum()),
        "n_d0_equals_edge_descendants": int(count_matches.sum()),
        "n_recurrent_overwrite_mismatch": int((~count_matches).sum()),
        "all_loci_discrete_randomized_pit": uniform_summary(
            data["randomized_pit"][finite]
        ),
        "count_matched_loci_discrete_randomized_pit": uniform_summary(
            data["randomized_pit"][finite & count_matches]
        ),
    }
    width_edges = np.quantile(data["edge_width"][count_matches], [0, .25, .5, .75, 1])
    result["count_matched_width_quartiles"] = []
    for j in range(4):
        take = (count_matches & finite & (data["edge_width"] >= width_edges[j])
                & ((data["edge_width"] <= width_edges[j + 1]) if j == 3
                   else (data["edge_width"] < width_edges[j + 1])))
        row = {"bin": j + 1, "width_lower": float(width_edges[j]),
               "width_upper": float(width_edges[j + 1])}
        row.update(uniform_summary(data["randomized_pit"][take]))
        result["count_matched_width_quartiles"].append(row)
    result["count_matched_d0_strata"] = []
    for label, lower, upper in (("1", 1, 1), ("2", 2, 2), ("3-5", 3, 5),
                                ("6-12", 6, 12), ("13-25", 13, 25)):
        take = count_matches & finite & (data["d0"] >= lower) & (data["d0"] <= upper)
        if np.any(take):
            row = {"d0": label}
            row.update(uniform_summary(data["randomized_pit"][take]))
            result["count_matched_d0_strata"].append(row)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    with (args.output_dir / "mutation_edge_positions.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(data.tolist())

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    values = data["randomized_pit"][finite & count_matches]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axes[0].hist(values, bins=np.linspace(0, 1, 21), density=True,
                 color="0.35", edgecolor="white")
    axes[0].axhline(1, color="C1", ls="--")
    axes[0].set(xlabel="Randomized within-edge position", ylabel="Density")
    ordered = np.sort(values)
    axes[1].plot(ordered, np.arange(1, len(ordered) + 1) / len(ordered))
    axes[1].plot([0, 1], [0, 1], color="C1", ls="--")
    axes[1].set(xlabel="Randomized within-edge position", ylabel="Empirical CDF")
    fig.savefig(args.output_dir / "edge_position_uniformity.png", dpi=180)
    plt.close(fig)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
