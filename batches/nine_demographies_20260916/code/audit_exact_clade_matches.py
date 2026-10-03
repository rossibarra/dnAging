#!/usr/bin/env python3
"""Audit exact inferred/true mutation clades and build a strict common mask."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import tskit
import tszip

from betabinom_real_data import orientation
from insertion_real_data import call_files, chromosome_layout


def _sample_ranks(ts: tskit.TreeSequence) -> np.ndarray:
    ranks = np.full(ts.num_nodes, -1, dtype=np.int32)
    ranks[ts.samples()] = np.arange(ts.num_samples, dtype=np.int32)
    return ranks


def _descendant_ranks(tree: tskit.Tree, node: int, ranks: np.ndarray) -> np.ndarray:
    values = ranks[np.fromiter(tree.samples(node), dtype=np.int32)]
    if np.any(values < 0):
        raise ValueError("mutation clade contains a sample outside the modern panel")
    return np.sort(values)


def _oriented_carriers(
    tree: tskit.Tree, node: int, ranks: np.ndarray, derived_is_alt: bool, n: int
) -> np.ndarray:
    carriers = _descendant_ranks(tree, node, ranks)
    if derived_is_alt:
        return carriers
    keep = np.ones(n, dtype=bool)
    keep[carriers] = False
    return np.flatnonzero(keep).astype(np.int32)


def _position_maps(input_dir: Path, chromosomes: list[str]) -> dict[str, dict[int, int]]:
    result = {}
    for chrom in chromosomes:
        path = input_dir / f"{chrom}.position_map.tsv"
        mapping = {}
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                mapping[int(row["singer_pos"])] = int(row["site_id"])
        result[chrom] = mapping
    return result


def _summary(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"n": 0, "mean": None, "median": None, "q025": None, "q975": None}
    x = np.asarray(values, dtype=float)
    return {
        "n": int(x.size),
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "q025": float(np.quantile(x, 0.025)),
        "q975": float(np.quantile(x, 0.975)),
    }


def audit_draw(args: argparse.Namespace) -> None:
    inferred = tszip.decompress(args.tree)
    layout = chromosome_layout(inferred)
    inferred_ranks = _sample_ranks(inferred)
    if inferred.num_samples != args.n_modern:
        raise ValueError(f"inferred ARG has {inferred.num_samples} samples, expected {args.n_modern}")

    calls = call_files(args.calls)
    chromosomes = sorted(f"chr{c}" if f"chr{c}" in layout else c for c in calls)
    maps = _position_maps(args.singer_input, chromosomes)
    positions = np.asarray(inferred.tables.sites.position)
    output = {"draw_id": args.draw_id, "tree": str(args.tree.resolve()), "chromosomes": {}}
    args.output.mkdir(parents=True, exist_ok=False)

    for chrom in chromosomes:
        suffix = chrom[3:] if chrom.startswith("chr") else chrom
        data = np.load(calls[suffix], allow_pickle=False)
        true_path = args.simulations / f"replicate_{int(suffix):03d}" / "known_modern_arg.trees"
        truth = tskit.load(true_path)
        if truth.num_samples != inferred.num_samples:
            raise ValueError(f"sample-count mismatch for {chrom}: {truth.num_samples} vs {inferred.num_samples}")
        truth_ranks = _sample_ranks(truth)
        offset, length = layout[chrom]
        exact_positions = []
        child_errors = []
        parent_errors = []
        width_errors = []
        jaccards = []
        count_equal = 0
        missing_truth = 0
        invalid = 0

        local_positions = np.asarray(data["position"], dtype=np.int64)
        if args.max_sites_per_chrom is not None:
            local_positions = local_positions[: args.max_sites_per_chrom]
        for row, local in enumerate(local_positions):
            site_id = maps[chrom].get(int(local))
            if site_id is None or site_id >= truth.num_sites:
                missing_truth += 1
                continue
            true_site = truth.site(site_id)
            if len(true_site.mutations) != 1:
                invalid += 1
                continue
            global_pos = float(local) + offset
            inferred_site_id = int(np.searchsorted(positions, global_pos))
            if inferred_site_id >= inferred.num_sites or positions[inferred_site_id] != global_pos:
                invalid += 1
                continue
            inferred_site = inferred.site(inferred_site_id)
            if len(inferred_site.mutations) != 1:
                invalid += 1
                continue

            true_mut = true_site.mutations[0]
            inferred_mut = inferred_site.mutations[0]
            true_tree = truth.at(true_site.position)
            inferred_tree = inferred.at(global_pos)
            true_parent = true_tree.parent(true_mut.node)
            inferred_parent = inferred_tree.parent(inferred_mut.node)
            if true_parent == tskit.NULL or inferred_parent == tskit.NULL:
                invalid += 1
                continue

            inferred_is_alt, _ = orientation(
                inferred_site,
                inferred_mut,
                str(data["ref"][row]),
                str(data["alt"][row]),
            )
            if inferred_is_alt is None:
                invalid += 1
                continue
            true_carriers = _descendant_ranks(true_tree, true_mut.node, truth_ranks)
            inferred_carriers = _oriented_carriers(
                inferred_tree,
                inferred_mut.node,
                inferred_ranks,
                inferred_is_alt,
                inferred.num_samples,
            )
            if true_carriers.size == inferred_carriers.size:
                count_equal += 1
            intersection = np.intersect1d(true_carriers, inferred_carriers, assume_unique=True).size
            union = true_carriers.size + inferred_carriers.size - intersection
            jaccards.append(1.0 if union == 0 else intersection / union)
            if np.array_equal(true_carriers, inferred_carriers):
                exact_positions.append(int(local))
                true_child = true_tree.time(true_mut.node)
                true_parent_time = true_tree.time(true_parent)
                inferred_child = inferred_tree.time(inferred_mut.node)
                inferred_parent_time = inferred_tree.time(inferred_parent)
                child_errors.append(inferred_child - true_child)
                parent_errors.append(inferred_parent_time - true_parent_time)
                width_errors.append(
                    (inferred_parent_time - inferred_child) - (true_parent_time - true_child)
                )

        exact = np.asarray(exact_positions, dtype=np.int64)
        np.save(args.output / f"exact_{chrom}.npy", exact)
        output["chromosomes"][chrom] = {
            "targets": int(local_positions.size),
            "exact": int(exact.size),
            "exact_fraction": float(exact.size / local_positions.size) if local_positions.size else None,
            "derived_count_equal": int(count_equal),
            "missing_truth": int(missing_truth),
            "invalid": int(invalid),
            "jaccard": _summary(jaccards),
            "exact_child_time_error": _summary(child_errors),
            "exact_parent_time_error": _summary(parent_errors),
            "exact_edge_width_error": _summary(width_errors),
        }
        print(chrom, output["chromosomes"][chrom], flush=True)

    (args.output / "audit.json").write_text(json.dumps(output, indent=2) + "\n")


def merge(args: argparse.Namespace) -> None:
    with args.manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("empty draw manifest")
    calls = call_files(args.calls)
    args.output.mkdir(parents=True, exist_ok=False)
    call_output = args.output / "calls"
    call_output.mkdir()
    summary = {"n_draws": len(rows), "draw_ids": [], "chromosomes": {}}

    for suffix, path in sorted(calls.items()):
        chrom = f"chr{suffix}"
        data = np.load(path, allow_pickle=False)
        common = set(np.asarray(data["position"], dtype=np.int64).tolist())
        per_draw = []
        for row in rows:
            draw_id = int(row["draw_id"])
            exact = np.load(args.audit_root / f"draw_{draw_id}" / f"exact_{chrom}.npy")
            exact_set = set(np.asarray(exact, dtype=np.int64).tolist())
            per_draw.append(len(exact_set))
            common.intersection_update(exact_set)
        positions = np.asarray(data["position"], dtype=np.int64)
        keep = np.isin(positions, np.asarray(sorted(common), dtype=np.int64))
        payload = {key: data[key][keep] if data[key].shape[:1] == (positions.size,) else data[key]
                   for key in data.files}
        np.savez_compressed(call_output / f"chr{suffix}.npz", **payload)
        summary["chromosomes"][chrom] = {
            "input_sites": int(positions.size),
            "strict_exact_sites": int(np.sum(keep)),
            "per_draw_exact_min": int(min(per_draw)),
            "per_draw_exact_median": float(np.median(per_draw)),
            "per_draw_exact_max": int(max(per_draw)),
        }

    with (args.output / "draw_manifest.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary["draw_ids"] = [int(row["draw_id"]) for row in rows]
    (args.output / "prepare.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    audit = sub.add_parser("audit-draw")
    audit.add_argument("--tree", type=Path, required=True)
    audit.add_argument("--draw-id", type=int, required=True)
    audit.add_argument("--calls", type=Path, required=True)
    audit.add_argument("--simulations", type=Path, required=True)
    audit.add_argument("--singer-input", type=Path, required=True)
    audit.add_argument("--n-modern", type=int, default=100)
    audit.add_argument("--max-sites-per-chrom", type=int)
    audit.add_argument("--output", type=Path, required=True)
    audit.set_defaults(func=audit_draw)

    merge_parser = sub.add_parser("merge")
    merge_parser.add_argument("--manifest", type=Path, required=True)
    merge_parser.add_argument("--calls", type=Path, required=True)
    merge_parser.add_argument("--audit-root", type=Path, required=True)
    merge_parser.add_argument("--output", type=Path, required=True)
    merge_parser.set_defaults(func=merge)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
