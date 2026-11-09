#!/usr/bin/env python3
"""Reconstruct raw SINGER dating and compare it with POLEGON dating.

The singer-snakemake workflow retains the raw SINGER node-time files even when
the final tree sequence is assembled from POLEGON node times.  This utility
reuses the workflow's own merge script, changing only which node-time files are
supplied, so topology, mutations, masks, and chunk boundaries remain fixed.
"""
from __future__ import annotations

import argparse
import json
import os
import runpy
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import tszip


def chunk_inputs(chrom_dir: Path, chromosome: str, draw: int) -> dict[str, list[str]]:
    """Return aligned chunk inputs for one posterior draw."""
    chunk_dir = chrom_dir / "chunks"
    params = sorted(chunk_dir.glob("*.yaml"), key=lambda path: int(path.stem))
    if not params:
        raise ValueError(f"no chunk YAML files found in {chunk_dir}")
    result = {name: [] for name in ("params", "recombs", "nodes", "muts", "branches")}
    suffixes = {
        "recombs": f"_recombs_{draw}.txt",
        "nodes": f"_nodes_{draw}.txt",
        "muts": f"_muts_{draw}.txt",
        "branches": f"_branches_{draw}.txt",
    }
    for param in params:
        result["params"].append(str(param))
        for name, suffix in suffixes.items():
            path = chunk_dir / f"{param.stem}{suffix}"
            if not path.exists():
                raise FileNotFoundError(path)
            result[name].append(str(path))
        dated = chunk_dir / f"{param.stem}_polegon_{draw}.txt"
        if not dated.exists():
            raise FileNotFoundError(dated)
    return result


def reconstruct(args: argparse.Namespace) -> None:
    if args.output.exists():
        raise FileExistsError(args.output)
    chrom_dir = args.singer_output / args.chromosome
    inputs = chunk_inputs(chrom_dir, args.chromosome, args.draw)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    log = args.output.with_suffix(args.output.suffix + ".log")
    scripts = args.singer_repo / "workflow" / "scripts"
    merge_script = scripts / "merge_chunks.py"
    if not merge_script.exists():
        raise FileNotFoundError(merge_script)
    sys.path.insert(0, str(scripts))
    snakemake = SimpleNamespace(
        rule="merge_raw_singer_chunks",
        config={
            "diagnostic": "raw SINGER node times; POLEGON bypassed",
            "source_output": str(args.singer_output.resolve()),
            "chromosome": args.chromosome,
            "draw": args.draw,
        },
        input=SimpleNamespace(
            **inputs,
            chunks=str(chrom_dir / f"{args.chromosome}.chunks.p"),
            metadata=str(chrom_dir / f"{args.chromosome}.vcf.samples.p"),
            alleles=str(chrom_dir / f"{args.chromosome}.alleles.p"),
            inaccessible=str(chrom_dir / f"{args.chromosome}.inaccessible.p"),
            omitted=str(chrom_dir / f"{args.chromosome}.omitted.p"),
            node_masks=str(chrom_dir / f"{args.chromosome}.nodemask.p"),
        ),
        output=SimpleNamespace(trees=str(args.output)),
        log=SimpleNamespace(log=str(log)),
        params=SimpleNamespace(
            min_branch_length=args.min_branch_length,
            stratify=args.stratify_by,
            version={"pipeline": "raw-SINGER diagnostic", "polegon": "bypassed"},
            repolarise=True,
            delete_genealogical_gaps=True,
            delete_masked_intervals=False,
            delete_imputed_intervals=True,
        ),
    )
    runpy.run_path(str(merge_script), init_globals={"snakemake": snakemake})


def write_fai(args: argparse.Namespace) -> None:
    if args.output.exists():
        raise FileExistsError(args.output)
    tree_sequence = tszip.decompress(args.tree)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(f"{args.chromosome}\t{int(tree_sequence.sequence_length)}\n")


def ratio_summary(numerator, denominator) -> dict[str, float]:
    keep = np.logical_and(denominator > 0, numerator > 0)
    values = numerator[keep] / denominator[keep]
    return {
        "q05": float(np.quantile(values, 0.05)),
        "median": float(np.median(values)),
        "q95": float(np.quantile(values, 0.95)),
    }


def chunk_time_metrics(chrom_dir: Path, draw: int) -> dict[str, object]:
    """Compare times on the shared SINGER chunk node IDs and branch endpoints."""
    inputs = chunk_inputs(chrom_dir, chrom_dir.name, draw)
    raw_nodes = []
    dated_nodes = []
    raw_edges = []
    dated_edges = []
    spans = []
    nonpositive_raw = 0
    nonpositive_dated = 0
    for raw_path, branch_path in zip(inputs["nodes"], inputs["branches"]):
        raw = np.atleast_1d(np.loadtxt(raw_path))
        dated_path = Path(raw_path).with_name(Path(raw_path).name.replace("_nodes_", "_polegon_"))
        dated = np.atleast_1d(np.loadtxt(dated_path))
        if raw.shape != dated.shape:
            raise ValueError(f"raw and POLEGON node arrays differ in shape: {raw_path}")
        nonsample = raw > 0
        raw_nodes.append(raw[nonsample])
        dated_nodes.append(dated[nonsample])
        branches = np.atleast_2d(np.loadtxt(branch_path))
        branches = branches[branches[:, 2] >= 0]
        parent = branches[:, 2].astype(np.int64)
        child = branches[:, 3].astype(np.int64)
        raw_length = raw[parent] - raw[child]
        dated_length = dated[parent] - dated[child]
        nonpositive_raw += int(np.sum(raw_length <= 0))
        nonpositive_dated += int(np.sum(dated_length <= 0))
        keep = np.logical_and(raw_length > 0, dated_length > 0)
        raw_edges.append(raw_length[keep])
        dated_edges.append(dated_length[keep])
        spans.append((branches[:, 1] - branches[:, 0])[keep])
    raw_node = np.concatenate(raw_nodes)
    dated_node = np.concatenate(dated_nodes)
    raw_edge = np.concatenate(raw_edges)
    dated_edge = np.concatenate(dated_edges)
    span = np.concatenate(spans)
    return {
        "comparison_level": "retained SINGER chunk node IDs and shared branch endpoints",
        "chunks": len(inputs["params"]),
        "internal_node_observations": int(raw_node.size),
        "positive_edge_observations": int(raw_edge.size),
        "nonpositive_raw_edge_lengths": nonpositive_raw,
        "nonpositive_polegon_edge_lengths": nonpositive_dated,
        "node_time_ratio_polegon_over_raw": ratio_summary(dated_node, raw_node),
        "edge_length_ratio_polegon_over_raw": ratio_summary(dated_edge, raw_edge),
        "span_weighted_mean_edge_length": {
            "raw_singer": float(np.average(raw_edge, weights=span)),
            "polegon": float(np.average(dated_edge, weights=span)),
        },
    }


def tree_sequence_counts(raw, dated) -> dict[str, object]:
    """Report output dimensions without assuming stable internal node IDs."""
    return {
        "raw_singer": {
            "nodes": raw.num_nodes, "edges": raw.num_edges, "trees": raw.num_trees,
            "sites": raw.num_sites, "mutations": raw.num_mutations,
        },
        "polegon": {
            "nodes": dated.num_nodes, "edges": dated.num_edges, "trees": dated.num_trees,
            "sites": dated.num_sites, "mutations": dated.num_mutations,
        },
        "note": "Final node IDs are not compared; dating-dependent table sorting can change their representation.",
    }


def calibration_metrics(path: Path) -> dict[str, object]:
    result = json.loads((path / "run.json").read_text())
    joint = result["joint_calibration"]["joint_map"]
    cv = result["cross_validation_metrics"]["joint_scale_epsilon"]
    return {
        "time_scale": joint["time_scale"],
        "epsilon": joint["epsilon"],
        "mean_bias": cv["mean_bias"],
        "mae": cv["mae"],
        "rmse": cv["rmse"],
        "regression_slope": cv.get("regression_slope"),
    }


def summarize(args: argparse.Namespace) -> None:
    if args.output.exists():
        raise FileExistsError(args.output)
    summary = {
        "interpretation_gate": (
            "POLEGON over-correction is supported if raw SINGER prefers a scale near 1 "
            "while the fixed-topology POLEGON tree prefers a substantially smaller scale."
        ),
        "chunk_time_comparison": chunk_time_metrics(
            args.singer_output / args.chromosome, args.draw
        ),
        "tree_sequence_counts": tree_sequence_counts(
            tszip.decompress(args.raw_tree), tszip.decompress(args.polegon_tree)
        ),
        "calibration": {
            "raw_singer": calibration_metrics(args.raw_calibration),
            "polegon": calibration_metrics(args.polegon_calibration),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(summary, indent=2) + "\n")
    os.replace(temporary, args.output)
    print(json.dumps(summary, indent=2), flush=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("reconstruct")
    build.add_argument("--singer-output", type=Path, required=True)
    build.add_argument("--singer-repo", type=Path, default=Path("/home/jri/src/singer-snakemake"))
    build.add_argument("--chromosome", default="chr1")
    build.add_argument("--draw", type=int, default=50)
    build.add_argument("--min-branch-length", type=float, default=1e-7)
    build.add_argument("--stratify-by", default="population")
    build.add_argument("--output", type=Path, required=True)
    build.set_defaults(func=reconstruct)

    fai = subparsers.add_parser("write-fai")
    fai.add_argument("--tree", type=Path, required=True)
    fai.add_argument("--chromosome", default="chr1")
    fai.add_argument("--output", type=Path, required=True)
    fai.set_defaults(func=write_fai)

    report = subparsers.add_parser("summarize")
    report.add_argument("--singer-output", type=Path, required=True)
    report.add_argument("--chromosome", default="chr1")
    report.add_argument("--draw", type=int, default=50)
    report.add_argument("--raw-tree", type=Path, required=True)
    report.add_argument("--polegon-tree", type=Path, required=True)
    report.add_argument("--raw-calibration", type=Path, required=True)
    report.add_argument("--polegon-calibration", type=Path, required=True)
    report.add_argument("--output", type=Path, required=True)
    report.set_defaults(func=summarize)
    return parser.parse_args(argv)


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
