#!/usr/bin/env python3
"""Date existing true-topology ARGs with POLEGON and build joint likelihood grids."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import tskit

import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / d) for d in ("pipeline", "validation")]

from insertion_likelihood import PiecewiseConstantNe, derived_probability_uniform_edge_grid
from insertion_real_data import parse_epsilon_grid, parse_time_scale_grid
from radiocarbon_oracle_stage1 import injected_calls, simulation_dirs


def root_stems(ts: tskit.TreeSequence) -> np.ndarray:
    """Return left, right, -1, root rows for each constant-root interval."""
    rows = []
    left = 0.0
    root = tskit.NULL
    for tree in ts.trees():
        next_root = tree.root
        if next_root != root:
            if root != tskit.NULL and tree.interval.left > left:
                rows.append([left, tree.interval.left, tskit.NULL, root])
            left = tree.interval.left
            root = next_root
    if root != tskit.NULL and ts.sequence_length > left:
        rows.append([left, ts.sequence_length, tskit.NULL, root])
    return np.asarray(rows, dtype=float).reshape((-1, 4))


def polegon_tables(ts: tskit.TreeSequence, mutation_rate: float):
    if mutation_rate <= 0:
        raise ValueError("mutation rate must be positive")
    edges = ts.tables.edges
    branches = np.column_stack([
        edges.left * mutation_rate,
        edges.right * mutation_rate,
        edges.parent,
        edges.child,
    ])
    stems = root_stems(ts)
    stems[:, :2] *= mutation_rate
    branches = np.concatenate([branches, stems])
    mutation_edges = ts.mutations_edge
    keep = mutation_edges != tskit.NULL
    mutations = np.column_stack([
        ts.sites_position[ts.mutations_site[keep]] * mutation_rate,
        ts.edges_child[mutation_edges[keep]],
        ts.edges_parent[mutation_edges[keep]],
        np.ones(np.sum(keep), dtype=np.int8),
    ])
    return ts.tables.nodes.time.copy(), branches, mutations


def replace_node_times(ts: tskit.TreeSequence, node_times) -> tskit.TreeSequence:
    node_times = np.asarray(node_times, dtype=float)
    if node_times.shape != (ts.num_nodes,) or np.any(~np.isfinite(node_times)):
        raise ValueError("POLEGON returned invalid node-time vector")
    samples = ts.samples()
    if np.any(np.abs(node_times[samples]) > 1e-8):
        raise ValueError("POLEGON changed sample times away from zero")
    node_times[samples] = 0.0
    edges = ts.tables.edges
    lengths = node_times[edges.parent] - node_times[edges.child]
    if np.any(lengths <= 0):
        raise ValueError(
            f"POLEGON returned {int(np.sum(lengths <= 0))} nonpositive edge lengths"
        )
    tables = ts.dump_tables()
    tables.nodes.time = node_times
    tables.time_units = "generations"
    # Simulation truth stores explicit mutation times on the original edges.
    # They are invalid after redating, so follow singer-snakemake's merge rule:
    # discard them and recompute midpoint times on the newly dated branches.
    tables.mutations.time = np.full(tables.mutations.num_rows, tskit.UNKNOWN_TIME)
    # Edge order depends on parent time. Replacing node times invalidates the
    # source tree sequence's order even though topology and coordinates are
    # unchanged, so restore tskit's canonical order before reconstruction.
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    tables.provenances.add_row(json.dumps({
        "software": {"name": "POLEGON", "version": "0.1.3-alpha"},
        "parameters": {"diagnostic": "true topology with POLEGON node dating"},
    }))
    return tables.tree_sequence()


def date(args):
    if args.output.exists():
        raise FileExistsError(args.output)
    ts = tskit.load(args.tree)
    nodes, branches, mutations = polegon_tables(ts, args.mutation_rate)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    scratch_root = Path(os.environ.get("TMPDIR", "/tmp"))
    with tempfile.TemporaryDirectory(prefix="true-polegon-", dir=scratch_root) as work:
        prefix = Path(work) / "arg"
        np.savetxt(f"{prefix}_nodes.txt", nodes)
        np.savetxt(f"{prefix}_branches.txt", branches)
        np.savetxt(f"{prefix}_muts.txt", mutations)
        command = [
            str(args.polegon_binary), "-input", str(prefix), "-seed", str(args.seed),
            "-Ne", str(args.ne), "-burn_in", str(args.burn_in),
            "-max_step", str(args.max_step), "-num_samples", str(args.num_samples),
            "-scaling_rep", str(args.scaling_reps), "-m", "1.0", "-thin", str(args.thin),
        ]
        completed = subprocess.run(command, text=True, capture_output=True)
        log_path = args.output.with_suffix(args.output.suffix + ".polegon.log")
        log_path.write_text(completed.stdout + completed.stderr)
        if completed.returncode != 0:
            raise RuntimeError(f"POLEGON failed with exit code {completed.returncode}; see {log_path}")
        dated_times = np.loadtxt(f"{prefix}_new_nodes.txt")
    dated = replace_node_times(ts, dated_times)
    temporary = args.output.with_name(args.output.name + f".tmp.{os.getpid()}")
    dated.dump(temporary)
    os.replace(temporary, args.output)
    metadata = {
        "source_tree": str(args.tree.resolve()),
        "output_tree": str(args.output.resolve()),
        "topology": "true simulated topology",
        "dating": "POLEGON",
        "mutation_rate": args.mutation_rate,
        "Ne": args.ne,
        "seed": args.seed,
        "num_nodes": dated.num_nodes,
        "num_edges": dated.num_edges,
        "num_sites": dated.num_sites,
        "node_time_ratio_median": float(np.median(
            dated.tables.nodes.time[dated.tables.nodes.time > 0]
            / ts.tables.nodes.time[dated.tables.nodes.time > 0]
        )),
    }
    args.output.with_suffix(args.output.suffix + ".json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


def likelihood(args):
    scale_grid = parse_time_scale_grid(args.time_scale_grid)
    epsilon_grid = parse_epsilon_grid(args.epsilon_grid)
    age_grid = np.arange(args.age_min, args.age_max + args.age_step / 2, args.age_step)
    simulations = simulation_dirs(args.run_root)
    total = None
    samples = None
    counts = {}
    for chromosome_index, sim in enumerate(simulations, start=1):
        metadata = json.loads((sim / "metadata.json").read_text())
        tree_path = args.dated_trees / args.tree_pattern.format(
            chromosome=chromosome_index
        )
        ts = tskit.load(tree_path)
        truth = tskit.load(sim / "known_modern_arg.trees")
        names, calls, _, _ = injected_calls(
            sim, args.n_modern, 0.0, args.injection_seed, chromosome_index
        )
        calls_by_position = {
            truth.site(site_id).position: value for site_id, value in calls.items()
        }
        if samples is None:
            samples = names
            total = np.zeros((
                len(scale_grid), len(epsilon_grid), len(simulations),
                len(samples), len(age_grid)
            ))
        elif samples != names:
            raise ValueError("ancient sample order differs among chromosomes")
        ne = PiecewiseConstantNe.from_tsv(
            sim / "constant_ne_epochs.tsv", metadata["ne_table_series"]
        )
        used = skipped = 0
        for site in ts.sites():
            if len(site.mutations) != 1 or site.position not in calls_by_position:
                skipped += 1
                continue
            mutation = site.mutations[0]
            tree = ts.at(site.position)
            d0, observed = calls_by_position[site.position]
            if tree.num_samples(mutation.node) != d0:
                skipped += 1
                continue
            parent = tree.parent(mutation.node)
            if parent == tskit.NULL:
                skipped += 1
                continue
            probabilities = np.asarray([
                derived_probability_uniform_edge_grid(
                    tree, mutation.node, age_grid, tree.time(mutation.node),
                    tree.time(parent), ne, time_scale=scale
                ) for scale in scale_grid
            ])
            eps = epsilon_grid[None, :, None]
            q = np.clip(eps + (1 - 2 * eps) * probabilities[:, None, :], 1e-300, 1 - 1e-15)
            observed = np.asarray(observed, dtype=float)
            total[:, :, chromosome_index - 1] += observed[None, None, :, None] * np.log(q)[:, :, None, :]
            total[:, :, chromosome_index - 1] += (1 - observed)[None, None, :, None] * np.log1p(-q)[:, :, None, :]
            used += 1
        counts[f"chr{chromosome_index}"] = {"used": used, "skipped": skipped}
        print(f"chr{chromosome_index}: {used} sites", flush=True)
    if args.output.exists():
        raise FileExistsError(args.output)
    temporary = args.output.with_name(args.output.name + f".tmp.{os.getpid()}")
    temporary.mkdir(parents=True)
    np.save(temporary / "grid.npy", age_grid)
    np.save(temporary / "time_scale_grid.npy", scale_grid)
    np.save(temporary / "epsilon_grid.npy", epsilon_grid)
    np.save(temporary / "ll_by_chrom_scale_epsilon.npy", total)
    (temporary / "samples.txt").write_text("\n".join(samples) + "\n")
    (temporary / "run.json").write_text(json.dumps({
        "model": "true topology with POLEGON node dates",
        "chromosomes": [f"chr{i}" for i in range(1, len(simulations) + 1)],
        "counts_by_chromosome": counts,
        "time_scale_grid": scale_grid.tolist(),
        "epsilon_grid": epsilon_grid.tolist(),
    }, indent=2) + "\n")
    os.replace(temporary, args.output)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    dating = sub.add_parser("date")
    dating.add_argument("--tree", type=Path, required=True)
    dating.add_argument("--output", type=Path, required=True)
    dating.add_argument("--mutation-rate", type=float, default=1e-8)
    dating.add_argument("--ne", type=float, required=True)
    dating.add_argument("--seed", type=int, required=True)
    dating.add_argument("--burn-in", type=int, default=50)
    dating.add_argument("--max-step", type=float, default=10.0)
    dating.add_argument("--num-samples", type=int, default=100)
    dating.add_argument("--scaling-reps", type=int, default=5)
    dating.add_argument("--thin", type=int, default=10)
    dating.add_argument("--polegon-binary", type=Path, default=Path(
        "/home/jri/src/singer-snakemake/resources/polegon-0.1.3-alpha-linux-x86_64/polegon"
    ))
    dating.set_defaults(func=date)

    grid = sub.add_parser("likelihood")
    grid.add_argument("--run-root", type=Path, required=True)
    grid.add_argument("--dated-trees", type=Path, required=True)
    grid.add_argument(
        "--tree-pattern",
        default="chr{chromosome}.polegon_true_topology.trees",
        help="filename pattern under --dated-trees; use {chromosome} as the index",
    )
    grid.add_argument("--n-modern", type=int, default=100)
    grid.add_argument("--injection-seed", type=int, default=2026110301)
    grid.add_argument("--time-scale-grid", default="0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4")
    grid.add_argument("--epsilon-grid", default="0:0.06:0.01")
    grid.add_argument("--age-min", type=float, default=0)
    grid.add_argument("--age-max", type=float, default=10000)
    grid.add_argument("--age-step", type=float, default=100)
    grid.add_argument("--output", type=Path, required=True)
    grid.set_defaults(func=likelihood)
    return parser.parse_args(argv)


if __name__ == "__main__":
    options = parse_args()
    options.func(options)
