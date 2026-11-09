#!/usr/bin/env python3
"""Re-date true ARGs with the production chunk-by-chunk POLEGON workflow."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import tskit
import yaml

from true_topology_polegon import polegon_tables, replace_node_times


def chunk_parameters(parameter_dir: Path) -> list[tuple[Path, dict]]:
    paths = sorted(parameter_dir.glob("*.yaml"), key=lambda path: int(path.stem))
    if not paths:
        raise ValueError(f"no chunk YAML files found in {parameter_dir}")
    result = []
    previous_end = None
    for path in paths:
        parameters = yaml.safe_load(path.read_text())
        start = int(parameters["singer"]["start"])
        end = int(parameters["singer"]["end"])
        if start >= end:
            raise ValueError(f"invalid chunk interval in {path}: [{start}, {end})")
        if previous_end is not None and start != previous_end:
            raise ValueError(
                f"non-contiguous production chunks: previous end {previous_end}, "
                f"next start {start} in {path}"
            )
        previous_end = end
        result.append((path, parameters))
    return result


def extract_chunk(ts: tskit.TreeSequence, start: int, end: int):
    """Return a compact, coordinate-shifted chunk and its source-node map."""
    if not (0 <= start < end <= ts.sequence_length):
        raise ValueError(f"chunk [{start}, {end}) is outside sequence length {ts.sequence_length}")
    cropped = ts.keep_intervals([[start, end]], simplify=False).trim()
    chunk, node_map = cropped.simplify(
        samples=cropped.samples(),
        keep_unary=True,
        keep_input_roots=True,
        filter_sites=False,
        map_nodes=True,
    )
    if chunk.sequence_length != end - start:
        raise ValueError(
            f"trimmed chunk length {chunk.sequence_length} != requested length {end - start}"
        )
    if not np.array_equal(chunk.samples(), np.arange(ts.num_samples)):
        raise ValueError("chunk simplification did not preserve shared sample IDs")
    return chunk, node_map


def run_polegon(
    chunk: tskit.TreeSequence,
    parameters: dict,
    draw: int,
    mutation_rate: float,
    binary: Path,
    work_dir: Path,
    log_path: Path,
) -> tskit.TreeSequence:
    nodes, branches, mutations = polegon_tables(chunk, mutation_rate)
    prefix = work_dir / "arg"
    np.savetxt(f"{prefix}_nodes.txt", nodes)
    np.savetxt(f"{prefix}_branches.txt", branches)
    np.savetxt(f"{prefix}_muts.txt", mutations)
    polegon = dict(parameters["polegon"])
    seed = int(polegon.pop("seed")) + draw
    command = [str(binary), "-input", str(prefix), "-seed", str(seed)]
    for name, value in polegon.items():
        command.extend([f"-{name}", str(value)])
    command.extend(["-m", "1.0"])
    completed = subprocess.run(command, text=True, capture_output=True)
    log_path.write_text(
        " ".join(command) + "\n" + completed.stdout + completed.stderr
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"POLEGON failed with exit code {completed.returncode}; see {log_path}"
        )
    return replace_node_times(chunk, np.loadtxt(f"{prefix}_new_nodes.txt"))


def merge_chunks(
    chunks: list[tuple[int, tskit.TreeSequence]],
    sequence_length: float,
) -> tskit.TreeSequence:
    """Merge chunks like singer-snakemake: samples shared, internal nodes duplicated."""
    if not chunks:
        raise ValueError("cannot merge zero chunks")
    first = chunks[0][1]
    tables = tskit.TableCollection(sequence_length=sequence_length)
    tables.time_units = "generations"
    for _ in first.samples():
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    n_samples = first.num_samples
    for start, chunk in chunks:
        if chunk.num_samples != n_samples:
            raise ValueError("sample count differs among chunks")
        node_map = np.full(chunk.num_nodes, tskit.NULL, dtype=np.int32)
        node_map[chunk.samples()] = np.arange(n_samples, dtype=np.int32)
        for node in chunk.nodes():
            if not node.is_sample():
                node_map[node.id] = tables.nodes.add_row(time=node.time)
        for edge in chunk.edges():
            tables.edges.add_row(
                left=edge.left + start,
                right=edge.right + start,
                parent=node_map[edge.parent],
                child=node_map[edge.child],
            )
        for site in chunk.sites():
            site_id = tables.sites.add_row(
                position=site.position + start,
                ancestral_state=site.ancestral_state,
            )
            for mutation in site.mutations:
                tables.mutations.add_row(
                    site=site_id,
                    node=node_map[mutation.node],
                    derived_state=mutation.derived_state,
                    time=tskit.UNKNOWN_TIME,
                )
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.compute_mutation_times()
    tables.provenances.add_row(json.dumps({
        "software": {"name": "POLEGON", "version": "0.1.3-alpha"},
        "parameters": {
            "diagnostic": "true topology split, dated, and merged by production chunks",
            "chunks": len(chunks),
        },
    }))
    return tables.tree_sequence()


def dump_atomic(ts: tskit.TreeSequence, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    ts.dump(temporary)
    os.replace(temporary, output)


def run(args: argparse.Namespace) -> None:
    truth = tskit.load(args.tree)
    specs = chunk_parameters(args.parameter_dir)
    if len(specs) != args.expected_chunks:
        raise ValueError(f"expected {args.expected_chunks} chunks, found {len(specs)}")
    true_chunks = []
    dated_chunks = []
    ratios = []
    scratch_root = Path(os.environ.get("TMPDIR", "/tmp"))
    args.chunk_output.mkdir(parents=True, exist_ok=True)
    for index, (parameter_path, parameters) in enumerate(specs, start=1):
        start = int(parameters["singer"]["start"])
        end = int(parameters["singer"]["end"])
        chunk, _ = extract_chunk(truth, start, end)
        log_path = args.chunk_output / f"{index:06d}.polegon.log"
        with tempfile.TemporaryDirectory(
            prefix=f"true-polegon-{index:06d}-", dir=scratch_root
        ) as work:
            dated = run_polegon(
                chunk, parameters, args.draw, args.mutation_rate,
                args.polegon_binary, Path(work), log_path,
            )
        chunk_output = args.chunk_output / f"{index:06d}.polegon.trees"
        dump_atomic(dated, chunk_output)
        positive = chunk.tables.nodes.time > 0
        ratios.extend(
            (dated.tables.nodes.time[positive] / chunk.tables.nodes.time[positive]).tolist()
        )
        true_chunks.append((start, chunk))
        dated_chunks.append((start, dated))
        print(
            f"chunk {index}/{len(specs)} [{start}, {end}): "
            f"{chunk.num_nodes} nodes, {chunk.num_sites} sites",
            flush=True,
        )
    true_merged = merge_chunks(true_chunks, truth.sequence_length)
    dated_merged = merge_chunks(dated_chunks, truth.sequence_length)
    dump_atomic(true_merged, args.true_output)
    dump_atomic(dated_merged, args.polegon_output)
    ratio = np.asarray(ratios)
    metadata = {
        "source_tree": str(args.tree.resolve()),
        "production_parameter_dir": str(args.parameter_dir.resolve()),
        "draw": args.draw,
        "chunks": len(specs),
        "covered_interval": [
            int(specs[0][1]["singer"]["start"]),
            int(specs[-1][1]["singer"]["end"]),
        ],
        "true_merged": {
            "nodes": true_merged.num_nodes,
            "edges": true_merged.num_edges,
            "trees": true_merged.num_trees,
            "sites": true_merged.num_sites,
        },
        "polegon_merged": {
            "nodes": dated_merged.num_nodes,
            "edges": dated_merged.num_edges,
            "trees": dated_merged.num_trees,
            "sites": dated_merged.num_sites,
        },
        "node_time_ratio_polegon_over_truth": {
            "q05": float(np.quantile(ratio, 0.05)),
            "median": float(np.median(ratio)),
            "q95": float(np.quantile(ratio, 0.95)),
        },
    }
    args.metadata.write_text(json.dumps(metadata, indent=2) + "\n")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--parameter-dir", type=Path, required=True)
    parser.add_argument("--chunk-output", type=Path, required=True)
    parser.add_argument("--true-output", type=Path, required=True)
    parser.add_argument("--polegon-output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--draw", type=int, default=50)
    parser.add_argument("--expected-chunks", type=int, default=50)
    parser.add_argument("--mutation-rate", type=float, default=1e-8)
    parser.add_argument("--polegon-binary", type=Path, default=Path(
        "/home/jri/src/singer-snakemake/resources/"
        "polegon-0.1.3-alpha-linux-x86_64/polegon"
    ))
    return parser.parse_args(argv)


if __name__ == "__main__":
    run(parse_args())
