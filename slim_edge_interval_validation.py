#!/usr/bin/env python3
"""Infer historical sample ages from independent SLiM mutation edges."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import tskit

from direct_frequency_age_infer import discrete_quantile, normalize_log_likelihood
from posterior_sample_age_infer import load_table, phi_lookup


def read_samples(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 33:
        raise ValueError(f"expected 33 retained haplotypes in {path}, found {len(rows)}")
    return rows


def node_for_row(ts: tskit.TreeSequence, row: dict) -> int:
    pedigree_id = int(row["pedigree_id"])
    haplotype_index = int(row["haplotype_index"])
    matches = [ind for ind in ts.individuals()
               if ind.metadata.get("pedigree_id") == pedigree_id]
    if len(matches) != 1:
        raise ValueError(f"pedigree {pedigree_id} has {len(matches)} matches")
    return int(matches[0].nodes[haplotype_index])


def compact(args):
    ts = tskit.load(args.input_tree)
    rows = read_samples(args.sample_file)
    # The production ARG is inferred from modern samples only. Historical
    # genotypes are recorded in the manifest, but must not split its edges.
    nodes = [node_for_row(ts, row) for row in rows if row["role"] == "modern"]
    if len(nodes) != 26:
        raise ValueError(f"expected 26 modern haplotypes, found {len(nodes)}")
    compacted = ts.simplify(nodes, keep_input_roots=True)
    args.output_tree.parent.mkdir(parents=True, exist_ok=True)
    compacted.dump(args.output_tree)
    print(json.dumps({"input_nodes": ts.num_nodes, "output_nodes": compacted.num_nodes,
                      "output_edges": compacted.num_edges,
                      "output_mutations": compacted.num_mutations}))


def mutation_slim_ids(mutation) -> set[str]:
    # SLiM's tskit allele state is a comma-separated mutation stack.  For the
    # mutation-table row at which a new lineage arose, its SLiM mutation ID is
    # the first element; older stacked IDs follow it.
    return {mutation.derived_state.split(",", 1)[0]}


def infer(args):
    tab = load_table(args.frequency_table)
    grid = np.asarray(tab["Tgrid"], dtype=float)
    ll = np.zeros((7, len(grid)), dtype=float)
    counts = {"examined": 0, "no_panel_polymorphic_lineage": 0,
              "ambiguous_focal_lineage": 0, "invalid_edge": 0,
              "nonfinite_frequency": 0, "retained": 0}
    diagnostic_rows = []

    for replicate in range(args.start, args.stop):
        stem = f"replicate_{replicate:06d}.seed_{args.base_seed + replicate}"
        tree_path = args.input_dir / "trees" / f"{stem}.trees"
        sample_path = args.input_dir / "samples" / f"{stem}.samples.tsv"
        if not tree_path.is_file() or not sample_path.is_file():
            raise ValueError(f"missing input for replicate {replicate}")
        counts["examined"] += 1
        ts = tskit.load(tree_path)
        rows = read_samples(sample_path)
        ancient = [row for row in rows if row["role"] == "ancient"]
        modern = [row for row in rows if row["role"] == "modern"]
        if len(ancient) != 7 or len(modern) != 26:
            raise ValueError(f"wrong sample roles in {sample_path}")

        modern_ids = [row["mutation_id"] for row in modern]
        allele_ids, allele_counts = np.unique(
            [x for x in modern_ids if x != "NA"], return_counts=True)
        candidates = [(mid, int(n)) for mid, n in zip(allele_ids, allele_counts)
                      if 0 < n < len(modern)]
        if not candidates:
            counts["no_panel_polymorphic_lineage"] += 1
            continue
        # A recurrent-mutation replicate can leave several panel-polymorphic
        # alleles. Choose one reproducibly without favoring low (older) SLiM
        # mutation IDs.
        rng = np.random.default_rng(
            np.uint64(args.base_seed + replicate) ^ np.uint64(0x5EED_A6E)
        )
        focal_id, d0 = candidates[int(rng.integers(len(candidates)))]
        focal_mutations = [m for m in ts.mutations() if focal_id in mutation_slim_ids(m)]
        if len(focal_mutations) != 1:
            counts["ambiguous_focal_lineage"] += 1
            continue
        mutation = focal_mutations[0]
        tree = ts.at(ts.site(mutation.site).position)
        parent = tree.parent(mutation.node)
        if parent == tskit.NULL:
            counts["invalid_edge"] += 1
            continue
        t_lo = float(ts.node(mutation.node).time)
        t_hi = float(ts.node(parent).time)
        if not t_hi > t_lo:
            counts["invalid_edge"] += 1
            continue
        lookup_lo, lookup_hi = (float(mutation.time), float(mutation.time)) \
            if args.mutation_age_source == "exact" else (t_lo, t_hi)
        phi = phi_lookup(tab, d0, lookup_lo, lookup_hi,
                         n_called=26, marginalise="uniform")
        if phi is None or not np.all(np.isfinite(phi)):
            counts["nonfinite_frequency"] += 1
            continue
        ancient_genotypes = np.asarray(
            [int(row["mutation_id"] == focal_id) for row in ancient], dtype=np.int8)
        with np.errstate(divide="ignore"):
            for j, genotype in enumerate(ancient_genotypes):
                ll[j] += np.log(phi) if genotype else np.log1p(-phi)
        diagnostic_rows.append((replicate, focal_id, d0, t_lo, t_hi,
                                *(int(x) for x in ancient_genotypes)))
        counts["retained"] += 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / f"chunk_{args.chunk_id:04d}.npz"
    with output.open("wb") as handle:
        np.savez_compressed(handle, log_likelihood=ll, candidate_ages=grid,
                            counts=json.dumps(counts), diagnostics=np.asarray(
                                diagnostic_rows, dtype=object))
    print(json.dumps(counts, indent=2))


def merge(args):
    paths = [args.chunk_dir / f"chunk_{i:04d}.npz" for i in range(args.n_chunks)]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise ValueError(f"missing {len(missing)} chunks; first {missing[0]}")
    ll = None
    counts = {}
    grid = None
    for path in paths:
        with np.load(path, allow_pickle=True) as data:
            if ll is None:
                ll = np.asarray(data["log_likelihood"], dtype=float)
                grid = np.asarray(data["candidate_ages"], dtype=float)
            else:
                ll += data["log_likelihood"]
            for key, value in json.loads(str(data["counts"])).items():
                counts[key] = counts.get(key, 0) + value
    posterior = normalize_log_likelihood(ll)
    true_ages = np.asarray([6000, 5000, 4000, 3000, 2000, 1000, 500], dtype=float)
    rows = []
    for j, truth in enumerate(true_ages):
        p = posterior[j]
        rows.append({"true_age": truth, "map_age": float(grid[np.argmax(p)]),
                     "posterior_mean_age": float(np.sum(grid * p)),
                     "ci95_lower": discrete_quantile(grid, p, 0.025),
                     "ci95_upper": discrete_quantile(grid, p, 0.975)})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "age_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    result = {"model": "independent SLiM loci; " + (
                  "exact SLiM mutation time" if args.mutation_age_source == "exact"
                  else "true mutation edge; uniform edge marginalisation"),
              "epsilon": 0.0, "counts": counts,
              "map_bias_estimated_minus_true": float(np.mean(
                  [row["map_age"] - row["true_age"] for row in rows])),
              "all_true_ages_in_ci95": bool(all(
                  row["ci95_lower"] <= row["true_age"] <= row["ci95_upper"]
                  for row in rows))}
    (args.output_dir / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    np.savez_compressed(args.output_dir / "age_posteriors.npz", candidate_ages=grid,
                        posterior=posterior, log_likelihood=ll, true_ages=true_ages)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(8, 8), constrained_layout=True)
    for j, truth in enumerate(true_ages):
        axes[0].plot(grid, posterior[j], label=f"true age {truth:g}")
        axes[0].axvline(truth, color=f"C{j}", lw=0.7, alpha=0.45)
    axes[0].set(xlabel="Candidate sample age (generations)",
                ylabel="Posterior density", xlim=(0, 7000))
    axes[0].legend(fontsize=8, ncol=2)
    maps = np.asarray([row["map_age"] for row in rows])
    axes[1].plot([0, 7000], [0, 7000], color="0.5", ls="--", lw=1)
    axes[1].scatter(true_ages, maps, s=35)
    axes[1].set(xlabel="True sample age (generations)",
                ylabel="Posterior MAP age (generations)",
                xlim=(0, 7000), ylim=(0, 7000))
    fig.savefig(args.output_dir / "age_validation.png", dpi=180)
    plt.close(fig)
    print(json.dumps(result, indent=2))


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    c = sub.add_parser("compact")
    c.add_argument("--input-tree", type=Path, required=True)
    c.add_argument("--sample-file", type=Path, required=True)
    c.add_argument("--output-tree", type=Path, required=True)
    c.set_defaults(func=compact)
    i = sub.add_parser("infer")
    i.add_argument("--input-dir", type=Path, required=True)
    i.add_argument("--frequency-table", type=Path, required=True)
    i.add_argument("--output-dir", type=Path, required=True)
    i.add_argument("--chunk-id", type=int, required=True)
    i.add_argument("--start", type=int, required=True)
    i.add_argument("--stop", type=int, required=True)
    i.add_argument("--base-seed", type=int, default=202609090000)
    i.add_argument("--mutation-age-source", choices=("edge", "exact"), default="edge")
    i.set_defaults(func=infer)
    m = sub.add_parser("merge")
    m.add_argument("--chunk-dir", type=Path, required=True)
    m.add_argument("--output-dir", type=Path, required=True)
    m.add_argument("--n-chunks", type=int, required=True)
    m.add_argument("--mutation-age-source", choices=("edge", "exact"), default="edge")
    m.set_defaults(func=merge)
    return p


if __name__ == "__main__":
    args = parser().parse_args()
    args.func(args)
