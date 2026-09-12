#!/usr/bin/env python3
"""Validate diffusion frequency estimates against one-site SLiM simulations.

Test A conditions on a mutation's exact SLiM origin age and a reproducibly
sampled present-day panel count.  It then replaces the known trajectory in the
age likelihood with E[p_T | d0, mutation age] from the production diffusion
table.  One replicate contributes at most one SNP.  When recurrent mutation
leaves multiple alleles, the modern panel is drawn jointly across all alleles
and one panel-polymorphic mutation lineage is selected reproducibly.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from direct_frequency_age_infer import (
    FREQUENCY_SUFFIX,
    SAMPLE_SUFFIX,
    discover_pairs,
    discrete_quantile,
    normalize_log_likelihood,
)
from posterior_sample_age_infer import load_table, phi_lookup


PANEL_RNG_XOR = 0x5EED_A6E


def read_sample_mutation_ids(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no genotype rows in {path}")
    replicates = {int(r["replicate"]) for r in rows}
    seeds = {int(r["seed"]) for r in rows}
    if len(replicates) != 1 or len(seeds) != 1:
        raise ValueError(f"replicate or seed changes within {path}")
    generations = np.asarray([int(r["generation"]) for r in rows], dtype=np.int64)
    mutation_ids = np.asarray([r["mutation_id"] for r in rows], dtype=str)
    if len(np.unique(generations)) != len(generations):
        raise ValueError(f"duplicate sampling generation in {path}")
    order = np.argsort(generations)
    return replicates.pop(), seeds.pop(), generations[order], mutation_ids[order]


def read_final_trajectories(path: Path, generations: np.ndarray, end_generation: int):
    """Return every final mutation and its frequencies at requested times.

    The full file is streamed once.  Mutation IDs are strings because SLiM's
    IDs need not fit assumptions made about a particular integer width.
    """
    wanted = {int(x) for x in generations}
    wanted.add(int(end_generation))
    by_generation = {g: {} for g in wanted}
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        ci = {name: i for i, name in enumerate(header)}
        required = {
            "generation", "allele_state", "mutation_id",
            "origin_generation", "population_frequency",
        }
        if not required.issubset(ci):
            raise ValueError(f"required columns absent from {path}")
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if fields[ci["allele_state"]] != "derived":
                continue
            generation = int(fields[ci["generation"]])
            if generation not in wanted:
                continue
            mutation_id = fields[ci["mutation_id"]]
            frequency = float(fields[ci["population_frequency"]])
            origin = int(fields[ci["origin_generation"]])
            if mutation_id in by_generation[generation]:
                raise ValueError(f"duplicate mutation {mutation_id} in {path}:{line_number}")
            by_generation[generation][mutation_id] = (frequency, origin)

    trajectories = []
    for mutation_id, (present_frequency, origin_generation) in sorted(
        by_generation[end_generation].items()
    ):
        if present_frequency <= 0.0:
            continue
        true_frequency = np.asarray(
            [by_generation[int(g)].get(
                mutation_id, (0.0, origin_generation)
            )[0] for g in generations],
            dtype=np.float64,
        )
        trajectories.append(
            (mutation_id, origin_generation, present_frequency, true_frequency)
        )
    return trajectories


def sample_panel_focal(seed: int, trajectories, n_panel: int):
    """Jointly draw the modern panel and select at most one observed mutation.

    A multinomial draw is required when recurrent mutation leaves multiple
    mutually exclusive alleles at the one-base locus.  Choosing one candidate
    only after that draw mirrors ascertainment of a panel-polymorphic SNP and
    guarantees that a replicate contributes no more than one likelihood term.
    """
    if not trajectories:
        return None
    derived = np.asarray([x[2] for x in trajectories], dtype=np.float64)
    ancestral = max(0.0, 1.0 - float(derived.sum()))
    probabilities = np.concatenate(([ancestral], derived))
    total = probabilities.sum()
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("invalid final allele frequencies")
    probabilities /= total
    rng = np.random.default_rng(np.uint64(seed) ^ np.uint64(PANEL_RNG_XOR))
    counts = rng.multinomial(n_panel, probabilities)[1:]
    candidates = np.flatnonzero((counts > 0) & (counts < n_panel))
    if not len(candidates):
        return None
    chosen = int(rng.choice(candidates))
    return trajectories[chosen], int(counts[chosen])


def run_chunk(args: argparse.Namespace) -> None:
    pairs = discover_pairs(args.input_dir)
    start = args.chunk_id * len(pairs) // args.n_chunks
    stop = (args.chunk_id + 1) * len(pairs) // args.n_chunks
    pairs = pairs[start:stop]
    if not pairs:
        raise ValueError(f"empty chunk {args.chunk_id}")

    tab = load_table(args.frequency_table)
    if tab["n_sample"] != args.n_panel:
        raise ValueError(
            f"frequency table has n_sample={tab['n_sample']}, expected {args.n_panel}"
        )
    candidate_ages = np.asarray(tab["Tgrid"], dtype=np.float64)
    log_likelihood = None
    diagnostics = []
    counts = {"examined": 0, "no_final_lineage": 0,
              "no_panel_polymorphic_lineage": 0,
              "mutation_age_outside_table": 0,
              "nonfinite_estimate": 0, "retained": 0}
    sample_generations_ref = None

    for pair_i, (frequency_path, sample_path) in enumerate(pairs, start=1):
        counts["examined"] += 1
        replicate, seed, sample_generations, sampled_mutation_ids = \
            read_sample_mutation_ids(sample_path)
        if sample_generations_ref is None:
            sample_generations_ref = sample_generations
            log_likelihood = np.zeros(
                (len(sample_generations), len(candidate_ages)), dtype=np.float64
            )
        elif not np.array_equal(sample_generations, sample_generations_ref):
            raise ValueError(f"sampling generations differ in {sample_path}")

        trajectories = read_final_trajectories(
            frequency_path, sample_generations, args.end_generation
        )
        if not trajectories:
            counts["no_final_lineage"] += 1
            continue
        sampled_focal = sample_panel_focal(seed, trajectories, args.n_panel)
        if sampled_focal is None:
            counts["no_panel_polymorphic_lineage"] += 1
            continue
        focal, d0 = sampled_focal
        focal_id, origin_generation, present_frequency, true_frequency = focal
        mutation_age = args.end_generation - origin_generation
        if mutation_age < float(tab["age"][0]) or mutation_age > float(tab["age"][-1]):
            counts["mutation_age_outside_table"] += 1
            continue
        estimated = phi_lookup(
            tab, d0, mutation_age, mutation_age,
            n_called=args.n_panel, marginalise="uniform",
        )
        if estimated is None or not np.all(np.isfinite(estimated)):
            counts["nonfinite_estimate"] += 1
            continue
        genotypes = (sampled_mutation_ids == focal_id).astype(np.int8)
        with np.errstate(divide="ignore"):
            for sample_i, genotype in enumerate(genotypes):
                log_likelihood[sample_i] += (
                    np.log(estimated) if genotype else np.log1p(-estimated)
                )

        estimated_at_truth = np.interp(
            args.end_generation - sample_generations,
            candidate_ages,
            estimated,
        )
        for generation, genotype, p_true, p_est in zip(
            sample_generations, genotypes, true_frequency, estimated_at_truth
        ):
            diagnostics.append((
                replicate, seed, focal_id, origin_generation, mutation_age,
                present_frequency, args.n_panel, d0, int(generation),
                args.end_generation - int(generation), int(genotype), p_true, p_est,
            ))
        counts["retained"] += 1
        if pair_i % 10 == 0 or pair_i == len(pairs):
            print(f"chunk {args.chunk_id}: {pair_i}/{len(pairs)} files, "
                  f"{counts['retained']} retained", flush=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / f"chunk_{args.chunk_id:04d}.npz"
    tmp = output.with_suffix(f".npz.tmp.{args.chunk_id}")
    diagnostic_dtype = np.dtype([
        ("replicate", "i8"), ("seed", "i8"), ("mutation_id", "U32"),
        ("origin_generation", "i8"), ("mutation_age", "i8"),
        ("present_frequency", "f8"), ("n_panel", "i4"), ("d0", "i4"),
        ("sample_generation", "i8"), ("true_sample_age", "i8"),
        ("haploid_genotype", "i1"), ("true_frequency", "f8"),
        ("estimated_frequency", "f8"),
    ])
    with tmp.open("wb") as handle:
        np.savez_compressed(
            handle,
            log_likelihood=log_likelihood,
            sample_generations=sample_generations_ref,
            candidate_ages=candidate_ages,
            diagnostics=np.asarray(diagnostics, dtype=diagnostic_dtype),
            count_names=np.asarray(list(counts), dtype="U40"),
            count_values=np.asarray(list(counts.values()), dtype=np.int64),
            chunk_id=np.asarray(args.chunk_id), n_chunks=np.asarray(args.n_chunks),
            n_total_replicates=np.asarray(len(discover_pairs(args.input_dir))),
            n_panel=np.asarray(args.n_panel), end_generation=np.asarray(args.end_generation),
        )
    tmp.replace(output)
    print(f"wrote {output}")


def _summaries(sample_generations, ages, log_likelihood, end_generation):
    posterior = normalize_log_likelihood(log_likelihood)
    rows = []
    for i, generation in enumerate(sample_generations):
        p = posterior[i]
        true_age = int(end_generation - generation)
        rows.append({
            "sample_generation": int(generation), "true_age": true_age,
            "map_age": float(ages[np.argmax(p)]),
            "posterior_mean_age": float(np.sum(ages * p)),
            "ci95_lower": discrete_quantile(ages, p, 0.025),
            "ci95_upper": discrete_quantile(ages, p, 0.975),
        })
    return posterior, rows


def merge_chunks(args: argparse.Namespace) -> None:
    paths = [args.chunk_dir / f"chunk_{i:04d}.npz" for i in range(args.n_chunks)]
    missing = [p for p in paths if not p.is_file()]
    if missing:
        raise ValueError(f"missing {len(missing)} chunks; first: {missing[0]}")
    log_likelihood = None
    diagnostics, total_counts = [], {}
    sample_generations = candidate_ages = None
    metadata = None
    for expected, path in enumerate(paths):
        with np.load(path) as z:
            if int(z["chunk_id"]) != expected or int(z["n_chunks"]) != args.n_chunks:
                raise ValueError(f"chunk metadata mismatch in {path}")
            current = (z["sample_generations"], z["candidate_ages"],
                       int(z["n_total_replicates"]), int(z["n_panel"]),
                       int(z["end_generation"]))
            if metadata is None:
                sample_generations, candidate_ages = current[:2]
                metadata = current[2:]
                log_likelihood = np.zeros_like(z["log_likelihood"])
            elif (not np.array_equal(current[0], sample_generations)
                  or not np.array_equal(current[1], candidate_ages)
                  or current[2:] != metadata):
                raise ValueError(f"incompatible metadata in {path}")
            log_likelihood += z["log_likelihood"]
            diagnostics.append(z["diagnostics"].copy())
            for name, value in zip(z["count_names"], z["count_values"]):
                total_counts[str(name)] = total_counts.get(str(name), 0) + int(value)
    if total_counts["examined"] != metadata[0]:
        raise ValueError("chunks did not examine every replicate exactly once")
    diagnostic = np.concatenate(diagnostics)
    if len(diagnostic) != total_counts["retained"] * len(sample_generations):
        raise ValueError("diagnostic row count disagrees with retained-site count")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    posterior, summaries = _summaries(
        sample_generations, candidate_ages, log_likelihood, metadata[2]
    )
    for row in summaries:
        row["ci95_contains_true_age"] = int(
            row["ci95_lower"] <= row["true_age"] <= row["ci95_upper"]
        )
    with (args.output_dir / "age_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader(); writer.writerows(summaries)
    with (args.output_dir / "frequency_diagnostics.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(diagnostic.dtype.names)
        writer.writerows(diagnostic.tolist())
    np.savez_compressed(
        args.output_dir / "age_posteriors.npz", candidate_ages=candidate_ages,
        sample_generations=sample_generations, posterior=posterior,
        log_likelihood=log_likelihood,
    )

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(sample_generations), 1,
                             figsize=(9, 1.7 * len(sample_generations)), sharex=True)
    for axis, row, p in zip(np.atleast_1d(axes), summaries, posterior):
        axis.plot(candidate_ages, p, lw=1.2)
        axis.axvline(row["true_age"], color="black", ls="--", lw=0.9)
        axis.set_ylabel(f"gen {row['sample_generation']}")
    axes[-1].set_xlabel("Sample age (generations before generation 40,000)")
    axes[-1].set_xlim(candidate_ages[0], candidate_ages[-1])
    fig.supylabel("Posterior probability"); fig.tight_layout()
    fig.savefig(args.output_dir / "age_posteriors.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.hexbin(diagnostic["true_frequency"], diagnostic["estimated_frequency"],
              gridsize=45, mincnt=1, bins="log")
    ax.plot([0, 1], [0, 1], color="black", ls="--", lw=1)
    ax.set(xlabel="True SLiM frequency", ylabel="Diffusion estimated frequency",
           xlim=(0, 1), ylim=(0, 1))
    fig.tight_layout(); fig.savefig(args.output_dir / "frequency_comparison.png", dpi=180)
    plt.close(fig)

    diff = diagnostic["estimated_frequency"] - diagnostic["true_frequency"]
    run = {
        "model": "Test A: exact mutation age + observed modern count -> diffusion E[p_T]",
        "n_input_replicates": metadata[0], "n_retained_snps": total_counts["retained"],
        "n_panel": metadata[1], "end_generation": metadata[2],
        "selection_counts": total_counts,
        "panel_count_rng": f"numpy default_rng(seed XOR {PANEL_RNG_XOR})",
        "frequency_bias_estimated_minus_true": float(np.mean(diff)),
        "frequency_mae": float(np.mean(np.abs(diff))),
        "frequency_rmse": float(np.sqrt(np.mean(diff * diff))),
        "ci95_coverage_count": int(sum(r["ci95_contains_true_age"] for r in summaries)),
        "ci95_coverage_total": len(summaries),
        "recurrent_mutation_handling": (
            "The modern panel is sampled jointly from ancestral plus every extant "
            "mutation lineage. One panel-polymorphic mutation is chosen uniformly, "
            "and that focal lineage is contrasted against every other allelic state. "
            "Each replicate contributes at most one SNP."
        ),
    }
    with (args.output_dir / "run.json").open("w") as handle:
        json.dump(run, handle, indent=2)
    print(json.dumps(run, indent=2))


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    chunk = sub.add_parser("chunk")
    chunk.add_argument("--input-dir", type=Path, required=True)
    chunk.add_argument("--frequency-table", type=Path, required=True)
    chunk.add_argument("--output-dir", type=Path, required=True)
    chunk.add_argument("--chunk-id", type=int, required=True)
    chunk.add_argument("--n-chunks", type=int, required=True)
    chunk.add_argument("--n-panel", type=int, default=26)
    chunk.add_argument("--end-generation", type=int, default=40000)
    chunk.set_defaults(func=run_chunk)
    merge = sub.add_parser("merge")
    merge.add_argument("--chunk-dir", type=Path, required=True)
    merge.add_argument("--output-dir", type=Path, required=True)
    merge.add_argument("--n-chunks", type=int, required=True)
    merge.set_defaults(func=merge_chunks)
    return parser


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
