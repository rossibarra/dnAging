#!/usr/bin/env python3
"""Infer sampling ages using known simulated allele-frequency trajectories.

Each independent SLiM replicate is treated as one biallelic SNP.  All mutation
lineages at the one-base locus are collapsed into a single derived state, so
the known derived frequency is one minus the recorded ancestral frequency.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


FREQUENCY_SUFFIX = ".frequencies.tsv"
SAMPLE_SUFFIX = ".samples.tsv"


def replicate_stem(path: Path, suffix: str) -> str:
    if not path.name.endswith(suffix):
        raise ValueError(f"{path} does not end with {suffix}")
    return path.name[: -len(suffix)]


def discover_pairs(input_dir: Path) -> list[tuple[Path, Path]]:
    frequency_dir = input_dir / "frequencies"
    sample_dir = input_dir / "samples"
    sample_paths = sorted(sample_dir.glob(f"*{SAMPLE_SUFFIX}"))
    if not sample_paths:
        raise ValueError(f"no sample files found in {sample_dir}")

    pairs = []
    for sample_path in sample_paths:
        stem = replicate_stem(sample_path, SAMPLE_SUFFIX)
        frequency_path = frequency_dir / f"{stem}{FREQUENCY_SUFFIX}"
        if not frequency_path.is_file():
            raise ValueError(f"missing frequency file paired with {sample_path}")
        pairs.append((frequency_path, sample_path))
    return pairs


def read_sample_genotypes(path: Path) -> tuple[int, int, np.ndarray, np.ndarray]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no genotype rows in {path}")

    replicates = {int(row["replicate"]) for row in rows}
    seeds = {int(row["seed"]) for row in rows}
    if len(replicates) != 1 or len(seeds) != 1:
        raise ValueError(f"replicate or seed changes within {path}")

    generations = np.asarray([int(row["generation"]) for row in rows], dtype=np.int64)
    genotypes = np.asarray([int(row["haploid_genotype"]) for row in rows], dtype=np.int8)
    if np.any((genotypes != 0) & (genotypes != 1)):
        raise ValueError(f"non-binary haploid genotype in {path}")
    if len(np.unique(generations)) != len(generations):
        raise ValueError(f"duplicate sampling generation in {path}")
    order = np.argsort(generations)
    return replicates.pop(), seeds.pop(), generations[order], genotypes[order]


def read_derived_frequency(path: Path, end_generation: int) -> np.ndarray:
    """Read p(derived) from the ancestral row at every generation."""
    derived = np.full(end_generation, np.nan, dtype=np.float64)
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        required = {"generation", "allele_state", "population_frequency"}
        if not required.issubset(header):
            raise ValueError(f"required columns absent from {path}")
        generation_i = header.index("generation")
        state_i = header.index("allele_state")
        frequency_i = header.index("population_frequency")

        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if fields[state_i] != "ancestral":
                continue
            generation = int(fields[generation_i])
            if generation < 1 or generation > end_generation:
                raise ValueError(
                    f"generation {generation} outside 1..{end_generation} "
                    f"in {path}:{line_number}"
                )
            if not math.isnan(derived[generation - 1]):
                raise ValueError(f"duplicate ancestral row in {path}:{line_number}")
            ancestral = float(fields[frequency_i])
            if ancestral < -1e-10 or ancestral > 1.0 + 1e-10:
                raise ValueError(f"invalid ancestral frequency in {path}:{line_number}")
            derived[generation - 1] = min(1.0, max(0.0, 1.0 - ancestral))

    missing = np.flatnonzero(np.isnan(derived))
    if len(missing):
        preview = ",".join(str(x + 1) for x in missing[:5])
        raise ValueError(f"missing ancestral rows in {path}; first generations: {preview}")
    return derived


def add_bernoulli_log_likelihood(
    log_likelihood: np.ndarray, derived_frequency: np.ndarray, genotypes: np.ndarray
) -> None:
    with np.errstate(divide="ignore", invalid="ignore"):
        log_derived = np.log(derived_frequency)
        log_ancestral = np.log1p(-derived_frequency)
    for sample_i, genotype in enumerate(genotypes):
        log_likelihood[sample_i] += log_derived if genotype == 1 else log_ancestral


def run_chunk(args: argparse.Namespace) -> None:
    pairs = discover_pairs(args.input_dir)
    n_replicates = len(pairs)
    start = args.chunk_id * n_replicates // args.n_chunks
    stop = (args.chunk_id + 1) * n_replicates // args.n_chunks
    selected = pairs[start:stop]
    if not selected:
        raise ValueError(
            f"chunk {args.chunk_id} is empty; {n_replicates} replicates "
            f"cannot fill {args.n_chunks} chunks"
        )

    sample_generations = None
    log_likelihood = None
    replicate_ids = []
    seeds = []
    for pair_i, (frequency_path, sample_path) in enumerate(selected, start=1):
        replicate, seed, generations, genotypes = read_sample_genotypes(sample_path)
        if sample_generations is None:
            sample_generations = generations
            log_likelihood = np.zeros(
                (len(generations), args.end_generation), dtype=np.float64
            )
        elif not np.array_equal(generations, sample_generations):
            raise ValueError(f"sampling generations differ in {sample_path}")

        derived_frequency = read_derived_frequency(frequency_path, args.end_generation)
        add_bernoulli_log_likelihood(log_likelihood, derived_frequency, genotypes)
        replicate_ids.append(replicate)
        seeds.append(seed)
        if pair_i % 10 == 0 or pair_i == len(selected):
            print(
                f"chunk {args.chunk_id}: {pair_i}/{len(selected)} replicates",
                flush=True,
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / f"chunk_{args.chunk_id:04d}.npz"
    tmp = output.with_suffix(f".npz.tmp.{args.chunk_id}")
    with tmp.open("wb") as handle:
        np.savez_compressed(
            handle,
            log_likelihood=log_likelihood,
            sample_generations=sample_generations,
            replicate_ids=np.asarray(replicate_ids, dtype=np.int64),
            seeds=np.asarray(seeds, dtype=np.int64),
            end_generation=np.asarray(args.end_generation, dtype=np.int64),
            n_total_replicates=np.asarray(n_replicates, dtype=np.int64),
            chunk_id=np.asarray(args.chunk_id, dtype=np.int64),
            n_chunks=np.asarray(args.n_chunks, dtype=np.int64),
        )
    tmp.replace(output)
    print(f"wrote {output}")


def normalize_log_likelihood(log_likelihood: np.ndarray) -> np.ndarray:
    posterior = np.zeros_like(log_likelihood, dtype=np.float64)
    for row_i, row in enumerate(log_likelihood):
        finite = np.isfinite(row)
        if not np.any(finite):
            raise ValueError(f"all candidate ages have zero likelihood for sample {row_i}")
        maximum = np.max(row[finite])
        weights = np.zeros_like(row)
        weights[finite] = np.exp(row[finite] - maximum)
        posterior[row_i] = weights / weights.sum()
    return posterior


def discrete_quantile(values: np.ndarray, probabilities: np.ndarray, q: float) -> float:
    index = int(np.searchsorted(np.cumsum(probabilities), q, side="left"))
    return float(values[min(index, len(values) - 1)])


def write_results(
    output_dir: Path,
    sample_generations: np.ndarray,
    end_generation: int,
    log_likelihood: np.ndarray,
    posterior_by_generation: np.ndarray,
    n_replicates: int,
) -> None:
    candidate_generations = np.arange(1, end_generation + 1, dtype=np.int64)
    candidate_ages_descending = end_generation - candidate_generations
    age_order = np.argsort(candidate_ages_descending)
    candidate_ages = candidate_ages_descending[age_order]
    posterior_by_age = posterior_by_generation[:, age_order]
    log_likelihood_by_age = log_likelihood[:, age_order]

    summaries = []
    for sample_i, sample_generation in enumerate(sample_generations):
        probabilities = posterior_by_age[sample_i]
        true_age = int(end_generation - sample_generation)
        map_age = int(candidate_ages[np.argmax(probabilities)])
        mean_age = float(np.sum(candidate_ages * probabilities))
        lower = discrete_quantile(candidate_ages, probabilities, 0.025)
        upper = discrete_quantile(candidate_ages, probabilities, 0.975)
        true_index = int(np.searchsorted(candidate_ages, true_age))
        pit = float(np.sum(probabilities[: true_index + 1]))
        summaries.append(
            {
                "sample_generation": int(sample_generation),
                "true_age": true_age,
                "map_age": map_age,
                "posterior_mean_age": mean_age,
                "ci95_lower": lower,
                "ci95_upper": upper,
                "ci95_contains_true_age": int(lower <= true_age <= upper),
                "posterior_cdf_at_true_age": pit,
            }
        )

    summary_path = output_dir / "age_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summaries[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summaries)

    posterior_path = output_dir / "age_posteriors.tsv"
    with posterior_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample_generation",
                "true_age",
                "candidate_generation",
                "candidate_age",
                "log_likelihood",
                "posterior_probability",
            ]
        )
        for sample_i, sample_generation in enumerate(sample_generations):
            true_age = int(end_generation - sample_generation)
            for age_i, candidate_age in enumerate(candidate_ages):
                writer.writerow(
                    [
                        int(sample_generation),
                        true_age,
                        int(end_generation - candidate_age),
                        int(candidate_age),
                        f"{log_likelihood_by_age[sample_i, age_i]:.17g}",
                        f"{posterior_by_age[sample_i, age_i]:.17g}",
                    ]
                )

    np.savez_compressed(
        output_dir / "age_posteriors.npz",
        sample_generations=sample_generations,
        true_ages=end_generation - sample_generations,
        candidate_ages=candidate_ages,
        log_likelihood=log_likelihood_by_age,
        posterior=posterior_by_age,
        n_replicates=np.asarray(n_replicates, dtype=np.int64),
    )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(
        len(sample_generations), 1, figsize=(9, 1.7 * len(sample_generations)), sharex=True
    )
    axes = np.atleast_1d(axes)
    for sample_i, (axis, summary) in enumerate(zip(axes, summaries)):
        axis.plot(candidate_ages, posterior_by_age[sample_i], linewidth=1.2)
        axis.axvline(summary["true_age"], color="black", linestyle="--", linewidth=0.9)
        axis.set_ylabel(f"gen {summary['sample_generation']}")
    axes[-1].set_xlabel("Sample age (generations before generation 40,000)")
    axes[-1].set_xlim(0, 7000)
    figure.supylabel("Posterior probability")
    figure.tight_layout()
    figure.savefig(output_dir / "age_posteriors.svg")
    figure.savefig(output_dir / "age_posteriors.png", dpi=180)
    plt.close(figure)

    true_ages = np.asarray([row["true_age"] for row in summaries])
    map_ages = np.asarray([row["map_age"] for row in summaries])
    lower = np.asarray([row["ci95_lower"] for row in summaries])
    upper = np.asarray([row["ci95_upper"] for row in summaries])
    figure, axis = plt.subplots(figsize=(6, 6))
    limits = [0, max(true_ages.max(), upper.max()) * 1.03]
    axis.plot(limits, limits, color="0.5", linestyle="--", linewidth=1)
    axis.vlines(true_ages, lower, upper, linewidth=1.2)
    axis.scatter(true_ages, map_ages, zorder=3)
    axis.set(xlabel="True age", ylabel="Posterior MAP age", xlim=limits, ylim=limits)
    figure.tight_layout()
    figure.savefig(output_dir / "age_calibration.svg")
    figure.savefig(output_dir / "age_calibration.png", dpi=180)
    plt.close(figure)

    with (output_dir / "run.json").open("w") as handle:
        json.dump(
            {
                "model": "known-frequency independent-SNP Bernoulli likelihood",
                "n_replicates": n_replicates,
                "n_samples": len(sample_generations),
                "sample_generations": sample_generations.tolist(),
                "true_ages": (end_generation - sample_generations).tolist(),
                "end_generation": end_generation,
                "prior": "discrete uniform over candidate generations 1..END_GEN",
                "derived_frequency": "1 - recorded ancestral frequency",
                "genotyping_error": 0.0,
            },
            handle,
            indent=2,
        )


def merge_chunks(args: argparse.Namespace) -> None:
    chunk_paths = [args.chunk_dir / f"chunk_{i:04d}.npz" for i in range(args.n_chunks)]
    missing = [str(path) for path in chunk_paths if not path.is_file()]
    if missing:
        raise ValueError(f"missing {len(missing)} chunk files; first: {missing[0]}")

    log_likelihood = None
    sample_generations = None
    end_generation = None
    replicate_ids = []
    seeds = []
    expected_total = None
    for expected_chunk_id, path in enumerate(chunk_paths):
        with np.load(path) as chunk:
            chunk_id = int(chunk["chunk_id"])
            if chunk_id != expected_chunk_id or int(chunk["n_chunks"]) != args.n_chunks:
                raise ValueError(f"chunk metadata mismatch in {path}")
            generations = chunk["sample_generations"]
            end = int(chunk["end_generation"])
            total = int(chunk["n_total_replicates"])
            if sample_generations is None:
                sample_generations = generations.copy()
                end_generation = end
                expected_total = total
                log_likelihood = np.zeros_like(chunk["log_likelihood"])
            elif (
                not np.array_equal(generations, sample_generations)
                or end != end_generation
                or total != expected_total
            ):
                raise ValueError(f"incompatible chunk metadata in {path}")
            log_likelihood += chunk["log_likelihood"]
            replicate_ids.extend(chunk["replicate_ids"].tolist())
            seeds.extend(chunk["seeds"].tolist())

    if len(replicate_ids) != expected_total:
        raise ValueError(f"expected {expected_total} replicates, merged {len(replicate_ids)}")
    if len(set(replicate_ids)) != len(replicate_ids):
        raise ValueError("duplicate replicate IDs across chunks")
    if len(set(seeds)) != len(seeds):
        raise ValueError("duplicate random seeds across chunks")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    posterior = normalize_log_likelihood(log_likelihood)
    write_results(
        args.output_dir,
        sample_generations,
        end_generation,
        log_likelihood,
        posterior,
        len(replicate_ids),
    )
    print(f"merged {len(replicate_ids)} independent SNPs into {args.output_dir}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    chunk = subparsers.add_parser("chunk", help="calculate a partial log likelihood")
    chunk.add_argument("--input-dir", type=Path, required=True)
    chunk.add_argument("--output-dir", type=Path, required=True)
    chunk.add_argument("--chunk-id", type=int, required=True)
    chunk.add_argument("--n-chunks", type=int, required=True)
    chunk.add_argument("--end-generation", type=int, default=40000)
    chunk.set_defaults(func=run_chunk)

    merge = subparsers.add_parser("merge", help="sum chunks and write posteriors")
    merge.add_argument("--chunk-dir", type=Path, required=True)
    merge.add_argument("--output-dir", type=Path, required=True)
    merge.add_argument("--n-chunks", type=int, required=True)
    merge.set_defaults(func=merge_chunks)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if hasattr(args, "chunk_id") and not 0 <= args.chunk_id < args.n_chunks:
        raise SystemExit("--chunk-id must be in [0, --n-chunks)")
    args.func(args)


if __name__ == "__main__":
    main()
