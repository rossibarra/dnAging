#!/usr/bin/env python3
"""Cross-validate chromosome-specific scale/epsilon calibration before age combination."""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import joint_scale_epsilon_calibration as joint
import radiocarbon_calibration as base


def combined_age_loglik(sample_ll, chromosome_log_weights):
    """Marginalize each chromosome's calibration, then combine one shared age."""
    sample_ll = np.asarray(sample_ll)
    weights = np.asarray(chromosome_log_weights)
    expected_weights = (sample_ll.shape[2], sample_ll.shape[0], sample_ll.shape[1])
    if sample_ll.ndim != 4 or weights.shape != expected_weights:
        raise ValueError("expected scale x epsilon x chromosome x age and matching weights")
    return np.sum([
        joint.marginal_age_loglik_weights(sample_ll[:, :, c, :], weights[c])
        for c in range(sample_ll.shape[2])
    ], axis=0)


def run(args):
    source = args.likelihood_dir
    age_grid = np.load(source / "grid.npy")
    scale_grid = np.load(source / "time_scale_grid.npy")
    epsilon_grid = np.load(source / "epsilon_grid.npy")
    ll = np.load(source / "ll_by_chrom_scale_epsilon.npy", mmap_mode="r")
    samples = (source / "samples.txt").read_text().split()
    metadata = json.loads((source / "run.json").read_text())
    chromosomes = metadata["chromosomes"]
    expected = (len(scale_grid), len(epsilon_grid), len(chromosomes), len(samples), len(age_grid))
    if ll.shape != expected:
        raise ValueError(f"likelihood shape {ll.shape}; expected {expected}")
    rc = base.read_radiocarbon(args.radiocarbon)
    dated = [sample for sample in samples if sample in rc]
    index = {sample: i for i, sample in enumerate(samples)}
    fold_map = base.automatic_folds(dated, rc, args.n_folds)
    folds = sorted(set(fold_map.values()))

    profiles = np.empty((len(chromosomes), len(dated), len(scale_grid), len(epsilon_grid)))
    for c in range(len(chromosomes)):
        for d, sample in enumerate(dated):
            profiles[c, d] = joint.joint_radiocarbon_log_profile(
                age_grid, ll[:, :, c, index[sample], :], *rc[sample]
            )

    all_log_weights = np.empty((len(chromosomes), len(scale_grid), len(epsilon_grid)))
    chromosome_summaries = {}
    for c, chromosome in enumerate(chromosomes):
        profile = np.sum(profiles[c], axis=0)
        all_log_weights[c] = joint.joint_log_posterior(profile, scale_grid, epsilon_grid)
        chromosome_summaries[chromosome] = joint.joint_summary(
            profile, scale_grid, epsilon_grid
        )

    rows = []
    dated_lookup = {sample: d for d, sample in enumerate(dated)}
    for fold in folds:
        training = [s for s in dated if fold_map[s] != fold]
        held_out = [s for s in dated if fold_map[s] == fold]
        fold_weights = np.empty_like(all_log_weights)
        for c in range(len(chromosomes)):
            training_profile = np.sum(
                [profiles[c, dated_lookup[s]] for s in training], axis=0
            )
            fold_weights[c] = joint.joint_log_posterior(
                training_profile, scale_grid, epsilon_grid
            )
        for sample in held_out:
            genetic = combined_age_loglik(ll[:, :, :, index[sample], :], fold_weights)
            estimate = base.summarize_density(age_grid, genetic)
            rc_age, rc_log_weight = rc[sample]
            rc_mean = float(np.sum(rc_age * np.exp(rc_log_weight)))
            rc_lo = base.weighted_quantile(rc_age, rc_log_weight, 0.025)
            rc_hi = base.weighted_quantile(rc_age, rc_log_weight, 0.975)
            row = {
                "sample": sample, "fold": fold, "n_training": len(training),
                "radiocarbon_mean_T": rc_mean,
                "radiocarbon_ci95_lower_T": rc_lo,
                "radiocarbon_ci95_upper_T": rc_hi,
                **estimate,
            }
            row["signed_error_T"] = row["map_T"] - rc_mean
            row["abs_error_T"] = abs(row["signed_error_T"])
            row["map_in_radiocarbon_95"] = rc_lo <= row["map_T"] <= rc_hi
            row["radiocarbon_mean_in_genetic_95"] = (
                row["ci95_lower_T"] <= rc_mean <= row["ci95_upper_T"]
            )
            rows.append(row)

    if args.output.exists():
        raise FileExistsError(args.output)
    temporary = args.output.with_name(args.output.name + f".tmp.{os.getpid()}")
    temporary.mkdir(parents=True)
    fields = list(rows[0])
    with (temporary / "cross_validated_ages.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    np.save(temporary / "chromosome_log_calibration_weights.npy", all_log_weights)
    result = {
        "method": "independent chromosome scale/epsilon calibration, then shared-age likelihood product",
        "likelihood_dir": str(source.resolve()),
        "radiocarbon": str(args.radiocarbon.resolve()),
        "chromosomes": chromosomes,
        "n_dated": len(dated),
        "folds": fold_map,
        "chromosome_calibration": chromosome_summaries,
        "cross_validation_metrics": base.accuracy_metrics(rows),
        "prior_accounting": "scale/epsilon prior applied once per chromosome; age prior applied once; held-out radiocarbon excluded from genetic age posterior",
    }
    (temporary / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    os.replace(temporary, args.output)
    print(json.dumps(result, indent=2), flush=True)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--likelihood-dir", type=Path, required=True)
    parser.add_argument("--radiocarbon", type=Path, required=True)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    run(parse_args())
