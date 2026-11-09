#!/usr/bin/env python3
"""Profile joint time-scale/epsilon calibration by chromosome subset."""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))

import joint_scale_epsilon_calibration as joint
import radiocarbon_calibration as base


def subset_profiles(ll, age_grid, samples, rc_profiles, chromosome_indices):
    """Integrate radiocarbon densities after summing selected chromosomes."""
    selected = np.sum(ll[:, :, chromosome_indices, :, :], axis=2)
    sample_index = {sample: index for index, sample in enumerate(samples)}
    profiles = []
    for sample, (rc_age, rc_log_weight) in rc_profiles.items():
        if sample not in sample_index:
            raise ValueError(f"radiocarbon sample absent from likelihood: {sample}")
        profiles.append(joint.joint_radiocarbon_log_profile(
            age_grid, selected[:, :, sample_index[sample], :], rc_age, rc_log_weight
        ))
    return np.asarray(profiles)


def profile_row(label, chromosomes, profiles, scale_grid, epsilon_grid):
    profile = np.sum(profiles, axis=0)
    summary = joint.joint_summary(profile, scale_grid, epsilon_grid)
    raw_index = np.unravel_index(np.argmax(profile), profile.shape)
    scale_log = base.logsumexp(
        joint.joint_log_posterior(profile, scale_grid, epsilon_grid), axis=1
    )
    return {
        "subset": label,
        "chromosomes": ",".join(chromosomes),
        "n_chromosomes": len(chromosomes),
        "n_dated": len(profiles),
        "joint_scale_map": summary["joint_map"]["time_scale"],
        "joint_epsilon_map": summary["joint_map"]["epsilon"],
        "raw_scale_mle": float(scale_grid[raw_index[0]]),
        "raw_epsilon_mle": float(epsilon_grid[raw_index[1]]),
        "scale_posterior_mean": summary["time_scale"]["mean"],
        "scale_ci95_lower": summary["time_scale"]["ci95"][0],
        "scale_ci95_upper": summary["time_scale"]["ci95"][1],
        "scale_boundary_mass": summary["time_scale"]["boundary_mass"],
        "posterior_mass_scale_0p7": float(
            np.sum(np.exp(scale_log)[np.isclose(scale_grid, 0.7)])
        ),
        "posterior_mass_scale_0p9": float(
            np.sum(np.exp(scale_log)[np.isclose(scale_grid, 0.9)])
        ),
    }


def run(args):
    metadata = json.loads((args.draw_dir / "run.json").read_text())
    chromosomes = list(metadata["chromosomes"])
    age_grid = np.load(args.draw_dir / "grid.npy")
    scale_grid = np.load(args.draw_dir / "time_scale_grid.npy")
    epsilon_grid = np.load(args.draw_dir / "epsilon_grid.npy")
    ll = np.load(args.draw_dir / "ll_by_chrom_scale_epsilon.npy", mmap_mode="r")
    samples = (args.draw_dir / "samples.txt").read_text().split()
    expected = (
        len(scale_grid), len(epsilon_grid), len(chromosomes), len(samples), len(age_grid)
    )
    if ll.shape != expected:
        raise ValueError(f"likelihood shape {ll.shape}; expected {expected}")
    rc_profiles = base.read_radiocarbon(args.radiocarbon)

    subsets = []
    for size in range(1, len(chromosomes) + 1):
        for indices in itertools.combinations(range(len(chromosomes)), size):
            chosen = [chromosomes[index] for index in indices]
            label = "+".join(chosen)
            profiles = subset_profiles(ll, age_grid, samples, rc_profiles, list(indices))
            subsets.append(profile_row(
                label, chosen, profiles, scale_grid, epsilon_grid
            ))

    if args.output.exists():
        raise FileExistsError(args.output)
    temporary = args.output.with_name(args.output.name + f".tmp.{os.getpid()}")
    temporary.mkdir(parents=True)
    with (temporary / "chromosome_scale_profiles.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(subsets[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(subsets)
    summary = {
        "draw_dir": str(args.draw_dir.resolve()),
        "radiocarbon": str(args.radiocarbon.resolve()),
        "chromosomes": chromosomes,
        "n_subsets": len(subsets),
        "profiles": subsets,
    }
    (temporary / "run.json").write_text(json.dumps(summary, indent=2) + "\n")
    os.replace(temporary, args.output)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draw-dir", type=Path, required=True)
    parser.add_argument("--radiocarbon", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    run(parse_args())
