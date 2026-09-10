#!/usr/bin/env python3
"""Profile a global ancient-call error rate from saved dated-sample data."""

import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np


parts_root = Path(sys.argv[1])
metadata_path = Path(sys.argv[2])
output_dir = Path(sys.argv[3])
output_dir.mkdir(parents=True, exist_ok=True)

with metadata_path.open(newline="") as handle:
    metadata = {r["sample_id"]: r for r in csv.DictReader(handle, delimiter="\t")}

part_paths = [parts_root / str(chrom) / "epsilon_calibration_data.npz"
              for chrom in range(1, 11)]
parts = [np.load(path, allow_pickle=False) for path in part_paths]
samples = parts[0]["samples"].tolist()
if any(part["samples"].tolist() != samples for part in parts[1:]):
    raise SystemExit("sample order differs among chromosome calibration archives")
if any(part["phi_alt"].shape[1] < 1 for part in parts):
    raise SystemExit("a chromosome calibration archive contains no ARG draws")
n_arg_draws = int(parts[0]["phi_alt"].shape[1])
if any(part["phi_alt"].shape[1] != n_arg_draws for part in parts[1:]):
    raise SystemExit("ARG draw count differs among chromosome calibration archives")

raw_age = np.array([float(metadata[s]["radiocarbon_age_bp"]) for s in samples])
cal_mid = np.full(len(samples), np.nan)
for i, sample in enumerate(samples):
    value = metadata[sample]["cal_bp_95"]
    numbers = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", value)]
    if len(numbers) == 2:
        cal_mid[i] = sum(numbers) / 2

eps_grid = np.linspace(0.0, 0.49, 491)


def logmeanexp(values, axis):
    peak = np.max(values, axis=axis, keepdims=True)
    return (np.squeeze(peak, axis=axis)
            + np.log(np.mean(np.exp(values - peak), axis=axis)))


def interpolate_site_phi(part, age):
    grid = part["Tgrid"]
    age = float(np.clip(age, grid[0], grid[-1]))
    hi = int(np.searchsorted(grid, age, side="right"))
    if hi == 0:
        return part["phi_alt"][:, :, 0].astype(np.float64)
    if hi >= len(grid):
        return part["phi_alt"][:, :, -1].astype(np.float64)
    lo = hi - 1
    w = (age - grid[lo]) / (grid[hi] - grid[lo])
    return ((1 - w) * part["phi_alt"][:, :, lo]
            + w * part["phi_alt"][:, :, hi]).astype(np.float64)


def sample_profile(sample_index, age):
    profile = np.zeros(len(eps_grid))
    for part in parts:
        p = interpolate_site_phi(part, age)       # site x ARG draw
        alt = part["observed_alt"][:, sample_index].astype(bool)
        called = part["called"][:, sample_index].astype(bool)
        ref = called & ~alt
        for start in range(0, len(eps_grid), 32):
            stop = min(start + 32, len(eps_grid))
            eps = eps_grid[start:stop, None, None]
            q = np.clip(eps + (1 - 2 * eps) * p[None, :, :], 1e-300, 1 - 1e-15)
            by_draw = (np.log(q[:, alt, :]).sum(axis=1)
                       + np.log1p(-q[:, ref, :]).sum(axis=1))
            profile[start:stop] += logmeanexp(by_draw, axis=1)
    return profile


raw_profiles = np.vstack([sample_profile(i, age) for i, age in enumerate(raw_age)])
cal_profiles = np.full_like(raw_profiles, np.nan)
for i, age in enumerate(cal_mid):
    if np.isfinite(age):
        cal_profiles[i] = sample_profile(i, age)


def result_for(profiles, mask, label):
    used = profiles[mask]
    total = used.sum(axis=0)
    best = int(np.argmax(total))
    rng = np.random.default_rng(20260908)
    boot = np.empty(2000)
    for b in range(len(boot)):
        sampled = rng.integers(0, len(used), len(used))
        boot[b] = eps_grid[int(np.argmax(used[sampled].sum(axis=0)))]
    return {
        "label": label,
        "n_samples": int(mask.sum()),
        "epsilon_mle": float(eps_grid[best]),
        "sample_bootstrap_95": [float(x) for x in np.quantile(boot, [0.025, 0.975])],
        "at_search_boundary": bool(best in (0, len(eps_grid) - 1)),
        "total_log_likelihood": total,
    }


cal_mask = np.isfinite(cal_mid)
cal = result_for(cal_profiles, cal_mask, "calibrated-95%-range midpoint")
raw = result_for(raw_profiles, np.ones(len(samples), dtype=bool),
                 "raw conventional radiocarbon point age")

with (output_dir / "epsilon_profile.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["epsilon", "loglik_calibrated_midpoint_36", "loglik_raw_14C_50"])
    for e, a, b in zip(eps_grid, cal["total_log_likelihood"],
                       raw["total_log_likelihood"]):
        writer.writerow([f"{e:.3f}", f"{a:.12g}", f"{b:.12g}"])

summary = {
    "format_version": 1,
    "epsilon_search_grid": [float(eps_grid[0]), float(eps_grid[-1]), 0.001],
    "primary": {k: v for k, v in cal.items() if k != "total_log_likelihood"},
    "sensitivity": {k: v for k, v in raw.items() if k != "total_log_likelihood"},
    "bootstrap_unit": "sample",
    "n_bootstrap": 2000,
    "caveats": [
        "Primary calibrated ages use the midpoint of each reported 95% range, not the full calibration density.",
        "The primary subset omits samples lacking a calibrated range and does not span the oldest raw radiocarbon ages.",
        f"This composite-likelihood estimate marginalizes {n_arg_draws} ARG draw(s).",
        "Use the next independently dated sample set for validation, not further tuning.",
    ],
}
(output_dir / "epsilon_estimate.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))

for part in parts:
    part.close()
