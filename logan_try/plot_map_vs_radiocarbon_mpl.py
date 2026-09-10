#!/usr/bin/env python3
"""Render the chromosome-1 MAP versus radiocarbon comparison with Matplotlib."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RADIOCARBON = ROOT / "radiocarbon_selection.tsv"
ESTIMATES = ROOT / "results" / "smoke_restart" / "1" / "ages_table.tsv"
OUTPUT = ESTIMATES.parent / "map_vs_radiocarbon.png"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


rc = {row["sample_id"]: row for row in read_tsv(RADIOCARBON)}
estimates = read_tsv(ESTIMATES)
missing = [row["sample"] for row in estimates if row["sample"] not in rc]
if missing:
    raise SystemExit(f"Missing radiocarbon ages for: {', '.join(missing)}")

x = np.array([float(rc[row["sample"]]["radiocarbon_age_bp"]) for row in estimates])
xerr = np.array(
    [
        float(rc[row["sample"]]["radiocarbon_error_1sigma"])
        if rc[row["sample"]]["radiocarbon_error_1sigma"]
        else np.nan
        for row in estimates
    ]
)
y = np.array([float(row["map_T"]) for row in estimates])
ylo = np.array([float(row["ci95_lower_T"]) for row in estimates])
yhi = np.array([float(row["ci95_upper_T"]) for row in estimates])

r = float(np.corrcoef(x, y)[0, 1])
rmse = float(np.sqrt(np.mean((y - x) ** 2)))
at_zero = y == 0

fig, ax = plt.subplots(figsize=(8.4, 7.4), constrained_layout=True)
ax.vlines(x, ylo, yhi, color="#087e8b", alpha=0.18, linewidth=0.8, zorder=1)
has_err = np.isfinite(xerr)
ax.errorbar(
    x[has_err], y[has_err], xerr=xerr[has_err], fmt="none",
    ecolor="#087e8b", alpha=0.35, linewidth=0.8, capsize=1.5, zorder=2,
)
ax.scatter(x[~at_zero], y[~at_zero], s=38, color="#087e8b", edgecolor="white",
           linewidth=0.7, alpha=0.88, zorder=3, label="MAP estimate")
ax.scatter(x[at_zero], y[at_zero], s=48, color="#d95f02", edgecolor="white",
           linewidth=0.7, zorder=4, label=f"MAP at T=0 (n={at_zero.sum()})")
ax.plot([0, 6000], [0, 6000], linestyle="--", color="#6b7280", linewidth=1.2,
        label="1:1")

ax.set(xlim=(-100, 6100), ylim=(-100, 6100),
       xlabel="Conventional radiocarbon age (¹⁴C years BP)",
       ylabel="Inferred MAP age (generations; 1 year/generation)")
ax.set_title("Chromosome 1 MAP age vs. radiocarbon age", weight="bold", pad=13)
ax.text(0.5, 1.01, "50 ancient maize samples; one ARG draw", transform=ax.transAxes,
        ha="center", va="bottom", fontsize=9.5, color="#59636e")
ax.text(0.025, 0.965, f"Pearson r = {r:.2f}\nRMSE = {rmse:,.0f} years",
        transform=ax.transAxes, va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#c7cdd4"))
ax.grid(True, color="#dfe3e8", linewidth=0.7)
ax.set_axisbelow(True)
ax.legend(loc="lower right", frameon=True)
fig.text(
    0.5, 0.005,
    "Vertical lines: inferred 95% intervals; horizontal lines: radiocarbon ±1σ where reported. "
    "Radiocarbon ages are not calendar-calibrated.",
    ha="center", fontsize=8, color="#68727d",
)
fig.savefig(OUTPUT, dpi=220, facecolor="white")
print(f"wrote {OUTPUT}")
print(f"Pearson r={r:.4f}; RMSE={rmse:.1f}; n={len(x)}")
