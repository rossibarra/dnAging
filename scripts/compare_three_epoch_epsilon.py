"""Compare estimated-ARG epsilon sensitivity with truth and true-ARG control."""
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("msprime_three_epoch/insertion_three_chrom")


def read(path):
    with path.open() as handle:
        return {row["sample"]: row for row in csv.DictReader(handle, delimiter="\t")}


def main():
    ages = json.loads(Path("msprime_three_epoch/simulations/replicate_001/metadata.json").read_text())["true_ancient_ages"]
    names = [f"ancient_{i:02d}" for i in range(1, 11)]
    truth = np.asarray(ages)
    out = ROOT / "epsilon_comparison_with_10pct"
    out.mkdir(exist_ok=False)
    stats = []

    def load(tag, merged):
        base = ROOT / tag / ("merged" if merged else "")
        rows = read(base / "ages_table.tsv")
        boots = read(base / "block_bootstrap_ci.tsv")
        maps = np.array([float(rows[n]["map_T"]) for n in names])
        lower = np.array([float(boots[n]["bootstrap_map_ci95_lower"]) for n in names])
        upper = np.array([float(boots[n]["bootstrap_map_ci95_upper"]) for n in names])
        error = maps - truth
        stats.append({"condition": tag, "bias": float(error.mean()), "mae": float(abs(error).mean()),
                      "rmse": float(np.sqrt(np.mean(error**2))),
                      "bootstrap_coverage": float(np.mean((lower <= truth) & (truth <= upper)))})
        return maps, lower, upper

    control = load("true_arg_true_ne", False)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True, sharey=True)
    for ax, ne in zip(axes, ("true_ne", "estimated_ne")):
        ax.plot([0, 10000], [0, 10000], "k--", lw=1)
        ax.scatter(truth, control[0], color="0.5", marker="x", label="True ARG + true Ne, ε=0")
        for eps, color in zip(("0", "0.001", "0.01", "0.10"), ("#d95f02", "#1b9e77", "#7570b3", "#e7298a")):
            tag = f"estimated_arg_{ne}" + (f"_eps_{eps}_4g" if eps != "0" else "")
            maps, lower, upper = load(tag, True)
            # Draw interval endpoints directly; percentile intervals need not contain the original MAP.
            ax.vlines(truth, lower, upper, color=color, alpha=.65, lw=1.5)
            ax.scatter(truth, maps, color=color, s=25, label=f"Estimated ARG, ε={eps}", zorder=3)
        ax.set(title="True Ne" if ne == "true_ne" else "Estimated Ne", xlabel="True age (generations)",
               xlim=(0, 10000), ylim=(0, 10000))
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("MAP age (generations); 95% block-bootstrap intervals")
    fig.tight_layout()
    fig.savefig(out / "map_vs_true.png", dpi=200)
    (out / "metrics.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
