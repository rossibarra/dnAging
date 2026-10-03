"""Overlay true and inferred Ne over the most recent 20,000 generations."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path("msprime_three_epoch")
fig, ax = plt.subplots(figsize=(8, 4.8))
for path, label, color in [
    (root / "demography/three_epoch.tsv", "True Ne", "#222222"),
    (root / "argtest/ne/estimated_ne.tsv", "Estimated Ne", "#377eb8"),
]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    values, edges = [], []
    for row in rows:
        left, right = float(row["time_left"]), float(row["time_right"])
        if left >= 20000:
            break
        edges.append(left)
        values.append(float(row["effective_population_size"]))
        end = min(right, 20000)
    edges.append(end)
    assert np.all(np.diff(edges) > 0)
    ax.stairs(values, edges, baseline=None, label=label, color=color, linewidth=2)
ax.set(xlim=(0, 20000), yscale="log", xlabel="Generations before present",
       ylabel="Diploid effective population size (Ne)", title="True and estimated Ne: most recent 20,000 generations")
ax.set_xticks([0, 1000, 5000, 10000, 15000, 20000])
ax.tick_params(axis="x", rotation=30)
ax.grid(axis="y", which="both", alpha=.2)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(root / "demography/true_vs_estimated_ne_20k_log.png", dpi=200)
fig.savefig(root / "demography/true_vs_estimated_ne_20k_log.pdf")
