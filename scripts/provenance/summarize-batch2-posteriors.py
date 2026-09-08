import os
from pathlib import Path

import numpy as np

batch = Path(os.environ.get("TRUTH_ROOT", "/Users/jeffreyross-ibarra/Projects/mutrates/haploid_simulations_100mb_batch 2"))
truth = np.loadtxt(batch / "ancient_sample_ages.tsv", skiprows=1, dtype=int)
out = Path(os.environ.get("SUMMARY_OUT", "/Users/jeffreyross-ibarra/src/dnAging/results/batch2_100mb_true_vs_estimated.tsv"))
input_prefix = os.environ.get("POSTERIOR_PREFIX", "/tmp/dnAging-batch2-simulation_")
epsilon_key = os.environ.get("POSTERIOR_EPSILON", "0.01")


def quantile(grid, density, probs):
    cdf = np.zeros(len(grid))
    cdf[1:] = np.cumsum(np.diff(grid) * (density[:-1] + density[1:]) / 2)
    cdf /= cdf[-1]
    return np.interp(probs, cdf, grid)


def shortest_interval(grid, density, mass=0.95):
    starts = np.linspace(0, 1 - mass, 2001)
    lo = quantile(grid, density, starts)
    hi = quantile(grid, density, starts + mass)
    i = np.argmin(hi - lo)
    return float(lo[i]), float(hi[i])


rows = []
for sim, true_age in truth:
    d = np.load(f"{input_prefix}{sim:02d}.npz")
    grid = d["grid"]
    density = d[f"density_epsilon_{epsilon_key}"]
    density /= np.trapezoid(density, grid)
    mean = float(np.trapezoid(grid * density, grid))
    map_age = float(grid[np.argmax(density)])
    lo, hi = shortest_interval(grid, density)
    rows.append((sim, true_age, mean, map_age, lo, hi, lo <= true_age <= hi))

with out.open("w") as fh:
    fh.write("simulation\ttrue_age\tposterior_mean\tMAP\tHPD95_lower\tHPD95_upper\tcovered\n")
    for row in rows:
        fh.write("%d\t%d\t%.6f\t%.1f\t%.6f\t%.6f\t%d\n" % row)

true = np.array([r[1] for r in rows])
est = np.array([r[2] for r in rows])
print(out)
print(f"bias={np.mean(est-true):+.2f}")
print(f"RMSE={np.sqrt(np.mean((est-true)**2)):.2f}")
print(f"MAE={np.mean(np.abs(est-true)):.2f}")
print(f"coverage={sum(r[-1] for r in rows)}/{len(rows)}")
