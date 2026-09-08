import json
import os
import sys
from pathlib import Path

import numpy as np
import tskit

sys.path.insert(0, "/Users/jeffreyross-ibarra/src/dnAging")
from posterior_sample_age_infer import load_table
from precompute_freq_trajectory_moments import MomentEngine

ROOT = Path(os.environ.get("SIM_ROOT", "/Users/jeffreyross-ibarra/Desktop/bigsims"))
NAME = os.environ.get("SIM_NAME", "simulation_01")
TABLE = Path(os.environ.get("SIM_TABLE", "/tmp/dnAging-archive2/frequency_table_ne50k_age100.npz"))
OUT = Path(os.environ.get("SIM_OUT", "/tmp/dnAging-bigsims-result.npz"))
NE = float(os.environ.get("SIM_NE", "50000"))
CUTOFF = float(os.environ.get("SIM_CUTOFF", str(6 * NE)))
EPSILONS = tuple(float(x) for x in os.environ.get("SIM_EPSILONS", "0.0,0.01").split(","))
N_BLOCKS = 100
ALIGN_BY_SITE_ID = bool(int(os.environ.get("SIM_ALIGN_BY_SITE_ID", "0")))
MIN_DERIVED_FREQ = float(os.environ.get("SIM_MIN_DERIVED_FREQ", "0"))
DROP_RARE_CARRIED = float(os.environ.get("SIM_DROP_RARE_CARRIED", "0"))


def vcf_alt_counts(path):
    out = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            calls = [x.split(":", 1)[0] for x in f[9:]]
            if all(x in ("0", "1") for x in calls):
                key = len(out) if ALIGN_BY_SITE_ID else int(f[1])
                out[key] = sum(x == "1" for x in calls)
    return out


def vcf_carried(path):
    out = set()
    with path.open() as fh:
        row = 0
        for line in fh:
            if not line.startswith("#"):
                f = line.rstrip().split("\t")
                if f[9].split(":", 1)[0] == "1":
                    out.add(row if ALIGN_BY_SITE_ID else int(f[1]))
                row += 1
    return out


def denominator_rows(age, n=26):
    eng = MomentEngine(n)
    out = np.full((n, len(age)), np.nan)
    x0 = 1 / (2 * NE)
    for ia, a in enumerate(age):
        moments = eng._moms(a / (2 * NE), x0)
        for d in range(1, n):
            value = sum(c * moments[k] for k, c in eng.coeff[d].items())
            if value > 0:
                out[d - 1, ia] = value
    return out


def interpolate_phi(plane, source_age, target_age, Tgrid):
    x = np.interp(np.log(np.maximum(target_age, source_age[0])),
                  np.log(source_age), np.arange(len(source_age)))
    i0 = np.floor(x).astype(int)
    i1 = np.minimum(i0 + 1, len(source_age) - 1)
    w = x - i0
    p = (1 - w[:, None]) * plane[i0] + w[:, None] * plane[i1]
    p[Tgrid[None, :] >= target_age[:, None]] = 0
    return p


def cumulative_trapezoid(y, x):
    dx = np.diff(x)
    pieces = 0.5 * dx.reshape((-1,) + (1,) * (y.ndim - 1)) * (y[:-1] + y[1:])
    return np.concatenate([np.zeros((1,) + y.shape[1:]), np.cumsum(pieces, axis=0)])


def eval_cumulative(cum, grid, points):
    points = np.asarray(points)
    i1 = np.clip(np.searchsorted(grid, points, side="right"), 1, len(grid) - 1)
    i0 = i1 - 1
    w = (points - grid[i0]) / (grid[i1] - grid[i0])
    if cum.ndim == 1:
        return (1 - w) * cum[i0] + w * cum[i1]
    return (1 - w[:, None]) * cum[i0] + w[:, None] * cum[i1]


def posterior_summary(ll, grid):
    density = np.exp(ll - ll.max())
    density /= np.trapezoid(density, grid)
    cdf = np.zeros(len(grid))
    cdf[1:] = np.cumsum(0.5 * np.diff(grid) * (density[:-1] + density[1:]))
    cdf /= cdf[-1]
    mean = np.trapezoid(grid * density, grid)
    var = np.trapezoid((grid - mean) ** 2 * density, grid)
    qs = np.interp([0.025, 0.5, 0.975], cdf, grid)
    return density, {
        "map": float(grid[ll.argmax()]), "mean": float(mean), "sd": float(np.sqrt(var)),
        "median": float(qs[1]), "ci95": [float(qs[0]), float(qs[2])],
    }


tab = load_table(TABLE)
grid = tab["Tgrid"]
panel = vcf_alt_counts(ROOT / f"{NAME}_modern.vcf")
carried = vcf_carried(ROOT / f"{NAME}_ancient.vcf")
ts = tskit.load(ROOT / f"{NAME}.trees")
records = {d: [] for d in range(1, 26)}
counts = {"tree_sites": ts.num_sites, "panel_vcf_sites": len(panel),
          "ancient_variant_records": len(carried), "used": 0,
          "multiple_mutations": 0, "missing_panel": 0, "monomorphic": 0,
          "root_mutation": 0, "age_filtered": 0}
for tree in ts.trees():
    for site in tree.sites():
        pos = site.id if ALIGN_BY_SITE_ID else int(site.position) + 1
        if pos not in panel:
            counts["missing_panel"] += 1
            continue
        if len(site.mutations) != 1:
            counts["multiple_mutations"] += 1
            continue
        mut = site.mutations[0]
        parent = tree.parent(mut.node)
        if parent == tskit.NULL:
            counts["root_mutation"] += 1
            continue
        d0 = panel[pos]
        if not (1 <= d0 < 26):
            counts["monomorphic"] += 1
            continue
        is_carried = pos in carried
        if d0 / 26 < MIN_DERIVED_FREQ:
            continue
        if is_carried and d0 / 26 < DROP_RARE_CARRIED:
            continue
        lo, hi = tree.time(mut.node), tree.time(parent)
        if lo >= CUTOFF:
            counts["age_filtered"] += 1
            continue
        block = min(int(site.position / (ts.sequence_length / N_BLOCKS)), N_BLOCKS - 1)
        records[d0].append((lo, min(hi, CUTOFF), is_carried, block))
        counts["used"] += 1

# Dense cumulative-integration grid. Every candidate T is an integration knot.
covered = tab["age"][tab["age"] <= CUTOFF]
parts = [np.geomspace(a, b, 9)[:-1] for a, b in zip(covered[:-1], covered[1:])]
fine_age = np.unique(np.concatenate([[0.0, tab["age"][0]], *parts,
                                     [covered[-1], CUTOFF],
                                     grid[(grid >= tab["age"][0]) & (grid <= CUTOFF)]]))
fine_age.sort()
den_source = denominator_rows(tab["age"])
block_ll = {eps: np.zeros((N_BLOCKS, len(grid))) for eps in EPSILONS}

for d0 in range(1, 26):
    source = den_source[d0 - 1]
    good = np.isfinite(source) & (source > 0)
    den = np.exp(np.interp(np.log(np.maximum(fine_age, tab["age"][0])),
                           np.log(tab["age"][good]), np.log(source[good])))
    den[0] = (1 / (2 * NE)) ** d0 * (1 - 1 / (2 * NE)) ** (26 - d0)
    plane = np.array(tab["table"][0, d0 - 1], copy=True)
    if np.isnan(plane).any():
        for ia in np.unique(np.argwhere(np.isnan(plane))[:, 0]):
            ok = np.isfinite(plane[ia])
            if ok.any():
                plane[ia, ~ok] = np.interp(grid[~ok], grid[ok], plane[ia, ok])
        # Moment cancellation can invalidate an entire very-young age row.
        # Fill any residual holes along the age axis before interpolation.
        log_age = np.log(tab["age"])
        for it in np.unique(np.argwhere(np.isnan(plane))[:, 1]):
            ok = np.isfinite(plane[:, it])
            if ok.any():
                plane[~ok, it] = np.interp(log_age[~ok], log_age[ok], plane[ok, it])
    phi = interpolate_phi(plane, tab["age"], fine_age, grid)
    cden = cumulative_trapezoid(den, fine_age)
    cnum = cumulative_trapezoid(den[:, None] * phi, fine_age)
    rec = records[d0]
    for start in range(0, len(rec), 1000):
        batch = rec[start:start + 1000]
        lo = np.array([x[0] for x in batch])
        hi = np.array([x[1] for x in batch])
        obs = np.array([x[2] for x in batch], dtype=bool)
        blocks = np.array([x[3] for x in batch])
        denom = eval_cumulative(cden, fine_age, hi) - eval_cumulative(cden, fine_age, lo)
        numer = eval_cumulative(cnum, fine_age, hi) - eval_cumulative(cnum, fine_age, lo)
        p = numer / denom[:, None]
        for eps in EPSILONS:
            q = np.clip(eps + (1 - 2 * eps) * p, 0.0, 1.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                site_ll = np.where(obs[:, None], np.log(q), np.log1p(-q))
            for block in np.unique(blocks):
                block_ll[eps][block] += site_ll[blocks == block].sum(axis=0)
    print(f"finished d0={d0}: {len(rec)} sites", flush=True)

result = {"counts": counts, "settings": {"Ne": NE, "cutoff": CUTOFF,
          "conditioning": "ratio of integrals", "blocks": N_BLOCKS}, "estimates": {}}
save = {"grid": grid}
rng = np.random.default_rng(20260905)
for eps in EPSILONS:
    total = block_ll[eps].sum(axis=0)
    density, summary = posterior_summary(total, grid)
    boot = np.empty(500)
    for b in range(len(boot)):
        sampled = rng.integers(0, N_BLOCKS, N_BLOCKS)
        boot[b] = grid[block_ll[eps][sampled].sum(axis=0).argmax()]
    summary["block_bootstrap_sd"] = float(boot.std(ddof=1))
    summary["block_bootstrap_95"] = [float(x) for x in np.percentile(boot, [2.5, 97.5])]
    result["estimates"][str(eps)] = summary
    save[f"ll_epsilon_{eps}"] = total
    save[f"density_epsilon_{eps}"] = density
    save[f"bootstrap_map_epsilon_{eps}"] = boot

np.savez_compressed(OUT, **save)
print(json.dumps(result, indent=2), flush=True)
print(f"saved {OUT}")
