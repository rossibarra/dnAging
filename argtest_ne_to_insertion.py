#!/usr/bin/env python3
"""Convert argtest's coalescence-Ne estimates into insertion's demography format.

`coalescence_ne_plots_from_ts.py` writes `*coalescence-ne-estimates.tsv` with the
columns `series / replicate_index / tree_file / bin_index / time_left /
time_right / ... / effective_population_size`, carrying one `replicate` series
per post-burnin ARG draw plus a `posterior_mean` series.  The column names
already match what `insertion_likelihood.PiecewiseConstantNe.from_tsv` reads, but
the file is not directly loadable, for three reasons:

1. **It does not start at zero.**  The quantile/log grid pads zero and infinity
   onto the time bins so tskit accepts it, then *drops* those padded intervals
   before writing.  `PiecewiseConstantNe.validate` requires `left[0] == 0`.
2. **It carries several series.**  `from_tsv(series=...)` can select one, but the
   caller has to know which; the posterior mean is the estimate of record.
3. **The final epoch is bounded.**  The ancestral epoch should extend to
   infinity, which `validate` permits only in the last position.

Bins where the estimate is not finite and positive (no coalescence mass in that
window) cannot be represented at all, so they are filled from the nearest valid
neighbour and the count is reported rather than silently absorbed.

The output is written with `series: estimated` so it cannot be confused with the
simulation's own `constant_ne_epochs.tsv`, whose series names the demography and
holds the *true* one.  Both are in the same schema, so the same loader
reads either and an estimated-versus-true comparison is a one-line swap.

Every file written here is loaded back through insertion's own parser before the
script exits, so a file that this produces is one that insertion accepts.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from insertion_likelihood import PiecewiseConstantNe


def read_series(path: Path, series: str):
    with path.open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t")
                if r.get("series") == series]
    if not rows:
        raise SystemExit(f"no rows with series={series!r} in {path}")
    rows.sort(key=lambda r: float(r["time_left"]))
    left = np.array([float(r["time_left"]) for r in rows])
    right = np.array([float(r["time_right"]) for r in rows])
    size = np.array([float(r["effective_population_size"]) for r in rows])
    return left, right, size


def to_insertion_epochs(left, right, size, floor, ceiling):
    """Make the grid contiguous, zero-anchored, and infinite at the top."""
    valid = np.isfinite(size) & (size > 0)
    if not valid.any():
        raise SystemExit("no finite positive Ne estimates to convert")
    filled = int((~valid).sum())
    if filled:
        # Nearest valid neighbour, preferring the younger side: an empty bin is
        # a statement about missing coalescence mass, not about Ne being small.
        idx = np.where(valid, np.arange(len(size)), -1)
        fwd = np.maximum.accumulate(idx)
        bwd = np.minimum.accumulate(np.where(valid, np.arange(len(size)),
                                             len(size))[::-1])[::-1]
        take = np.where(fwd >= 0, fwd, bwd)
        size = size[np.clip(take, 0, len(size) - 1)]

    left = left.copy()
    right = right.copy()
    left[0] = 0.0                       # anchor at the present
    left[1:] = right[:-1]               # enforce contiguity
    right[-1] = ceiling                 # ancestral epoch
    keep = right > left
    left, right, size = left[keep], right[keep], size[keep]
    size = np.clip(size, floor, None)
    return left, right, size, filled


def write_tsv(path: Path, left, right, size, series="estimated"):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["series", "time_left", "time_right",
                         "effective_population_size"])
        for a, b, n in zip(left, right, size):
            writer.writerow([series, f"{a:.10g}",
                             "inf" if not np.isfinite(b) else f"{b:.10g}",
                             f"{n:.10g}"])


def compare_to_truth(truth_path: Path, estimate: PiecewiseConstantNe):
    """Harmonic-mean Ne over the sampled-age window, estimated versus true.

    The harmonic mean is the quantity that matters: coalescence, and therefore
    the diffusion clock, integrates 1/(2Ne), so an arithmetic comparison would
    flatter a trajectory that gets the deep past wrong.
    """
    # The truth file's series names the demography ("variable" when drawn per
    # replicate, e.g. "three_epoch" when fixed), so detect it instead of assuming.
    with Path(truth_path).open() as handle:
        names = {r.get("series") for r in csv.DictReader(handle, delimiter="\t")}
    names.discard(None)
    truth = PiecewiseConstantNe.from_tsv(
        truth_path, series=(names.pop() if len(names) == 1 else "variable"))
    grid = np.logspace(0, np.log10(10_000), 400)
    t_est = estimate.at(grid)
    t_true = truth.at(grid)
    return {
        "harmonic_mean_ne_estimated_0_10k": float(1.0 / np.mean(1.0 / t_est)),
        "harmonic_mean_ne_true_0_10k": float(1.0 / np.mean(1.0 / t_true)),
        "median_log2_ratio_0_10k": float(np.median(np.log2(t_est / t_true))),
        "max_abs_log2_ratio_0_10k": float(np.max(np.abs(np.log2(t_est / t_true)))),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--estimates", type=Path, required=True,
                   help="argtest *coalescence-ne-estimates.tsv")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--series", default="posterior_mean")
    p.add_argument("--truth", type=Path, default=None,
                   help="the simulation's constant_ne_epochs.tsv, for QC")
    p.add_argument("--floor", type=float, default=1.0,
                   help="minimum representable Ne")
    p.add_argument("--ceiling", type=float, default=float("inf"),
                   help="right edge of the ancestral epoch (default infinite)")
    args = p.parse_args()

    left, right, size = read_series(args.estimates, args.series)
    left, right, size, filled = to_insertion_epochs(left, right, size,
                                                    args.floor, args.ceiling)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output, left, right, size)

    # Load it back through insertion's parser: the file is only useful if the
    # consumer accepts it, so assert that rather than assume it.
    loaded = PiecewiseConstantNe.from_tsv(args.output, series="estimated")
    summary = {
        "estimates": str(args.estimates),
        "output": str(args.output),
        "series_used": args.series,
        "epochs": int(len(loaded.size)),
        "bins_filled_from_neighbour": filled,
        "time_range": [float(loaded.left[0]), float(loaded.right[-1])],
        "ne_range": [float(loaded.size.min()), float(loaded.size.max())],
        "loads_with_insertion_parser": True,
    }
    if args.truth is not None and args.truth.exists():
        summary.update(compare_to_truth(args.truth, loaded))
    args.output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
