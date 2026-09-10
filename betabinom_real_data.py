#!/usr/bin/env python3
"""ARG-conditioned beta-binomial age inference for the real maize data.

This wires the model on the ``betabinom`` branch to a full Singer tree sequence,
pre-filtered pseudo-haploid ancient calls, and a piecewise-variable Ne history.
One invocation evaluates one whole-genome ARG posterior draw.  ``--merge`` then
adds chromosomes within draws and marginalises the chromosome-coupled draws.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import tszip


CHROM_OFFSETS = np.array([
    0, 308_452_472, 552_127_664, 790_145_432, 1_040_475_893,
    1_266_829_343, 1_448_186_578, 1_633_995_495, 1_816_406_698,
    1_979_411_443,
], dtype=np.int64)
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


class VariableNePhi:
    """Interpolate beta-binomial numerator/denominator in diffusion duration."""

    def __init__(self, path):
        self.path = Path(path)
        self.u = np.load(self.path / "u.npy", mmap_mode="r")
        self.logu = np.log(self.u)
        self.left = np.load(self.path / "window_left.npy", mmap_mode="r")
        self.right = np.load(self.path / "window_right.npy", mmap_mode="r")
        self.ne = np.load(self.path / "window_ne.npy", mmap_mode="r")
        self.log_num = np.load(self.path / "log_num.npy", mmap_mode="r")
        self.log_den = np.load(self.path / "log_den.npy", mmap_mode="r")
        seg = (self.right - self.left) / (2.0 * self.ne)
        self.cum_left = np.concatenate([[0.0], np.cumsum(seg)])[:-1]

    def window(self, t):
        return np.clip(np.searchsorted(self.right, t, side="right"),
                       0, len(self.ne) - 1)

    def tau(self, t):
        values = np.asarray(t, dtype=float)
        idx = self.window(values)
        return self.cum_left[idx] + np.maximum(values - self.left[idx], 0.0) / (
            2.0 * self.ne[idx])

    def time_at_tau(self, target):
        cum_right = self.cum_left + (self.right - self.left) / (2.0 * self.ne)
        idx = int(np.clip(np.searchsorted(cum_right, target, side="right"),
                          0, len(self.ne) - 1))
        return float(self.left[idx] + (target - self.cum_left[idx]) * 2.0 * self.ne[idx])

    def _lookup(self, table, k, nt, sample_t, origin_t):
        origin = np.asarray(origin_t, dtype=float)
        sample = np.broadcast_to(np.asarray(sample_t, dtype=float), origin.shape)
        idx = self.window(origin)
        du = np.maximum(self.tau(origin) - self.tau(sample), self.u[0])
        lu = np.log(np.clip(du, self.u[0], self.u[-1]))
        out = np.empty(origin.shape, dtype=float)
        flat_out = out.ravel()
        flat_idx = idx.ravel()
        flat_lu = lu.ravel()
        for iw in np.unique(flat_idx):
            take = np.flatnonzero(flat_idx == iw)
            row = np.asarray(table[int(iw), :, k, nt], dtype=float)
            if not np.all(np.isfinite(row)):
                raise ValueError(f"missing beta-binomial table row k={k}, nT={nt}, window={iw}")
            flat_out[take] = np.interp(flat_lu[take], self.logu, row)
        return out

    @staticmethod
    def _logsumexp(a, axis=-1):
        peak = np.max(a, axis=axis, keepdims=True)
        return np.squeeze(peak, axis=axis) + np.log(
            np.sum(np.exp(a - peak), axis=axis))

    def integrate(self, k, nt, sample_t, origin_t, weights):
        ln = self._lookup(self.log_num, k, nt, sample_t, origin_t)
        ld = self._lookup(self.log_den, k, nt, sample_t, origin_t)
        lw = np.log(np.broadcast_to(np.asarray(weights, float), ln.shape))
        value = np.exp(self._logsumexp(lw + ln) - self._logsumexp(lw + ld))
        return np.clip(value, 0.0, 1.0)


def orientation(site, mutation, ref, alt):
    """Return whether the ARG-derived state is ancient ALT, including strand flip."""
    ancestral = site.ancestral_state.upper()
    derived = mutation.derived_state.upper()
    if len(ancestral) != 1 or len(derived) != 1:
        return None, "non_snv_tree_state"
    if {ancestral, derived} == {ref, alt}:
        return derived == alt, "exact"
    cref, calt = COMPLEMENT[ref], COMPLEMENT[alt]
    if {ancestral, derived} == {cref, calt}:
        return derived == calt, "strand_complement"
    return None, "allele_mismatch"


def run_draw(args):
    phi = VariableNePhi(args.moment_table)
    grid_data = np.load(args.grid_table, allow_pickle=True)
    grid = np.asarray(grid_data["Tgrid"], dtype=float)
    cutoff = phi.time_at_tau(args.max_mutation_tau)
    weights = np.ones(args.n_quad, dtype=float)
    weights[[0, -1]] = 0.5

    calls_by_chrom = {}
    samples = None
    for chrom in range(1, 11):
        data = np.load(args.calls / f"chr{chrom}.npz", allow_pickle=False)
        current = list(data["samples"].astype(str))
        if samples is None:
            samples = current
        elif samples != current:
            raise SystemExit(f"sample order differs in chr{chrom} calls")
        calls_by_chrom[chrom] = data

    print(f"decompressing {args.tree}", flush=True)
    ts = tszip.decompress(args.tree)
    if ts.num_samples != 26:
        raise SystemExit(f"tree sequence has {ts.num_samples} samples, expected 26")
    site_positions = np.asarray(ts.tables.sites.position)
    ll_by_chrom = np.zeros((10, len(samples), len(grid)), dtype=np.float64)
    used = {}
    counts = {}
    tree = ts.first()

    for chrom in range(1, 11):
        calls = calls_by_chrom[chrom]
        local_pos = np.asarray(calls["position"], dtype=np.int64)
        global_pos = local_pos + CHROM_OFFSETS[chrom - 1]
        refs, alts = calls["ref"].astype(str), calls["alt"].astype(str)
        observed_alt = np.asarray(calls["observed_alt"], dtype=np.uint8)
        called = np.asarray(calls["called"], dtype=np.uint8)
        chrom_ll = ll_by_chrom[chrom - 1]
        kept = []
        stat = {
            "targets": len(local_pos), "used": 0, "absent_from_tree": 0,
            "multiple_mutations": 0, "root_mutation": 0, "allele_mismatch": 0,
            "strand_complement": 0, "age_filtered": 0, "invalid_topology": 0,
        }
        for row, (lp, gp) in enumerate(zip(local_pos, global_pos)):
            si = int(np.searchsorted(site_positions, float(gp)))
            if si >= len(site_positions) or site_positions[si] != float(gp):
                stat["absent_from_tree"] += 1
                continue
            site = ts.site(si)
            if len(site.mutations) != 1:
                stat["multiple_mutations"] += 1
                continue
            mutation = site.mutations[0]
            is_alt, match = orientation(site, mutation, refs[row], alts[row])
            if is_alt is None:
                stat["allele_mismatch"] += 1
                continue
            if match == "strand_complement":
                stat["strand_complement"] += 1
            tree.seek(float(gp))
            parent = tree.parent(mutation.node)
            if parent == -1:
                stat["root_mutation"] += 1
                continue
            tc, tp = tree.time(mutation.node), tree.time(parent)
            if tc >= cutoff:
                stat["age_filtered"] += 1
                continue
            tp = min(tp, cutoff)
            if tp <= tc:
                stat["age_filtered"] += 1
                continue

            internal = np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
            ntv = ts.num_samples - np.searchsorted(internal, grid, side="right")
            descendants = np.sort([
                tree.time(u) for u in tree.nodes(root=mutation.node)
                if tree.is_internal(u)
            ])
            nl = tree.num_samples(mutation.node)
            kv = nl - np.searchsorted(descendants, grid, side="right")
            probability = np.zeros(len(grid), dtype=float)

            above = np.flatnonzero(
                (grid < tc) & (ntv >= 2) & (kv >= 1) & (kv < ntv))
            if above.size:
                origin = np.linspace(tc, tp, args.n_quad)
                key = kv[above] * 100 + ntv[above]
                for value in np.unique(key):
                    take = above[key == value]
                    k, nt = int(value // 100), int(value % 100)
                    origins = np.broadcast_to(origin, (len(take), args.n_quad))
                    times = np.broadcast_to(grid[take, None], origins.shape)
                    probability[take] = phi.integrate(k, nt, times, origins, weights)

            crossing = np.flatnonzero((grid >= tc) & (grid < tp) & (ntv >= 2))
            for it in crossing:
                origin = np.linspace(grid[it], tp, args.n_quad)
                value = float(phi.integrate(
                    1, int(ntv[it]), np.full(args.n_quad, grid[it]), origin, weights))
                probability[it] = value * (tp - grid[it]) / (tp - tc)

            if not np.all(np.isfinite(probability)):
                stat["invalid_topology"] += 1
                continue
            r = np.clip(args.epsilon + (1.0 - 2.0 * args.epsilon) * probability,
                        1e-300, 1.0 - 1e-15)
            derived = observed_alt[row] if is_alt else called[row] - observed_alt[row]
            chrom_ll += (derived[:, None] * np.log(r)[None, :]
                         + (called[row] - derived)[:, None] * np.log1p(-r)[None, :])
            kept.append(int(lp))
            stat["used"] += 1

        used[chrom] = np.asarray(kept, dtype=np.int64)
        counts[str(chrom)] = stat
        print(f"chr{chrom}: used {stat['used']}/{stat['targets']}", flush=True)

    args.output.mkdir(parents=True, exist_ok=False)
    np.save(args.output / "grid.npy", grid)
    np.save(args.output / "ll_by_chrom.npy", ll_by_chrom)
    (args.output / "samples.txt").write_text("\n".join(samples) + "\n")
    for chrom, positions in used.items():
        np.save(args.output / f"used_chr{chrom}.npy", positions)
    (args.output / "run.json").write_text(json.dumps({
        "model": "ARG-conditioned beta-binomial",
        "tree": str(args.tree.resolve()),
        "draw_id": args.draw_id,
        "epsilon": args.epsilon,
        "n_quad": args.n_quad,
        "max_mutation_tau": args.max_mutation_tau,
        "max_mutation_age_generations": cutoff,
        "moment_table": str(args.moment_table.resolve()),
        "counts_by_chromosome": counts,
    }, indent=2) + "\n")


def trapz(y, x):
    return np.trapezoid(y, x)


def summarize(grid, logp):
    weight = np.exp(logp - np.max(logp))
    density = weight / trapz(weight, grid)
    cdf = np.concatenate([[0.0], np.cumsum(
        (density[:-1] + density[1:]) * 0.5 * np.diff(grid))])
    cdf /= cdf[-1]
    quantile = lambda p: float(np.interp(p, cdf, grid))
    return (float(grid[np.argmax(logp)]), float(trapz(density * grid, grid)),
            quantile(0.5), quantile(0.025), quantile(0.975))


def merge_draws(args):
    parts = [Path(path) for path in args.merge]
    samples = (parts[0] / "samples.txt").read_text().split()
    grid = np.load(parts[0] / "grid.npy")
    draw_ll = []
    reference_used = None
    run_meta = []
    for part in parts:
        if (part / "samples.txt").read_text().split() != samples:
            raise SystemExit(f"sample order differs in {part}")
        if not np.array_equal(np.load(part / "grid.npy"), grid):
            raise SystemExit(f"grid differs in {part}")
        used = [np.load(part / f"used_chr{chrom}.npy") for chrom in range(1, 11)]
        if reference_used is None:
            reference_used = used
        elif any(not np.array_equal(a, b) for a, b in zip(reference_used, used)):
            raise SystemExit(
                f"retained sites differ in {part}; build a common-draw mask before merging")
        ll = np.load(part / "ll_by_chrom.npy")
        if ll.shape != (10, len(samples), len(grid)) or not np.all(np.isfinite(ll)):
            raise SystemExit(f"invalid likelihood array in {part}: {ll.shape}")
        draw_ll.append(ll.sum(axis=0))
        run_meta.append(json.loads((part / "run.json").read_text()))
    draw_ll = np.stack(draw_ll)
    peak = np.max(draw_ll, axis=0)
    marginal = peak + np.log(np.mean(np.exp(draw_ll - peak[None, :, :]), axis=0))

    args.output.mkdir(parents=True, exist_ok=False)
    np.save(args.output / "grid.npy", grid)
    np.save(args.output / "ll_marginal.npy", marginal)
    (args.output / "samples.txt").write_text("\n".join(samples) + "\n")
    with (args.output / "ages_table.tsv").open("w") as handle:
        handle.write("sample\tmap_T\tmean_T\tmedian_T\tci95_lower_T\tci95_upper_T\n")
        for sample, lp in zip(samples, marginal):
            values = summarize(grid, lp)
            handle.write(sample + "\t" + "\t".join(f"{v:.6g}" for v in values) + "\n")
    (args.output / "run.json").write_text(json.dumps({
        "model": "ARG-conditioned beta-binomial",
        "n_arg_draws": len(parts),
        "draw_ids": [meta["draw_id"] for meta in run_meta],
        "epsilon": run_meta[0]["epsilon"],
        "sites_by_chromosome": {
            str(chrom): len(reference_used[chrom - 1]) for chrom in range(1, 11)
        },
        "marginalisation": "sum chromosomes within each ARG draw, then uniform draw mixture",
    }, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=Path)
    parser.add_argument("--draw-id", type=int)
    parser.add_argument("--calls", type=Path)
    parser.add_argument("--moment-table", type=Path)
    parser.add_argument("--grid-table", type=Path)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--n-quad", type=int, default=12)
    parser.add_argument("--max-mutation-tau", type=float, default=3.0)
    parser.add_argument("--merge", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.merge:
        merge_draws(args)
    else:
        missing = [name for name in ("tree", "draw_id", "calls", "moment_table", "grid_table")
                   if getattr(args, name) is None]
        if missing:
            parser.error("draw run missing: " + ", ".join("--" + x.replace("_", "-") for x in missing))
        run_draw(args)


if __name__ == "__main__":
    main()
