#!/usr/bin/env python3
"""T2/T3: how much does the ARG-draw mixture actually integrate?

Reads the per-site, per-draw ALT frequencies that `--save-epsilon-data` already
writes (`epsilon_calibration_data.npz`: phi_alt with axes site x ARG_draw x
sample-age grid, plus the ancient calls), and reconstructs the per-draw
log-likelihoods without rerunning inference.  From those it reports two things.

T2, the effective sample size of the draw mixture.  Eq. (11) marginalises the
ARG posterior as a mixture over G draws *outside* the site product, which is the
correct order -- but sum_i log ell_ig is O(n_sites), so the mixture weights
w_g = exp(sum_i log ell_ig - max) can collapse onto a single draw.  ESS =
(sum w)^2 / sum w^2 measures that: ESS ~ G means the tree posterior is being
integrated, ESS ~ 1 means the answer is conditional on the modal draw and the
correct ordering buys nothing beyond picking that draw.  ESS is a function of T,
so it is reported at the MAP.

T3, the size of the effect the ordering is there to capture.  The same data give
both orders directly:

  correct   log L_c(T) = logsumexp_g sum_i log ell_ig(T) - log G
  swapped   log L_c(T) = sum_i log ell_i(T; mean_g phi_ig(T))

The gap between the two posteriors is the empirical magnitude of the across-site
age dependence carried by the draw index -- the thing that cannot be recovered
by marginalising per site.  (For pseudo-haploid calls ell is linear in phi, so
the swapped form is exactly the old per-site draw average, not an approximation
of it.)

Draws are mixed per chromosome and chromosome log-likelihoods are then summed,
following eq. (12); the ESS is therefore a per-chromosome quantity.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

GROWTH_SITES = (1, 3, 10, 30, 100, 300, 1000, 3000)


def load_part(path: Path):
    d = np.load(path, allow_pickle=True)
    meta = json.loads(str(d["meta"]))
    return {"chrom": str(meta.get("chrom", path.parent.name)),
            "draw_ids": meta.get("store_draw_ids"),
            "ploidy": int(meta.get("ploidy", 1)),
            "Tgrid": d["Tgrid"],
            "phi": d["phi_alt"].astype(np.float64),      # (S, G, nT)
            "alt": d["observed_alt"].astype(np.float64),  # (S, n_samples)
            "called": d["called"].astype(np.float64),
            "samples": [str(s) for s in d["samples"]]}


def per_draw_loglik(part, eps, site_index=None):
    """(n_samples, G, nT) sum over sites of the per-draw log-likelihood."""
    phi = part["phi"] if site_index is None else part["phi"][site_index]
    alt = part["alt"] if site_index is None else part["alt"][site_index]
    called = part["called"] if site_index is None else part["called"][site_index]
    qA = np.clip(eps + (1 - 2 * eps) * phi, 1e-300, 1.0)
    S, G, nT = qA.shape
    logA = np.log(qA).reshape(S, G * nT)
    logR = np.log(np.clip(1.0 - qA, 1e-300, 1.0)).reshape(S, G * nT)
    # a_i in {0,1} and c_i - a_i in {0,1}: one Bernoulli per called site.
    a = alt * (called >= 1)
    r = np.clip(called - alt, 0.0, None) * (called >= 1)
    out = a.T @ logA + r.T @ logR
    return out.reshape(-1, G, nT)


def swapped_loglik(part, eps, site_index=None):
    """(n_samples, nT): phi averaged over draws per site, then the site product."""
    phi = part["phi"] if site_index is None else part["phi"][site_index]
    alt = part["alt"] if site_index is None else part["alt"][site_index]
    called = part["called"] if site_index is None else part["called"][site_index]
    qA = np.clip(eps + (1 - 2 * eps) * phi.mean(axis=1), 1e-300, 1.0)
    logA = np.log(qA)
    logR = np.log(np.clip(1.0 - qA, 1e-300, 1.0))
    a = alt * (called >= 1)
    r = np.clip(called - alt, 0.0, None) * (called >= 1)
    return a.T @ logA + r.T @ logR


def mix_draws(ll_by_draw):
    """logsumexp over the draw axis, minus log G: eq. (11)."""
    peak = ll_by_draw.max(axis=1)
    return peak + np.log(np.mean(np.exp(ll_by_draw - peak[:, None, :]), axis=1))


def ess_at(ll_by_draw, t_index):
    """Effective sample size of the draw weights, per sample, at one T."""
    x = ll_by_draw[:, :, t_index]
    w = np.exp(x - x.max(axis=1, keepdims=True))
    return (w.sum(axis=1) ** 2) / np.sum(w ** 2, axis=1)


def summarize_posterior(grid, logl):
    """MAP, posterior mean and central 90% interval under a uniform prior."""
    w = np.exp(logl - logl.max())
    area = np.trapezoid(w, grid)
    if not np.isfinite(area) or area <= 0:
        return dict(map_T=float("nan"), mean_T=float("nan"),
                    lo90=float("nan"), hi90=float("nan"))
    dens = w / area
    cdf = np.concatenate(([0.0], np.cumsum(np.diff(grid) *
                                           0.5 * (dens[1:] + dens[:-1]))))
    cdf /= cdf[-1]
    return dict(map_T=float(grid[int(np.argmax(logl))]),
                mean_T=float(np.trapezoid(dens * grid, grid)),
                lo90=float(np.interp(0.05, cdf, grid)),
                hi90=float(np.interp(0.95, cdf, grid)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("npz", type=Path, nargs="+",
                   help="epsilon_calibration_data.npz, one per chromosome part")
    p.add_argument("--epsilon", type=float, default=0.01)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--seed", type=int, default=20260910)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    parts = [load_part(path) for path in args.npz]
    for part in parts:
        if part["ploidy"] != 1:
            raise SystemExit(f"chrom {part['chrom']}: only ploidy 1 handled")
    grid = parts[0]["Tgrid"]
    samples = parts[0]["samples"]
    for part in parts[1:]:
        if not np.array_equal(part["Tgrid"], grid) or part["samples"] != samples:
            raise SystemExit("parts disagree on the T grid or the sample order")
    if len({tuple(part["draw_ids"]) for part in parts}) != 1:
        print("NOTE: draw ids differ between parts; draws are mixed per "
              "chromosome so this is admissible, but they are not matched.")

    n_draws = parts[0]["phi"].shape[1]
    total_correct = np.zeros((len(samples), len(grid)))
    total_swapped = np.zeros((len(samples), len(grid)))
    per_chrom = []
    growth = {n: [] for n in GROWTH_SITES}

    for part in parts:
        ll_by_draw = per_draw_loglik(part, args.epsilon)
        correct = mix_draws(ll_by_draw)
        swapped = swapped_loglik(part, args.epsilon)
        total_correct += correct
        total_swapped += swapped

        n_sites = part["phi"].shape[0]
        for s, name in enumerate(samples):
            t_map = int(np.argmax(correct[s]))
            per_chrom.append({
                "chrom": part["chrom"], "sample": name, "sites": n_sites,
                "map_T_correct": float(grid[t_map]),
                "map_T_swapped": float(grid[int(np.argmax(swapped[s]))]),
                "ess_at_map": float(ess_at(ll_by_draw, t_map)[s]),
                "ess_min_over_T": float(min(ess_at(ll_by_draw, t)[s]
                                            for t in range(0, len(grid), 8))),
                "spread_loglik_across_draws":
                    float(ll_by_draw[s, :, t_map].max()
                          - ll_by_draw[s, :, t_map].min()),
            })

        # ESS growth: how fast does the mixture collapse as sites accumulate?
        for n in GROWTH_SITES:
            if n > n_sites:
                continue
            idx = rng.choice(n_sites, size=n, replace=False)
            sub = per_draw_loglik(part, args.epsilon, site_index=idx)
            sub_map = np.argmax(mix_draws(sub), axis=1)
            growth[n].extend(float(ess_at(sub, int(t))[s])
                             for s, t in enumerate(sub_map))

    with (args.outdir / "per_chromosome_ess.tsv").open("w") as fh:
        cols = list(per_chrom[0])
        fh.write("\t".join(cols) + "\n")
        for row in per_chrom:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    with (args.outdir / "ess_growth.tsv").open("w") as fh:
        fh.write("n_sites\tn_observations\tess_mean\tess_median\tess_p05\tess_p95\n")
        for n in GROWTH_SITES:
            v = np.asarray(growth[n])
            if v.size == 0:
                continue
            fh.write(f"{n}\t{v.size}\t{v.mean():.4f}\t{np.median(v):.4f}"
                     f"\t{np.percentile(v, 5):.4f}\t{np.percentile(v, 95):.4f}\n")

    rows = []
    for s, name in enumerate(samples):
        a = summarize_posterior(grid, total_correct[s])
        b = summarize_posterior(grid, total_swapped[s])
        pa = np.exp(total_correct[s] - total_correct[s].max())
        pb = np.exp(total_swapped[s] - total_swapped[s].max())
        pa /= np.trapezoid(pa, grid); pb /= np.trapezoid(pb, grid)
        rows.append({
            "sample": name,
            "map_correct": a["map_T"], "map_swapped": b["map_T"],
            "map_shift": b["map_T"] - a["map_T"],
            "mean_correct": a["mean_T"], "mean_swapped": b["mean_T"],
            "width90_correct": a["hi90"] - a["lo90"],
            "width90_swapped": b["hi90"] - b["lo90"],
            "width90_ratio": ((b["hi90"] - b["lo90"]) / (a["hi90"] - a["lo90"])
                              if a["hi90"] > a["lo90"] else float("nan")),
            "total_variation": float(0.5 * np.trapezoid(np.abs(pa - pb), grid)),
        })
    with (args.outdir / "order_comparison.tsv").open("w") as fh:
        cols = list(rows[0])
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    ess = np.array([r["ess_at_map"] for r in per_chrom])
    tv = np.array([r["total_variation"] for r in rows])
    shift = np.array([r["map_shift"] for r in rows])
    ratio = np.array([r["width90_ratio"] for r in rows])
    print(f"draws per chromosome           : {n_draws}")
    print(f"chromosome parts               : {len(parts)}")
    print(f"sites total                    : {sum(p['phi'].shape[0] for p in parts)}")
    print(f"samples                        : {len(samples)}")
    print(f"epsilon                        : {args.epsilon}")
    print(f"ESS at MAP  (median / min / max): {np.median(ess):.3f} / "
          f"{ess.min():.3f} / {ess.max():.3f}   of {n_draws}")
    print(f"fraction of chrom x sample with ESS < 1.05 : "
          f"{float(np.mean(ess < 1.05)):.3f}")
    print(f"swapped-vs-correct total variation (median/max): "
          f"{np.median(tv):.4f} / {tv.max():.4f}")
    print(f"MAP shift, swapped - correct (median/max abs) : "
          f"{np.median(shift):.1f} / {np.abs(shift).max():.1f} generations")
    print(f"90% width ratio swapped/correct (median)      : "
          f"{np.median(ratio):.4f}")
    print(f"wrote {args.outdir}/per_chromosome_ess.tsv, ess_growth.tsv, "
          f"order_comparison.tsv")


if __name__ == "__main__":
    main()
