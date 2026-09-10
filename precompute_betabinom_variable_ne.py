#!/usr/bin/env python3
"""Precompute the beta-binomial ARG likelihood under piecewise-variable Ne.

The beta-binomial model needs the neutral frequency moments between a mutation's
absolute origin time ``t_i`` and a candidate sample time ``T``.  With variable
population size these depend on two quantities:

    u = tau(t_i) - tau(T),  tau(t) = integral_0^t ds / (2 Ne(s))
    x0 = 1 / (2 Ne(t_i))

``Ne(t_i)`` is piecewise constant in the supplied coalescence-Ne table.  We
therefore tabulate log numerator and log denominator by Ne window and diffusion
duration u.  This is the variable-demography generalisation of ``PhiD`` on the
``betabinom`` branch; inference interpolates the two logs separately so mutation-
age marginalisation remains a ratio of integrals, as required by that branch.
"""

from __future__ import annotations

import argparse
import json
from math import comb
from pathlib import Path

import mpmath as mp
import numpy as np

from precompute_freq_trajectory_moments import ExactMomentEngine, load_demography


def build(args):
    _tau_of_t, _ne_of_t, (left, right, ne) = load_demography(
        args.ne, series=args.ne_series)
    u = np.geomspace(args.u_min, args.u_max, args.n_u)
    n = args.n_panel
    shape = (len(ne), len(u), n + 1, n + 1)
    log_num = np.full(shape, np.nan, dtype=np.float32)
    log_den = np.full(shape, np.nan, dtype=np.float32)

    # One sufficiently wide engine avoids rebuilding exact rational coefficient
    # tables for every small change in required precision.  At u_min=1e-9 and
    # n=26, required_dps is about 320; 360 retains a generous safety margin.
    engine = ExactMomentEngine(n, dps=args.dps)
    for iw, nev in enumerate(ne):
        x0 = 1.0 / (2.0 * float(nev))
        for iu, uv in enumerate(u):
            needed = engine.required_dps(float(uv))
            if needed > engine.dps:
                raise SystemExit(
                    f"--dps={engine.dps} is insufficient at u={uv:g}; need {needed}")
            moments = engine.moms(float(uv), x0)
            with mp.workdps(engine.dps):
                for nt in range(2, n + 1):
                    for k in range(1, nt):
                        de = mp.fsum(
                            mp.mpf(comb(nt - k, j)) * (-1) ** j * moments[k + j]
                            for j in range(nt - k + 1)
                        )
                        nu = mp.fsum(
                            mp.mpf(comb(nt - k, j)) * (-1) ** j * moments[k + 1 + j]
                            for j in range(nt - k + 1)
                        )
                        if de <= 0 or nu <= 0 or nu > de or not (
                                mp.isfinite(de) and mp.isfinite(nu)):
                            raise SystemExit(
                                f"invalid moments window={iw} u={uv:g} k={k} nT={nt}: "
                                f"num={nu} den={de}")
                        log_num[iw, iu, k, nt] = float(mp.log(nu))
                        log_den[iw, iu, k, nt] = float(mp.log(de))
        print(f"window {iw + 1}/{len(ne)} Ne={nev:.8g}", flush=True)

    args.output.mkdir(parents=True, exist_ok=False)
    np.save(args.output / "u.npy", u)
    np.save(args.output / "window_left.npy", left)
    np.save(args.output / "window_right.npy", right)
    np.save(args.output / "window_ne.npy", ne)
    np.save(args.output / "log_num.npy", log_num)
    np.save(args.output / "log_den.npy", log_den)
    (args.output / "metadata.json").write_text(json.dumps({
        "model": "ARG-conditioned beta-binomial",
        "demography": "piecewise variable Ne via diffusion-time change",
        "ne_file": str(args.ne.resolve()),
        "ne_series": args.ne_series,
        "n_panel": n,
        "u_min": args.u_min,
        "u_max": args.u_max,
        "n_u": args.n_u,
        "dps": args.dps,
        "n_ne_windows": len(ne),
    }, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ne", type=Path, required=True)
    parser.add_argument("--ne-series", default="posterior_mean")
    parser.add_argument("--n-panel", type=int, default=26)
    parser.add_argument("--u-min", type=float, default=1e-9)
    parser.add_argument("--u-max", type=float, default=3.0)
    parser.add_argument("--n-u", type=int, default=80)
    parser.add_argument("--dps", type=int, default=360)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.u_min < args.u_max:
        parser.error("require 0 < --u-min < --u-max")
    build(args)


if __name__ == "__main__":
    main()
