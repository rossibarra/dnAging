#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
posterior_sample_age_infer.py
=============================

Inference stage: two-sided age posteriors for MANY ancient samples (one multi-
sample VCF), using the precomputed age-conditioned frequency table
E[p_T | d0, t_i] (see precompute_freq_trajectory_moments.py).

For each ascertained site and each ARG posterior draw we need only two numbers,
both cheap:
  * t_i : the mutation age, from the interval store's [below, above] for that draw
  * d0  : the present derived-allele count in the n=26 ARG panel (from the panel
          VCF), oriented by the per-draw polarity

We look up phi = E[p_T | d0, t_i] (interpolated over the age grid and averaged
over the store's age interval), track it as the ALT-allele frequency so polarity
flips are handled consistently, average over draws -> phibar_alt(T), then form the
two-sided likelihood for every sample. No D(T)/L(T), no tree walk.

    P(observe ALT carried)  = (1-eps) phibar + eps (1-phibar)
    P(observe hom REF)      = (1-eps)(1-phibar) + eps phibar
    log p(T | D_s) = log p(T) + sum_i log l_{s,i}(T)

ARG draws are averaged into phibar per site; chromosomes are independent given T,
so run one chromosome per invocation and combine chromosomes by summing per-sample
log-marginals (--merge).

The T grid is taken from the frequency table (so precompute and inference share
one grid). Only sites present+called in the ancient VCF contribute for a sample;
the (1-phibar) term only on confident homozygous-REF calls. Sites monomorphic in
the 26 panel (d0 = 0 or n) are skipped (no informative trajectory).

Outputs (into --output/)
------------------------
    samples.txt, grid.npy, ll_marginal.npy (N_samples x grid),
    ages_table.tsv (per-sample MAP/mean/median/95% CI), run.json, cohort_ages.png

============================ REPO ADAPTER LAYER ============================
Only the functions in the ADAPTER section touch normalize_tes / VcfChunk / tskit-
free store internals. Verify those names against your repo; the model math below
is self-contained.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_BASE = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}
_MISS = 255


# =============================================================================
# ADAPTER LAYER  (verify against the repo)
# =============================================================================


def _import_repo():
    try:
        from normalize_tes.snp_age_store import open_snp_age_store
        from normalize_tes.build_draw_polarity import NO_CALL, open_draw_polarity
        from normalize_tes.individual_age_spectrum import read_vcf_chunks
    except Exception as exc:  # pragma: no cover
        raise SystemExit("Could not import normalize_tes; run inside the repo/env. "
                         f"({exc!r})")
    return open_snp_age_store, open_draw_polarity, read_vcf_chunks, NO_CALL


def _chunk_sites(chunk):
    chrom = np.asarray(getattr(chunk, "chromosomes", getattr(chunk, "chrom")))
    pos = np.asarray(getattr(chunk, "positions"), dtype=np.int64)
    rb = np.array([_BASE.get(str(a), _MISS) for a in np.asarray(getattr(chunk, "ref"))],
                  dtype=np.int16)
    ab = np.array([_BASE.get(str(a), _MISS) for a in np.asarray(getattr(chunk, "alt"))],
                  dtype=np.int16)
    return chrom, pos, rb, ab


def _chunk_codes(chunk, names, want):
    """Return the packed genotype codes (rows aligned to `want`) as (len(want) x sites)."""
    codes = np.asarray(getattr(chunk, "codes"))
    names = list(names)
    idx = np.array([names.index(s) for s in want], dtype=np.int64)
    return codes[idx, :] if codes.shape[0] == len(names) else codes[:, idx].T


def _resolve_rows(store, chrom, pos):
    from normalize_tes.snp_position_resolution import resolve_native_position_requests
    res = resolve_native_position_requests(store, np.asarray(chrom).astype(str),
                                           np.asarray(pos, dtype=np.int64), policy="skip")
    rows = np.asarray(getattr(res, "rows", getattr(res, "row_indices")), dtype=np.int64)
    mask = getattr(res, "eligible", None)
    if mask is not None:
        out = np.full(len(pos), -1, dtype=np.int64)
        out[np.asarray(mask, bool)] = rows
        return out
    return rows


def _ancestral_per_draw(polarity_table, row, n_draws):
    r = np.asarray(polarity_table[row]).astype(np.int16)
    if r.shape[0] != n_draws:
        fixed = np.full(n_draws, _MISS, dtype=np.int16)
        fixed[:min(n_draws, r.shape[0])] = r[:n_draws]
        r = fixed
    return r


def _row_intervals(store, row):
    b = store.intervals(np.asarray([row], dtype=np.int64))
    return (np.asarray(b.below, float), np.asarray(b.above, float),
            np.asarray(b.draw_id, np.int64))


# =============================================================================
# Frequency table + lookup
# =============================================================================


def load_table(path):
    d = np.load(path, allow_pickle=True)
    return {"table": d["table"], "d0": d["d0"], "age": d["age"],
            "Tgrid": d["Tgrid"], "n_sample": int(d["n_sample"])}


def phi_lookup(tab, d0, t_lo, t_hi, n_quad=16):
    """E[p_T | d0, t_i] as a T-grid vector, marginalised over the mutation age t_i.

    Under the infinite-sites model the mutation age is UNIFORM on its branch
    interval [t_lo, t_hi] (the store's [below, above]); we therefore average the
    table uniformly in time:  (1/(t_hi-t_lo)) * integral_{t_lo}^{t_hi} E[p_T|d0,t] dt,
    by trapezoidal quadrature at n_quad linearly-spaced ages. The table's age axis
    is log-spaced, so each node is interpolated in log-age. A degenerate branch
    (t_lo == t_hi) collapses to a single point.
    """
    if d0 < 1 or d0 > tab["n_sample"] - 1:      # monomorphic in panel -> no info
        return None
    age = tab["age"]
    la = np.log(np.clip(age, 1e-9, None))
    T = tab["table"][d0 - 1]                     # (n_age, n_T)

    def row_at(a):                               # table row at age a, log-age interp
        k = np.interp(np.log(max(a, 1e-9)), la, np.arange(len(age)))
        k0 = int(np.floor(k)); k1 = min(k0 + 1, len(age) - 1); w = k - k0
        return np.nan_to_num((1 - w) * T[k0] + w * T[k1], nan=0.0)

    lo = max(float(min(t_lo, t_hi)), age[0])
    hi = min(max(float(max(t_lo, t_hi)), lo), age[-1])
    if hi <= lo:                                 # point age (degenerate branch)
        return row_at(lo)
    nodes = np.linspace(lo, hi, n_quad)          # UNIFORM in time along the branch
    wts = np.full(n_quad, 1.0); wts[0] = wts[-1] = 0.5   # trapezoidal weights
    acc = np.zeros(T.shape[1])
    for a, wt in zip(nodes, wts):
        acc += wt * row_at(a)
    return acc / wts.sum()                        # = trapezoidal average over [lo,hi]


# =============================================================================
# VCF reading
# =============================================================================


def read_ancient(vcf_paths, samples, chrom, include, chunk_records, quiet):
    """sample_order, calls: pos -> (rb, ab, carries_alt[bool N], present[bool N])."""
    read_vcf_chunks = _import_repo()[2]
    order = None; calls = {}
    for vcf in vcf_paths:
        for names, chunk, _c, _d in read_vcf_chunks(
                Path(vcf), sample_filter=(list(samples) if samples else None),
                chunk_records=chunk_records, multiallelic="skip", progress=not quiet):
            names = list(names)
            want = list(samples) if samples else names
            if order is None:
                order = want
            chrom_a, pos, rb, ab = _chunk_sites(chunk)
            sub = _chunk_codes(chunk, names, order)
            alt = (sub >> 4).astype(np.int8); called = (sub & 15).astype(np.int8)
            for j in range(len(pos)):
                if chrom is not None and str(chrom_a[j]) != str(chrom):
                    continue
                if include is not None and int(pos[j]) not in include:
                    continue
                if rb[j] == _MISS or ab[j] == _MISS:
                    continue
                # store ALT dosage (0/1/2) and ploidy called (0/1/2): the
                # diploid Hardy-Weinberg genotype likelihood weights by allele
                # counts, and this also covers pseudo-haploid calls (called==1).
                calls[int(pos[j])] = (int(rb[j]), int(ab[j]),
                                      alt[:, j].astype(np.int8).copy(),
                                      called[:, j].astype(np.int8).copy())
    return (order or []), calls


def read_panel_alt(vcf_paths, chrom, chunk_records, quiet, n_expected):
    """pos -> ALT allele count among the panel (require n_expected called haplotypes)."""
    read_vcf_chunks = _import_repo()[2]
    c_alt = {}
    for vcf in vcf_paths:
        for names, chunk, _c, _d in read_vcf_chunks(
                Path(vcf), sample_filter=None, chunk_records=chunk_records,
                multiallelic="skip", progress=not quiet):
            _chrom, pos, rb, ab = _chunk_sites(chunk)
            codes = np.asarray(getattr(chunk, "codes"))
            # sum ALT alleles and called alleles across all panel samples
            if codes.shape[0] == len(list(names)):
                alt = (codes >> 4).astype(np.int64); cal = (codes & 15).astype(np.int64)
                tot_alt = alt.sum(axis=0); tot_called = cal.sum(axis=0)
            else:
                alt = (codes >> 4).astype(np.int64); cal = (codes & 15).astype(np.int64)
                tot_alt = alt.sum(axis=1); tot_called = cal.sum(axis=1)
            for j in range(len(pos)):
                if chrom is not None and str(_chrom[j]) != str(chrom):
                    continue
                if int(tot_called[j]) == n_expected:      # full panel called
                    c_alt[int(pos[j])] = int(tot_alt[j])
    return c_alt


# =============================================================================
# Inference
# =============================================================================


def run_chromosome(args, tab):
    grid = tab["Tgrid"]; n = tab["n_sample"]
    include = None
    if args.include_positions is not None:
        include = set()
        for ln in Path(args.include_positions).read_text().splitlines():
            if ln.strip():
                p = ln.split(); include.add(int(p[1]) if len(p) > 1 else int(p[0]))
    samples = None
    if args.samples_file:
        samples = [s.strip() for s in Path(args.samples_file).read_text().splitlines() if s.strip()]

    order, calls = read_ancient(args.vcf, samples, args.chrom, include,
                                args.chunk_records, args.quiet)
    if not calls:
        raise SystemExit("No called ancient sites matched --chrom/--include.")
    c_alt = read_panel_alt(args.panel_vcf, args.chrom, args.chunk_records,
                           args.quiet, n)

    open_store, open_pol, _rv, _nc = _import_repo()
    store = open_store(args.store); n_draws = int(store.n_posterior_draws)
    polarity, _pm = open_pol(args.draw_polarity, store)

    N = len(order); G = len(grid)
    ll = np.zeros((N, G))
    stats = {"n_samples": N, "sites_used": 0, "sites_no_panel": 0,
             "sites_monomorphic": 0, "chrom": args.chrom}

    # resolve store rows for all ancient positions
    positions = np.array(sorted(calls), dtype=np.int64)
    chroms = np.array([args.chrom] * len(positions))
    rows = _resolve_rows(store, chroms, positions)

    for pos, row in zip(positions, rows):
        pos = int(pos); row = int(row)
        if row < 0:
            continue
        if pos not in c_alt:
            stats["sites_no_panel"] += 1; continue
        ca = c_alt[pos]
        if ca <= 0 or ca >= n:
            stats["sites_monomorphic"] += 1; continue
        rb, ab, alt_ct, cl = calls[pos]

        below, above, draw_id = _row_intervals(store, row)
        anc = _ancestral_per_draw(polarity, row, n_draws)

        phi_sum = np.zeros(G); cnt = 0
        for d in np.unique(draw_id):
            a = int(anc[d]) if d < len(anc) else _MISS
            if a == _MISS or a not in (rb, ab):
                continue
            m = draw_id == d
            t_lo = float(below[m].min()); t_hi = float(above[m].max())
            if a == rb:                       # ALT is derived; d0 = ALT count
                phi = phi_lookup(tab, ca, t_lo, t_hi)
            else:                              # REF is derived; ALT freq = 1 - E[p_T|c_ref]
                phi = phi_lookup(tab, n - ca, t_lo, t_hi)
                if phi is not None:
                    phi = 1.0 - phi
            if phi is None:
                continue
            phi_sum += np.clip(phi, 0.0, 1.0); cnt += 1
        if cnt == 0:
            continue
        phibar = phi_sum / cnt
        eps = args.epsilon
        # per-ALT-allele observation probability q_alt = (1-eps) phibar + eps (1-phibar)
        qA = np.clip((1 - eps) * phibar + eps * (1 - phibar), 1e-300, 1.0)
        logA = np.log(qA); logR = np.log(np.clip(1.0 - qA, 1e-300, 1.0))
        # allele dosage a and ploidy c per sample; log-lik = a logqA + (c-a) log(1-qA)
        if args.ploidy == 1:
            # haploid / pseudo-haploid: one allele per called site (collapse hom calls)
            a_eff = (alt_ct >= 1).astype(np.float64)
            c_eff = (cl >= 1).astype(np.float64)
        else:
            # true diploid genotypes: observed ALT dosage (0/1/2) under Hardy-Weinberg
            a_eff = alt_ct.astype(np.float64)
            c_eff = cl.astype(np.float64)
        ll += np.outer(a_eff, logA) + np.outer(c_eff - a_eff, logR)
        stats["sites_used"] += 1

    return order, grid, ll, stats


def load_prior(args, grid):
    if args.prior_file:
        d = np.loadtxt(args.prior_file); tp, pp = d[:, 0], np.clip(d[:, 1], 1e-300, None)
        return np.log(np.interp(grid, tp, pp, left=pp[0], right=pp[-1]))
    return np.zeros_like(grid)


def merge(args, tab):
    grid = tab["Tgrid"]
    order = (Path(args.merge[0]) / "samples.txt").read_text().split()
    for p in args.merge[1:]:
        if (Path(p) / "samples.txt").read_text().split() != order:
            raise SystemExit("sample order differs across parts.")
    ll = np.sum([np.load(Path(p) / "ll_marginal.npy") for p in args.merge], axis=0)
    return order, grid, ll, {"mode": "merge-chroms", "n_parts": len(args.merge)}


def _trapz(y, x):
    return (getattr(np, "trapezoid", None) or getattr(np, "trapz"))(y, x)


def summarize(grid, lp):
    w = np.exp(lp - lp.max()); area = _trapz(w, grid); dens = w / area if area > 0 else w
    cdf = np.concatenate([[0.0], np.cumsum((dens[:-1] + dens[1:]) / 2 * np.diff(grid))])
    cdf /= cdf[-1] if cdf[-1] > 0 else 1.0
    q = lambda p: float(np.interp(p, cdf, grid))
    return {"map_T": float(grid[int(np.argmax(lp))]), "mean_T": float(_trapz(dens * grid, grid)),
            "median_T": q(.5), "ci95_lower_T": q(.025), "ci95_upper_T": q(.975)}, dens


def write_outputs(outdir, order, grid, log_prior, ll, stats, args):
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "samples.txt").write_text("\n".join(order) + "\n")
    np.save(outdir / "grid.npy", grid); np.save(outdir / "ll_marginal.npy", ll)
    lp = log_prior[None, :] + ll
    rows = []
    per = outdir / "posterior"
    if args.per_sample_tsv:
        per.mkdir(exist_ok=True)
    for i, s in enumerate(order):
        m, dens = summarize(grid, lp[i]); rows.append((s, m))
        if args.per_sample_tsv:
            with (per / f"{s}.tsv").open("w") as fh:
                fh.write("T\tlog_posterior\tdensity\n")
                for k in range(len(grid)):
                    fh.write(f"{grid[k]:.6g}\t{lp[i][k]:.6g}\t{dens[k]:.6g}\n")
    with (outdir / "ages_table.tsv").open("w") as fh:
        fh.write("sample\tmap_T\tmean_T\tmedian_T\tci95_lower_T\tci95_upper_T\n")
        for s, m in rows:
            fh.write(f"{s}\t{m['map_T']:.6g}\t{m['mean_T']:.6g}\t{m['median_T']:.6g}\t"
                     f"{m['ci95_lower_T']:.6g}\t{m['ci95_upper_T']:.6g}\n")
    (outdir / "run.json").write_text(json.dumps({"counts": stats,
        "settings": {"epsilon": args.epsilon, "chrom": args.chrom}}, indent=2))
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        maps = np.array([m["map_T"] for _, m in rows])
        fig, ax = plt.subplots(figsize=(7, 4)); ax.hist(maps, bins=40)
        ax.set_xlabel("per-sample MAP age T (ARG generations)"); ax.set_ylabel("samples")
        ax.set_title(f"cohort age distribution (N={len(rows)})")
        fig.tight_layout(); fig.savefig(outdir / "cohort_ages.png", dpi=150); plt.close(fig)
    except Exception:
        pass


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="posterior_sample_age_infer.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Age posteriors for many ancient samples via the precomputed "
                    "E[p_T|d0,t_i] table (no tree walk, no D/L). One chromosome/run.")
    io = p.add_argument_group("inputs / outputs")
    io.add_argument("--freq-table", type=Path, required=True,
                    help=".npz from precompute_freq_trajectory_moments.py")
    io.add_argument("--store", type=Path, help="interval store (t_i per site/draw).")
    io.add_argument("--draw-polarity", type=Path, help="per-draw polarity table.")
    io.add_argument("--panel-vcf", type=Path, nargs="+", default=[],
                    help="VCF of the 26 ARG-panel haplotypes (gives d0 = ALT count).")
    io.add_argument("--vcf", type=Path, nargs="+", default=[],
                    help="the multi-sample ancient VCF(s).")
    io.add_argument("--samples-file", type=Path, help="optional subset of samples.")
    io.add_argument("--chrom", type=str, help="chromosome label (matches the VCFs).")
    io.add_argument("--output", type=Path, required=True)
    io.add_argument("--include-positions", type=Path,
                    help="'chrom pos' site list (e.g. an approximately-neutral set).")
    io.add_argument("--per-sample-tsv", action="store_true")
    p.add_argument("--ploidy", type=int, choices=(1, 2), default=1,
                   help="ploidy of the ANCIENT-sample genotypes: 1 = haploid / "
                        "pseudo-haploid (one allele per called site; homozygous "
                        "calls collapsed to one observation) [default]; 2 = true "
                        "diploid genotypes (ALT dosage 0/1/2 used, Hardy-Weinberg). "
                        "Use 1 for pseudo-haploid aDNA even if written as hom "
                        "diploid, or 2 would double-count each site.")
    p.add_argument("--epsilon", type=float, default=1e-3,
                   help="per-ALLELE genotype-error probability [1e-3].")
    p.add_argument("--prior-file", type=Path, default=None)
    p.add_argument("--chunk-records", type=int, default=20000)
    p.add_argument("--merge", type=Path, nargs="+", default=None,
                   help="sum per-sample marginals across chromosome parts.")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)
    if args.merge is None:
        need = [n for n in ("store", "draw_polarity", "chrom") if getattr(args, n) is None]
        if not args.vcf or not args.panel_vcf or need:
            p.error("a run needs --store, --draw-polarity, --panel-vcf, --vcf and "
                    "--chrom (or use --merge).")
    return args


def main(argv=None):
    args = parse_args(argv)
    tab = load_table(args.freq_table)
    if args.merge is not None:
        order, grid, ll, stats = merge(args, tab)
    else:
        order, grid, ll, stats = run_chromosome(args, tab)
    log_prior = load_prior(args, grid)
    write_outputs(args.output, order, grid, log_prior, ll, stats, args)
    if not args.quiet:
        print(f"[infer] {len(order)} samples chrom={args.chrom or stats.get('mode')} "
              f"sites_used={stats.get('sites_used','-')} -> {args.output}/ages_table.tsv",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
