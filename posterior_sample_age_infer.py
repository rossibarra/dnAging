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
two-sided likelihood for every sample. No D(T)/L(T), no tree walk. With r the
per-allele probability of OBSERVING ALT, r = eps + (1-2eps) p_T,

    haploid  (--ploidy 1): P(ALT) = E[r],  P(REF) = 1 - E[r]
    diploid  (--ploidy 2): P(2) = E[r^2], P(1) = 2(E[r]-E[r^2]), P(0) = 1-2E[r]+E[r^2]
    log p(T | D_s) = log p(T) + sum_i log l_{s,i}(T)

P(dosage) is quadratic in r, so --ploidy 2 also needs the conditional SECOND moment
E[p_T^2 | d0, t_i] (the table's 'table2' plane): E[p_T^2] != E[p_T]^2.

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


def _attr(obj, *names):
    """First attribute of `obj` present among `names`.

    A nested getattr(obj, "a", getattr(obj, "b")) evaluates the fallback EAGERLY,
    so it raises even when "a" exists -- exactly the case the fallback is for.
    """
    for nm in names:
        if hasattr(obj, nm):
            return getattr(obj, nm)
    raise AttributeError(f"{type(obj).__name__} has none of the attributes "
                         f"{', '.join(map(repr, names))}")


def _chunk_sites(chunk):
    chrom = np.asarray(_attr(chunk, "chromosomes", "chrom"))
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
    rows = np.asarray(_attr(res, "rows", "row_indices"), dtype=np.int64)
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


def _is_multiply_mapped(draw_id):
    """Whether any ARG draw supplies multiple branch intervals for one site."""
    _draws, counts = np.unique(np.asarray(draw_id, np.int64), return_counts=True)
    return bool(np.any(counts > 1))


# =============================================================================
# Frequency table + lookup
# =============================================================================


def load_table(path):
    d = np.load(path, allow_pickle=True)
    t = {"table": d["table"], "d0": d["d0"], "age": d["age"],
         "Tgrid": d["Tgrid"], "n_sample": int(d["n_sample"])}
    t["n_panel"] = (d["n_panel"] if "n_panel" in d.files else
                    np.array([t["n_sample"]], dtype=np.int64))
    if "age_tau" in d.files:
        t["age_tau"] = d["age_tau"]
    if "table2" in d.files:                      # E[p_T^2|d0,t_i]: diploid only
        t["table2"] = d["table2"]
    return t


def phi_lookup(tab, d0, t_lo, t_hi, n_quad=16, key="table", n_called=None):
    """E[p_T | d0, t_i] as a T-grid vector, marginalised over the mutation age t_i.

    Under the infinite-sites model the mutation age is UNIFORM on its branch
    interval [t_lo, t_hi] (the store's [below, above]); we therefore average the
    table uniformly in time:  (1/(t_hi-t_lo)) * integral_{t_lo}^{t_hi} E[p_T|d0,t] dt,
    by trapezoidal quadrature at n_quad linearly-spaced ages. The table's age axis
    is log-spaced, so each node is interpolated in log-age. A degenerate branch
    (t_lo == t_hi) collapses to a single point.

    A branch reaching below the table's youngest age (the store's `below` can be 0)
    is integrated only over the covered part but normalised by the TRUE branch
    length, which is exact for T >= age[0]; see the comment at the renormalisation.

    key="table2" reads the second-moment plane E[p_T^2 | d0, t_i] instead; the
    branch average is linear in the tabulated quantity, so the same quadrature
    applies. (Nonzero values for T inside [t_lo, t_hi] are correct: the mutation
    age is uncertain within the branch.)
    """
    n_called = tab["n_sample"] if n_called is None else int(n_called)
    if d0 < 1 or d0 > n_called - 1:             # monomorphic in called panel -> no info
        return None
    age = tab["age"]; Tg = tab["Tgrid"]
    la = np.log(np.clip(age, 1e-9, None))
    plane = tab[key]
    if plane.ndim == 4:
        matches = np.flatnonzero(np.asarray(tab["n_panel"]) == n_called)
        if matches.size != 1:
            return None
        T = plane[int(matches[0]), d0 - 1]
    else:                                        # legacy fixed-n table
        if n_called != tab["n_sample"]:
            return None
        T = plane[d0 - 1]                        # (n_age, n_T)

    def row_at(a):                               # table row at age a, log-age interp
        k = np.interp(np.log(max(a, 1e-9)), la, np.arange(len(age)))
        k0 = int(np.floor(k)); k1 = min(k0 + 1, len(age) - 1); w = k - k0
        r = (1 - w) * T[k0] + w * T[k1]
        # blending a T>=t_i zero row with a nonzero one leaks probability across the
        # mutation-existence boundary, so re-impose p_T = 0 at the interpolated age
        return np.where(Tg >= a, 0.0, r)

    b_lo = float(min(t_lo, t_hi)); b_hi = float(max(t_lo, t_hi))   # true branch
    lo = max(b_lo, age[0])
    hi = min(max(b_hi, lo), age[-1])
    if hi <= lo:                                 # point age, or branch wholly below
        return row_at(lo)                        # age[0] (row_at zeroes T >= age[0])
    nodes = np.linspace(lo, hi, n_quad)          # UNIFORM in time along the branch
    wts = np.full(n_quad, 1.0); wts[0] = wts[-1] = 0.5   # trapezoidal weights
    acc = np.zeros(T.shape[1])
    for a, wt in zip(nodes, wts):
        acc += wt * row_at(a)
    avg = acc / wts.sum()                        # mean over the COVERED part [lo,hi]
    # Renormalise onto the TRUE branch length. Mutation ages in the uncovered
    # [b_lo, age[0]) are younger than any sample age T >= age[0], so they contribute
    # p_T = 0 exactly: the covered integral is already the whole numerator, and
    # dividing by (hi-lo) instead of (b_hi-b_lo) inflates the site by
    # (b_hi-b_lo)/(hi-lo). This is EXACT for T >= age[0] only; below age[0] the
    # uncovered ages in (T, age[0]) do contribute and this underestimates.
    # ASYMMETRIC ON PURPOSE -- do NOT mirror it at the upper end, where uncovered
    # ages are OLDER than the sample and genuinely contribute; that end is handled
    # by capping t_hi at --mutation-age-max before we are called.
    if lo > b_lo and hi >= b_hi and b_hi > b_lo:
        avg = avg * (hi - lo) / (b_hi - b_lo)
    return avg


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
    """pos -> (ALT count, called count, ref code, alt code) among the panel.

    The ref/alt codes are kept so inference can check the panel's orientation
    against the ancient VCF's: joining on position alone would silently read the
    count as n-count wherever the two files disagree on which allele is REF.
    """
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
                nc = int(tot_called[j])
                if 0 < nc <= n_expected:
                    c_alt[int(pos[j])] = (int(tot_alt[j]), nc,
                                          int(rb[j]), int(ab[j]))
    return c_alt


# =============================================================================
# Inference
# =============================================================================


def run_chromosome(args, tab):
    grid = tab["Tgrid"]; n = tab["n_sample"]
    if "age_tau" not in tab:
        raise SystemExit("--freq-table lacks the age_tau axis required by the default "
                         "diffusion-time mutation cutoff; rebuild it with the current "
                         "precompute_freq_trajectory_moments.py")
    age_tau = np.asarray(tab["age_tau"], float)
    if age_tau[-1] <= args.mutation_age_max:
        raise SystemExit(f"frequency table reaches tau={age_tau[-1]:.6g}, which does "
                         f"not extend beyond --mutation-age-max={args.mutation_age_max:g}")
    mutation_age_max_generations = float(np.interp(
        args.mutation_age_max, age_tau, np.asarray(tab["age"], float)))
    age_min_generations = float(np.asarray(tab["age"], float)[0])
    available_n = set(int(v) for v in np.asarray(tab["n_panel"]))
    needed_n = set(range(args.min_n, n + 1))
    if not needed_n.issubset(available_n):
        missing = sorted(needed_n - available_n)
        raise SystemExit(f"--freq-table lacks called-panel sizes required by "
                         f"--min-n={args.min_n}: {missing}; rebuild the table")
    need2 = args.ploidy == 2                     # diploid genotypes need E[p_T^2]
    if need2 and "table2" not in tab:
        raise SystemExit("--ploidy 2 needs the conditional second-moment plane "
                         "'table2', absent from --freq-table (pre-format-2 table). "
                         "Rebuild it with the current "
                         "precompute_freq_trajectory_moments.py; the squared mean is "
                         "NOT a valid substitute (E[p^2] != E[p]^2).")
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
             "sites_monomorphic": 0, "sites_allele_mismatch": 0,
             "sites_age_filtered": 0, "sites_numerical_failure": 0,
             "sites_multiple_mapped": 0, "sites_age_clipped_low": 0,
             "sites_panel_below_min_n": 0,
             "chrom": args.chrom}

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
        ca, n_called, p_rb, p_ab = c_alt[pos]
        if n_called < args.min_n:
            stats["sites_panel_below_min_n"] += 1
            continue
        rb, ab, alt_ct, cl = calls[pos]
        # harmonise the panel ALT count to the ANCIENT VCF's REF/ALT orientation
        # BEFORE testing it: the two files are joined on position alone
        if _MISS in (p_rb, p_ab, rb, ab):        # unknown base code on either side
            stats["sites_allele_mismatch"] += 1; continue
        if p_rb == rb and p_ab == ab:            # same orientation
            pass
        elif p_rb == ab and p_ab == rb:          # REF/ALT swapped between the VCFs
            ca = n_called - ca
        else:                                    # different alleles, or an unknown base
            stats["sites_allele_mismatch"] += 1; continue
        if ca <= 0 or ca >= n_called:
            stats["sites_monomorphic"] += 1; continue

        below, above, draw_id = _row_intervals(store, row)
        if _is_multiply_mapped(draw_id):
            stats["sites_multiple_mapped"] += 1
            continue
        anc = _ancestral_per_draw(polarity, row, n_draws)

        phi_sum = np.zeros(G); phi2_sum = np.zeros(G); cnt = 0
        age_rejected = numerical_rejected = low_clipped = False
        for d in np.unique(draw_id):
            a = int(anc[d]) if d < len(anc) else _MISS
            if a == _MISS or a not in (rb, ab):
                continue
            m = draw_id == d
            t_lo = float(below[m].min()); t_hi = float(above[m].max())
            if t_lo >= mutation_age_max_generations:
                age_rejected = True
                continue
            if t_hi > mutation_age_max_generations:
                t_hi = mutation_age_max_generations
            if t_lo < age_min_generations:       # branch reaches below the table
                low_clipped = True
            d0 = ca if a == rb else n_called - ca
            phi = phi_lookup(tab, d0, t_lo, t_hi, n_called=n_called)
            if phi is None:
                continue
            phi2 = (phi_lookup(tab, d0, t_lo, t_hi, key="table2", n_called=n_called)
                    if need2 else None)
            if not np.all(np.isfinite(phi)) or (phi2 is not None and
                                                not np.all(np.isfinite(phi2))):
                numerical_rejected = True
                continue
            if a != rb:                        # REF is derived; ALT freq = 1 - E[p_T|c_ref]
                if phi2 is not None:           # E[(1-X)^2] = 1 - 2 E[X] + E[X^2]
                    phi2 = 1.0 - 2.0 * phi + phi2
                phi = 1.0 - phi
            phi_sum += np.clip(phi, 0.0, 1.0); cnt += 1
            if phi2 is not None:
                phi2_sum += np.clip(phi2, 0.0, 1.0)
        if cnt == 0:
            if age_rejected:
                stats["sites_age_filtered"] += 1
            if numerical_rejected:
                stats["sites_numerical_failure"] += 1
            continue
        if low_clipped:            # site IS used; renormalised onto the true branch
            stats["sites_age_clipped_low"] += 1
        phibar = phi_sum / cnt
        eps = args.epsilon
        # per-ALT-allele observation probability r = (1-eps) X + eps (1-X); averaging
        # it over the ALT frequency X needs only the first moment: E[r] = eps + (1-2eps) E[X]
        qA = np.clip((1 - eps) * phibar + eps * (1 - phibar), 1e-300, 1.0)
        logA = np.log(qA); logR = np.log(np.clip(1.0 - qA, 1e-300, 1.0))
        if args.ploidy == 1:
            # haploid / pseudo-haploid: ONE Bernoulli(E[r]) per called site (collapse
            # hom calls); log-lik = a logqA + (c-a) log(1-qA) with a,c in {0,1}
            a_eff = (alt_ct >= 1).astype(np.float64)
            c_eff = (cl >= 1).astype(np.float64)
            ll += np.outer(a_eff, logA) + np.outer(c_eff - a_eff, logR)
        else:
            # true diploid genotypes: the two alleles are iid Bernoulli(r) only GIVEN
            # the latent frequency X, and P(dosage) is QUADRATIC in r, so plugging the
            # mean in would be wrong (E[r^2] != E[r]^2); use the second moment.
            phi2bar = phi2_sum / cnt
            Er2 = np.clip(eps ** 2 + 2 * eps * (1 - 2 * eps) * phibar
                          + (1 - 2 * eps) ** 2 * phi2bar, 0.0, 1.0)
            lg = np.log(np.clip(np.stack([1.0 - 2.0 * qA + Er2,   # P(dosage=0) = E[(1-r)^2]
                                          2.0 * (qA - Er2),       # P(dosage=1) = 2E[r(1-r)]
                                          Er2]), 1e-300, 1.0))    # P(dosage=2) = E[r^2]
            dip = cl >= 2                      # both alleles called: full genotype
            ll[dip] += lg[np.clip(alt_ct[dip].astype(np.int64), 0, 2)]
            half = cl == 1                     # partially-called: one Bernoulli(E[r])
            if half.any():
                ll[half] += np.where((alt_ct[half] >= 1)[:, None], logA, logR)
        stats["sites_used"] += 1

    return order, grid, ll, stats


def load_prior(args, grid):
    if args.prior_file:
        try:
            d = np.loadtxt(args.prior_file, ndmin=2)
        except (OSError, ValueError) as exc:
            raise SystemExit(f"Could not read --prior-file {args.prior_file}: {exc}") from exc
        if d.ndim != 2 or d.shape[1] != 2:
            raise SystemExit("--prior-file must contain exactly two columns: age density")
        if d.shape[0] < 2:
            raise SystemExit("--prior-file must contain at least two rows")
        if not np.all(np.isfinite(d)):
            raise SystemExit("--prior-file ages and densities must all be finite")
        tp, pp = d[:, 0], d[:, 1]
        if not np.all(np.diff(tp) > 0):
            raise SystemExit("--prior-file ages must be strictly increasing and unique")
        if np.any(pp < 0):
            raise SystemExit("--prior-file densities must be non-negative")
        if not np.any(pp > 0):
            raise SystemExit("--prior-file must contain at least one positive density")
        pp = np.clip(pp, 1e-300, None)
        return np.log(np.interp(grid, tp, pp, left=pp[0], right=pp[-1]))
    return np.zeros_like(grid)


def merge(args, tab):
    grid = tab["Tgrid"]
    order = (Path(args.merge[0]) / "samples.txt").read_text().split()
    expected_shape = (len(order), len(grid))
    arrays = []
    for p in args.merge:
        part = Path(p)
        if (part / "samples.txt").read_text().split() != order:
            raise SystemExit(f"sample order differs in merge part {part}")
        part_grid = np.load(part / "grid.npy", allow_pickle=False)
        if part_grid.shape != grid.shape or not np.array_equal(part_grid, grid):
            raise SystemExit(f"sample-age grid differs in merge part {part}")
        part_ll = np.load(part / "ll_marginal.npy", allow_pickle=False)
        if part_ll.shape != expected_shape:
            raise SystemExit(f"likelihood shape in merge part {part} is {part_ll.shape}; "
                             f"expected {expected_shape}")
        arrays.append(part_ll)
    ll = np.stack(arrays, axis=0).sum(axis=0)
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
        "settings": {"epsilon": args.epsilon, "chrom": args.chrom,
                     "mutation_age_max": args.mutation_age_max}}, indent=2))
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
                        "diploid genotypes (ALT dosage 0/1/2 used, Hardy-Weinberg; "
                        "needs a freq table carrying the second-moment plane). "
                        "Use 1 for pseudo-haploid aDNA even if written as hom "
                        "diploid, or 2 would double-count each site.")
    p.add_argument("--epsilon", type=float, default=0.01,
                   help="per-ALLELE ancient-VCF genotype-error probability [0.01].")
    p.add_argument("--min-n", type=int, default=20,
                   help="drop panel sites with fewer than this many called haplotypes [20].")
    p.add_argument("--mutation-age-max", type=float, default=3.0,
                   help="maximum mutation age in diffusion units tau; older mutation-"
                        "age mass is discarded and crossing intervals are truncated "
                        "[3.0]. At constant Ne=10000, tau=3 is 60000 generations.")
    p.add_argument("--prior-file", type=Path, default=None)
    p.add_argument("--chunk-records", type=int, default=20000)
    p.add_argument("--merge", type=Path, nargs="+", default=None,
                   help="sum per-sample marginals across chromosome parts.")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)
    if not 0.0 <= args.epsilon < 0.5:
        p.error("--epsilon must satisfy 0 <= epsilon < 0.5")
    if args.min_n < 2:
        p.error("--min-n must be at least 2")
    if args.mutation_age_max <= 0:
        p.error("--mutation-age-max must be positive")
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
