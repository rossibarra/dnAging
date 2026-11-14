"""Non-sharing-ratio dating of ancient haplotypes.

Implements the estimator described in "Dating ancient genomes with inferred ARGs"
(Nov 11 2026).  The ARG is used only to decide which mutations are old.  For a
mutation of true age a > T carried by a reference haplotype r, a target fails to
share it only if its lineage has not coalesced with r's lineage by time a:

    P_mod(not shared | a) = exp(-int_0^a c),   P_anc(not shared | a) = exp(-int_T^a c)

so P_anc / P_mod = exp(int_0^T c) for every a > T.  With f the fraction of
(site, reference) pairs where the reference carries the derived allele and the
target does not,

    H_hat = log(f_anc / mean f_mod),     int_0^T_hat c(t) dt = H_hat.

Everything here is pure numpy (haploid 0/1 genotypes, 1 = derived) so that it
can be unit tested without tskit.  Tree-sequence helpers are in the validation
driver and in `lineage_counts`.

Design rules from the report, which this module enforces by construction:

* modern controls must be held out of the ARG (callers pass them as `targets`,
  never as part of `panel`);
* sites fixed in the panel are included without dating (they are older than the
  panel's local common ancestor); sites fixed in every sample are dropped.
"""
from __future__ import annotations

import numpy as np


def select_sites(panel, targets, ages, a0, include_panel_fixed=True, called=None):
    """Boolean mask of sites used by the estimator.

    panel   : (n_panel, n_sites) 0/1, haplotypes used as references (the ARG panel)
    targets : (n_targets, n_sites) 0/1, ancient samples and held-out moderns
    ages    : (n_sites,) ARG-dated (or true) mutation ages; NaN where unknown.
              Only consulted for sites that segregate in the panel.
    a0      : cutoff; segregating sites dated strictly above a0 are kept
    called  : optional (n_targets, n_sites) bool; False = missing call.  Missing
              targets never count towards "fixed in every sample".
    """
    panel = np.asarray(panel)
    targets = np.asarray(targets)
    ages = np.asarray(ages, dtype=float)
    n_panel = panel.shape[0]
    count = panel.sum(axis=0)
    fixed = count == n_panel
    segregating = (count > 0) & ~fixed
    with np.errstate(invalid="ignore"):
        old = segregating & (ages > a0)
    keep = old | (fixed if include_panel_fixed else False)
    # Sites fixed in every sample add equally to every target's denominator.
    if called is None:
        called = np.ones(targets.shape, dtype=bool)
    called = np.asarray(called, dtype=bool)
    derived_called = (targets == 1) & called
    fixed_everywhere = (fixed & (called.sum(axis=0) > 0)
                        & (derived_called.sum(axis=0) == called.sum(axis=0)))
    return keep & ~fixed_everywhere


def nonsharing_counts(panel, targets, mask, called=None, block_ids=None, n_blocks=None,
                      reference_counts=None):
    """Numerator/denominator of f per block and target.

    num[b, t] = #(site, reference) pairs in block b where the reference carries the
    derived allele and target t is called ancestral; den[b, t] = #pairs where t is
    called at all.  f for any union of blocks is sum(num)/sum(den).
    """
    mask = np.asarray(mask, dtype=bool)
    panel = np.asarray(panel)[:, mask]
    targets = np.asarray(targets)[:, mask]
    called = (np.ones(targets.shape, dtype=bool) if called is None
              else np.asarray(called, dtype=bool)[:, mask])
    ref_derived = (panel.sum(axis=0).astype(float) if reference_counts is None else
                   np.asarray(reference_counts, dtype=float)[mask])
    if ref_derived.shape != (mask.sum(),):
        raise ValueError("reference_counts must have one entry per unmasked site")
    if block_ids is None:
        block_ids, n_blocks = np.zeros(mask.sum(), dtype=int), 1
    else:
        block_ids = np.asarray(block_ids)[mask]
    anc_called = ((targets == 0) & called) * ref_derived
    all_called = called * ref_derived
    num = np.zeros((n_blocks, targets.shape[0]))
    den = np.zeros((n_blocks, targets.shape[0]))
    np.add.at(num, block_ids, anc_called.T)
    np.add.at(den, block_ids, all_called.T)
    return num, den


def nonsharing_fraction(panel, targets, mask, called=None, reference_counts=None):
    """f per target: derived-in-reference pairs where the target is ancestral.

    Returns (f, n_pairs), both of shape (n_targets,).  Each (site, reference) pair
    counts once, so a site carried by k panel haplotypes contributes k pairs.
    """
    num, den = nonsharing_counts(panel, targets, mask, called,
                                 reference_counts=reference_counts)
    num, den = num[0], den[0]
    if np.any(den == 0):
        raise ValueError("no (site, reference) pairs after site selection "
                         "(or a target has no called sites)")
    return num / den, den


def estimate_H(f_anc, f_mod):
    """log(f_anc / mean f_mod).  `f_mod` is an array of control fractions."""
    f_mod = np.asarray(f_mod, dtype=float)
    return np.log(np.asarray(f_anc, dtype=float) / f_mod.mean())


def cumulative_hazard(edges, rates):
    """Cumulative int_0^t c at the epoch boundaries `edges` (len m+1, edges[0]=0).

    The last edge may be inf; the cumulative value there is then inf unless the
    last rate is 0.
    """
    edges = np.asarray(edges, dtype=float)
    rates = np.asarray(rates, dtype=float)
    if edges.shape[0] != rates.shape[0] + 1 or edges[0] != 0:
        raise ValueError("edges must have len(rates)+1 entries starting at 0")
    widths = np.diff(edges)
    with np.errstate(invalid="ignore"):
        steps = np.where(rates == 0, 0.0, rates * widths)
    return np.concatenate([[0.0], np.cumsum(steps)])


def invert_hazard(H, edges, rates):
    """Solve int_0^T c = H for T (piecewise-constant c).

    H <= 0 returns 0.  H beyond the total hazard returns the last finite edge
    (or inf if the last epoch is open with positive rate, which always reaches H).
    """
    edges = np.asarray(edges, dtype=float)
    rates = np.asarray(rates, dtype=float)
    cum = cumulative_hazard(edges, rates)
    H = np.atleast_1d(np.asarray(H, dtype=float))
    out = np.empty_like(H)
    for j, h in enumerate(H):
        if not np.isfinite(h) or h <= 0:
            out[j] = 0.0 if h <= 0 else np.nan
            continue
        i = int(np.searchsorted(cum, h, side="right")) - 1
        i = min(i, rates.shape[0] - 1)
        if rates[i] <= 0:
            out[j] = edges[i + 1]
        else:
            out[j] = edges[i] + (h - cum[i]) / rates[i]
    return out


def rate_from_lineage_counts(times, kappa, pair_counts=None):
    """Pairwise coalescence rate c per epoch from average lineage counts.

    times : (m+1,) increasing epoch boundaries starting at 0
    kappa : (m+1,) average number of lineages at each boundary

    Lineages are lost at rate C(k,2) c(t), so over an epoch
    c ~= (kappa_i - kappa_{i+1}) / (dt * mean C(k, 2)).  `pair_counts`, when
    supplied, is the genome-averaged C(k, 2) at each boundary.  Omitting it uses
    C(kbar, 2), which is exact only when k does not vary among local trees.
    Epochs with no lineage pairs or non-positive loss inherit the previous
    epoch's rate (first epoch: 0).
    """
    times = np.asarray(times, dtype=float)
    kappa = np.asarray(kappa, dtype=float)
    dt = np.diff(times)
    loss = kappa[:-1] - kappa[1:]
    if pair_counts is None:
        kbar = 0.5 * (kappa[:-1] + kappa[1:])
        pairs = kbar * (kbar - 1) / 2.0
    else:
        pair_counts = np.asarray(pair_counts, dtype=float)
        if pair_counts.shape != kappa.shape:
            raise ValueError("pair_counts must have the same shape as kappa")
        pairs = 0.5 * (pair_counts[:-1] + pair_counts[1:])
    c = np.zeros_like(dt)
    last = 0.0
    for i in range(dt.shape[0]):
        if pairs[i] > 0 and loss[i] > 0:
            last = loss[i] / (dt[i] * pairs[i])
        c[i] = last
    return c


def lineage_counts(ts, times):
    """Genome-averaged number of lineages at each time (tskit tree sequence).

    Counts edges with child_time <= t < parent_time, weighted by span.
    """
    tables = ts.tables
    node_time = tables.nodes.time
    edges = tables.edges
    child_t = node_time[edges.child]
    parent_t = node_time[edges.parent]
    span = edges.right - edges.left
    times = np.asarray(times, dtype=float)
    out = np.empty(times.shape[0])
    for i, t in enumerate(times):
        sel = (child_t <= t) & (t < parent_t)
        out[i] = span[sel].sum() / ts.sequence_length
    return out


def lineage_moments(ts, times):
    """Genome-averaged k and C(k, 2) at each time across local trees.

    The second moment is required because mean(C(k, 2)) is not generally
    C(mean(k), 2) when local genealogies have different lineage counts.
    """
    times = np.asarray(times, dtype=float)
    counts = np.zeros(times.shape)
    pairs = np.zeros(times.shape)
    for tree in ts.trees():
        k = np.fromiter((tree.num_lineages(float(t)) for t in times),
                        dtype=float, count=times.size)
        weight = tree.span / ts.sequence_length
        counts += weight * k
        pairs += weight * k * (k - 1) / 2.0
    return counts, pairs


def estimate_ages(panel, ancient, controls, ages, a0, edges, rates,
                  include_panel_fixed=True, called_ancient=None, called_controls=None):
    """End-to-end: returns dict with T_hat per ancient sample and diagnostics.

    ancient  : (n_anc, n_sites) 0/1;  controls : (n_ctrl, n_sites) 0/1 held-out moderns
    """
    ancient = np.asarray(ancient)
    controls = np.asarray(controls)
    targets = np.vstack([ancient, controls])
    called = None
    if called_ancient is not None or called_controls is not None:
        ca = np.ones(ancient.shape, bool) if called_ancient is None else np.asarray(called_ancient, bool)
        cc = np.ones(controls.shape, bool) if called_controls is None else np.asarray(called_controls, bool)
        called = np.vstack([ca, cc])
    mask = select_sites(panel, targets, ages, a0, include_panel_fixed, called)
    f, n_pairs = nonsharing_fraction(panel, targets, mask, called)
    f_anc, f_mod = f[: ancient.shape[0]], f[ancient.shape[0]:]
    H = estimate_H(f_anc, f_mod)
    T = invert_hazard(H, edges, rates)
    return {
        "T_hat": T, "H": H, "f_anc": f_anc, "f_mod": f_mod,
        "n_sites": int(mask.sum()), "n_pairs": float(np.mean(n_pairs)),
    }


def edge_midpoint_ages_slow(ts):
    """Dated age per site: midpoint of the edge above each mutation (oldest if several).

    Sites without mutations are NaN; a mutation on a local root has no edge above it
    and is given +inf (it predates the panel's common ancestor).
    """
    ages = np.full(ts.num_sites, np.nan)
    times = ts.tables.nodes.time
    for tree in ts.trees():
        for site in tree.sites():
            for m in site.mutations:
                p = tree.parent(m.node)
                a = np.inf if p == -1 else 0.5 * (times[m.node] + times[p])
                if np.isnan(ages[site.id]) or a > ages[site.id]:
                    ages[site.id] = a
    return ages


def panel_derived_matrix_slow(ts):
    """(n_samples, n_sites) 0/1 derived indicator for the ARG's samples.

    Sites with more than two alleles are marked -1 for every sample (callers should
    drop them); the derived allele is any non-ancestral state.
    """
    out = np.zeros((ts.num_samples, ts.num_sites), dtype=np.int8)
    for var in ts.variants():
        anc = var.site.ancestral_state
        if len(var.alleles) > 2 or (len(var.alleles) == 2 and anc not in var.alleles):
            out[:, var.site.id] = -1
            continue
        anc_idx = var.alleles.index(anc) if anc in var.alleles else 0
        out[:, var.site.id] = (var.genotypes != anc_idx).astype(np.int8)
    return out


def ages_from_counts(num, den, is_control, edges, rates):
    """T_hat per ancient target from block counts (blocks x targets).

    is_control : (n_targets,) bool; True columns are held-out moderns.
    Returns (T_hat, H, f) for the non-control targets, in column order.
    """
    num = np.asarray(num, dtype=float).sum(axis=0)
    den = np.asarray(den, dtype=float).sum(axis=0)
    is_control = np.asarray(is_control, dtype=bool)
    if not is_control.any() or is_control.all():
        raise ValueError("need at least one control and one ancient target")
    with np.errstate(divide="ignore", invalid="ignore"):
        f = num / den
        H = np.log(f[~is_control] / f[is_control].mean())
    return invert_hazard(H, edges, rates), H, f[~is_control]


FLAG_NONE, FLAG_F_ZERO, FLAG_H_NONPOS, FLAG_NO_DATA = "", "f0", "H<=0", "no_data"


def target_flags(f, H, den=None):
    """Per-target warning labels for the estimator's degenerate cases.

    "f0"      : the target never fails to share a derived reference allele (f = 0), so
                H = -inf and T_hat is reported as 0 only by convention
    "H<=0"    : the target shares at least as much as the controls do (H <= 0), so the
                estimate is truncated to T_hat = 0
    "no_data" : no (site, reference) pairs were called for the target, or the control
                mean f is not positive (H undefined; T_hat is NaN)
    Returns an array of strings, '' for a regular target.
    """
    f = np.asarray(f, dtype=float)
    H = np.asarray(H, dtype=float)
    out = np.full(f.shape, FLAG_NONE, dtype=object)
    with np.errstate(invalid="ignore"):
        out[H <= 0] = FLAG_H_NONPOS
        out[f == 0] = FLAG_F_ZERO
    out[~np.isfinite(f) | (np.isnan(H))] = FLAG_NO_DATA
    if den is not None:
        out[np.asarray(den, dtype=float) == 0] = FLAG_NO_DATA
    return out


def edge_child_ages(ts):
    """Dated age per site: time of the child node the mutation sits on (oldest if several).

    This is the lower bound of the mutation's edge, so `age > A0` means the mutation is
    older than A0 whatever the parent's time.  Sites without mutations are NaN.
    """
    tables = ts.tables
    ages = np.full(tables.sites.num_rows, np.nan)
    if tables.mutations.num_rows == 0:
        return ages
    t = np.asarray(tables.nodes.time)[np.asarray(tables.mutations.node)]
    site = np.asarray(tables.mutations.site)
    ages[np.unique(site)] = -np.inf
    np.maximum.at(ages, site, t)
    return ages


def edge_midpoint_ages(ts):
    """Dated age per site: midpoint of the edge above each mutation (oldest if several).

    Same result as `edge_midpoint_ages_slow` (kept as the reference), but vectorised over
    the edge table so a genome-wide draw (tens of millions of sites) takes seconds, not
    an hour.  Falls back to the slow loop when edge or site coordinates are not integers.
    """
    tables = ts.tables
    n = tables.sites.num_rows
    ages = np.full(n, np.nan)
    if tables.mutations.num_rows == 0:
        return ages
    left = np.asarray(tables.edges.left)
    spos = np.asarray(tables.sites.position)
    if (not np.all(left == np.floor(left)) or not np.all(spos == np.floor(spos))
            or ts.sequence_length >= 2 ** 32):
        return edge_midpoint_ages_slow(ts)
    child = np.asarray(tables.edges.child).astype(np.int64)
    key = (child << 32) | left.astype(np.int64)
    order = np.argsort(key, kind="stable")
    key_sorted = key[order]
    right = np.asarray(tables.edges.right)[order]
    parent = np.asarray(tables.edges.parent)[order]
    times = np.asarray(tables.nodes.time)
    mnode = np.asarray(tables.mutations.node).astype(np.int64)
    msite = np.asarray(tables.mutations.site)
    mpos = spos[msite].astype(np.int64)
    q = (mnode << 32) | mpos
    idx = np.searchsorted(key_sorted, q, side="right") - 1
    idc = np.clip(idx, 0, key_sorted.size - 1) if key_sorted.size else idx
    if key_sorted.size:
        has = (idx >= 0) & ((key_sorted[idc] >> 32) == mnode) & (right[idc] > mpos)
        par = np.where(has, parent[idc], -1)
    else:
        par = np.full(mnode.size, -1)
    age = np.where(par == -1, np.inf, 0.5 * (times[mnode] + times[np.maximum(par, 0)]))
    # oldest mutation per site; NaN where absent
    np.fmax.at(ages, msite, age)
    return ages


def _site_derived_state_counts(ts):
    """Number of distinct non-ancestral derived states per site (0 if no mutations)."""
    import tskit
    sites, muts = ts.tables.sites, ts.tables.mutations
    out = np.zeros(sites.num_rows, dtype=np.int32)
    if muts.num_rows == 0:
        return out
    anc = np.array(tskit.unpack_strings(sites.ancestral_state, sites.ancestral_state_offset), dtype=object)
    mder = np.array(tskit.unpack_strings(muts.derived_state, muts.derived_state_offset), dtype=object)
    uniq, code = np.unique(np.concatenate([mder, anc]), return_inverse=True)
    mcode, acode = code[: mder.size], code[mder.size:]
    msite = np.asarray(muts.site)
    keep = mcode != acode[msite]
    pairs = np.unique(np.stack([msite[keep], mcode[keep]], axis=1), axis=0)
    np.add.at(out, pairs[:, 0], 1)
    return out


def panel_derived_matrix(ts):
    """(n_samples, n_sites) 0/1 derived indicator for the ARG's samples.

    Sites with more than two alleles are marked -1 for every sample (callers should
    drop them); the derived allele is any non-ancestral state.  Equivalent to
    `panel_derived_matrix_slow` (the reference) but uses tskit's C genotype matrix.
    """
    out = (ts.genotype_matrix() != 0).T.astype(np.int8)       # allele 0 is the ancestral state
    multi = _site_derived_state_counts(ts) > 1
    out[:, multi] = -1
    return np.ascontiguousarray(out)
