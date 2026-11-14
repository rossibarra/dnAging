"""Pairwise coalescence rate c(t) for the non-sharing-ratio estimator.

Two sources, both returned as piecewise-constant rates on epoch boundaries
`edges` (len m+1, edges[0]=0, last edge inf) so they plug straight into
`nonsharing_ratio.invert_hazard`:

* ``rates_from_arg``  : from the ARG's own genome-averaged lineage counts
  kappa(t) (the default; lineages are lost at rate C(kappa,2) c(t)).
* ``rates_from_demes``: from a demes-format demography.  For a single panmictic
  deme c(t) = 1 / (ploidy * N(t)), with N(t) the deme's size in individuals.

Time units are generations throughout (the units of SINGER/tskit node times).
"""
from __future__ import annotations

import math
import warnings

import numpy as np

from nonsharing_ratio import lineage_moments, rate_from_lineage_counts


def make_edges(grid_max, grid_n):
    """Epoch boundaries 0..grid_max in `grid_n` equal steps, plus an open last epoch."""
    if grid_max <= 0 or grid_n < 1:
        raise ValueError("grid_max must be positive and grid_n >= 1")
    return np.concatenate([np.linspace(0.0, float(grid_max), int(grid_n) + 1), [np.inf]])


def _with_open_tail(edges_finite, rates_finite, tail=None):
    """Append the open final epoch; its rate defaults to the last positive finite rate."""
    rates_finite = np.asarray(rates_finite, dtype=float)
    if tail is None:
        positive = rates_finite[rates_finite > 0]
        if positive.size == 0:
            raise ValueError("coalescence rate is zero everywhere on the grid")
        tail = positive[-1]
    edges = np.concatenate([np.asarray(edges_finite, dtype=float), [np.inf]])
    return edges, np.concatenate([rates_finite, [tail]])


def rates_from_arg(ts, edges):
    """c(t) per epoch from the lineage counts of one tree sequence (one posterior draw).

    `edges` is the output of `make_edges`; the open last epoch inherits the last
    estimated rate.  Where the panel has coalesced to a single lineage the rate
    is not identifiable and is carried forward from the previous epoch.
    """
    finite = np.asarray(edges, dtype=float)
    finite = finite[np.isfinite(finite)]
    kappa, pairs = lineage_moments(ts, finite)
    c = rate_from_lineage_counts(finite, kappa, pairs)
    return _with_open_tail(finite, c)


# ---------------------------------------------------------------- demes ---

def _generations(graph, generation_time):
    """Factor converting the graph's time unit to generations."""
    unit = graph.time_units
    if unit == "generations":
        return 1.0
    gt = generation_time if generation_time is not None else graph.generation_time
    if gt is None:
        raise ValueError(f"demes time_units is {unit!r}; give a generation time "
                         "(graph.generation_time or generation_time=) to convert to generations")
    return 1.0 / float(gt)


def _deme_chain(graph, deme):
    """The named deme followed back through single-ancestor links, youngest first."""
    chain = [graph[deme]]
    while True:
        anc = chain[-1].ancestors
        if not anc:
            return chain
        if len(anc) > 1:
            raise ValueError(f"deme {chain[-1].name!r} has several ancestors; "
                             "a single panmictic lineage is required")
        chain.append(graph[anc[0]])


def _epoch_pieces(chain, to_gen):
    """(t_young, t_old, size_young, size_old, kind) per epoch in generations, young->old.

    Epochs of ancestor demes only count for times before the child existed.
    """
    pieces = []
    for d in chain:
        for ep in d.epochs:
            t_old, t_young = ep.start_time, ep.end_time
            if t_old <= t_young:
                continue
            pieces.append((t_young * to_gen, t_old * to_gen,
                           float(ep.end_size), float(ep.start_size), ep.size_function))
    pieces.sort(key=lambda p: p[0])
    return pieces


def _size_at(t, piece):
    """Deme size (individuals) at time t inside an epoch piece."""
    t_young, t_old, n_young, n_old, kind = piece
    if kind == "constant" or not math.isfinite(t_old):
        return n_old
    frac = (t_old - t) / (t_old - t_young)  # 0 at the old end, 1 at the young end
    if kind == "exponential":
        return n_old * (n_young / n_old) ** frac
    if kind == "linear":
        return n_old + (n_young - n_old) * frac
    raise ValueError(f"unsupported demes size_function {kind!r}")


_GL_X, _GL_W = np.polynomial.legendre.leggauss(16)


def _hazard_between(a, b, pieces, ploidy):
    """int_a^b dt / (ploidy N(t)) over the epoch pieces (exact for constant epochs)."""
    total = 0.0
    for piece in pieces:
        lo, hi = max(a, piece[0]), min(b, piece[1])
        if hi <= lo:
            continue
        if piece[4] == "constant" or not math.isfinite(piece[1]):
            total += (hi - lo) / (ploidy * piece[3])
            continue
        x = 0.5 * (hi - lo) * _GL_X + 0.5 * (hi + lo)
        inv = np.array([1.0 / _size_at(t, piece) for t in x])
        total += 0.5 * (hi - lo) * float(np.dot(_GL_W, inv)) / ploidy
    return total


def rates_from_demes(path, edges, deme=None, ploidy=2, generation_time=None):
    """c(t) per epoch from a demes YAML file (bin-averaged 1 / (ploidy N(t))).

    Only a single panmictic lineage is supported.  With several demes, name the
    one the samples belong to; migration and pulses are then ignored (a warning
    is emitted), so the result is only meaningful if they are negligible.
    """
    import demes

    graph = demes.load(str(path))
    if deme is None:
        if len(graph.demes) != 1:
            raise ValueError("demography has several demes; choose one with deme=")
        deme = graph.demes[0].name
    elif len(graph.demes) > 1 and (graph.migrations or graph.pulses):
        warnings.warn("demography has migrations/pulses, which are ignored when "
                      "deriving a single-deme c(t)")
    to_gen = _generations(graph, generation_time)
    pieces = _epoch_pieces(_deme_chain(graph, deme), to_gen)
    if not pieces or pieces[0][0] > 0:
        raise ValueError(f"deme {deme!r} does not exist at time 0")
    finite = np.asarray(edges, dtype=float)
    finite = finite[np.isfinite(finite)]
    c = np.empty(finite.size - 1)
    for i in range(c.size):
        c[i] = _hazard_between(finite[i], finite[i + 1], pieces, ploidy) / (finite[i + 1] - finite[i])
    # The open last epoch inherits the final bin's rate; only c on [0, T] matters.
    return _with_open_tail(finite, c)


def rates_from_ne_tsv(path, edges, series=None, ploidy=2):
    """c(t) per epoch from the repo's Ne trajectory TSV (bin-averaged 1 / (ploidy N(t))).

    The TSV is the one read by `insertion_likelihood.PiecewiseConstantNe.from_tsv`
    (columns time_left, time_right, effective_population_size and optionally series).
    Times are generations; N is the diploid effective size, so the pairwise
    coalescence rate is 1 / (2 N) with ploidy=2 (the same convention as the insertion
    code's tau = dt / (2 N)).  The last epoch is extended to infinity.  If the file has
    several series, pass `series` (a file with a `series` column and no choice is an
    error, as is a file with several replicate_index values within one series).
    """
    import csv
    from insertion_likelihood import PiecewiseConstantNe

    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if rows and "series" in rows[0]:
        names = sorted({r["series"] for r in rows})
        if series is None:
            if len(names) != 1:
                raise ValueError(f"{path} has several series {names}; choose one with series=")
            series = names[0]
        elif series not in names:
            raise ValueError(f"series {series!r} not in {path}: {names}")
        reps = {r.get("replicate_index") for r in rows if r["series"] == series}
        if len(reps) > 1:
            raise ValueError(f"series {series!r} in {path} has several replicates; "
                             "use a single-trajectory series such as posterior_mean")
    ne = PiecewiseConstantNe.from_tsv(path, series=series, extend_last=True)
    finite = np.asarray(edges, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite[-1] <= 0:
        raise ValueError("edges must extend beyond 0")
    c = np.empty(finite.size - 1)
    for i in range(c.size):
        a, b = finite[i], finite[i + 1]
        total = 0.0
        for left, right, size in zip(ne.left, ne.right, ne.size):
            lo, hi = max(a, left), min(b, right)
            if hi > lo:
                total += (hi - lo) / (ploidy * size)
        c[i] = total / (b - a)
    return _with_open_tail(finite, c)


def resolve_rates(source, edges, ts=None, deme=None, ploidy=2, generation_time=None):
    """Dispatch: source is 'arg' (needs `ts`) or a path to a demes file."""
    if source == "arg":
        if ts is None:
            raise ValueError("source 'arg' needs a tree sequence")
        return rates_from_arg(ts, edges)
    return rates_from_demes(source, edges, deme=deme, ploidy=ploidy,
                            generation_time=generation_time)
