#!/usr/bin/env python3
"""Direct ancient-lineage insertion probabilities on a fixed modern tree.

This prototype conditions on a local modern genealogy rather than replacing it
with a population-frequency diffusion. Backward from ancient sampling time T,
an added lineage coalesces with each extant modern lineage at rate 1/(2 Ne) in a
diploid population. It carries the derived allele exactly when its first
coalescence before mutation time is into ancestry below the mutation edge.
"""

from __future__ import annotations

import numpy as np
import tskit


def _lineage_intervals(tree: tskit.Tree, focal_node: int):
    """Return (lower, upper, is_focal_descendant) for every tree edge."""
    ts = tree.tree_sequence
    focal_descendants = set(tree.nodes(focal_node))
    intervals = []
    for parent in tree.nodes():
        for child in tree.children(parent):
            intervals.append((float(ts.node(child).time),
                              float(ts.node(parent).time),
                              child in focal_descendants))
    return intervals


def derived_probability_exact_time(
    tree: tskit.Tree,
    focal_node: int,
    sample_time: float,
    mutation_time: float,
    ne: float,
) -> float:
    """Probability an inserted ancient haplotype carries the focal mutation.

    The calculation is conditional on the fixed modern tree and assumes a
    diploid constant-size Kingman coalescent. Times are generations before the
    present. A sample older than the mutation is necessarily ancestral.
    """
    T = float(sample_time)
    m = float(mutation_time)
    ne = float(ne)
    if not np.isfinite([T, m, ne]).all() or T < 0 or m < 0 or ne <= 0:
        raise ValueError("sample time, mutation time and Ne must be finite and valid")
    if T >= m:
        return 0.0
    if focal_node == tskit.NULL or focal_node < 0:
        raise ValueError("focal_node must identify a node in the tree")

    intervals = _lineage_intervals(tree, focal_node)
    knots = {T, m}
    for lower, upper, _ in intervals:
        if T < lower < m:
            knots.add(lower)
        if T < upper < m:
            knots.add(upper)
    knots = sorted(knots)

    survival = 1.0
    probability = 0.0
    for lower, upper in zip(knots[:-1], knots[1:]):
        midpoint = 0.5 * (lower + upper)
        active = [(a, b, derived) for a, b, derived in intervals
                  if a <= midpoint < b]
        k = len(active)
        if k == 0:
            raise ValueError(
                f"modern tree has no ancestral lineage at time {midpoint:g}"
            )
        d = sum(derived for _, _, derived in active)
        interval_hazard = k * (upper - lower) / (2.0 * ne)
        coalescence_probability = -np.expm1(-interval_hazard)
        probability += survival * (d / k) * coalescence_probability
        survival *= np.exp(-interval_hazard)
    return float(np.clip(probability, 0.0, 1.0))


def derived_probability_uniform_edge(
    tree: tskit.Tree,
    focal_node: int,
    sample_time: float,
    edge_lower: float,
    edge_upper: float,
    ne: float,
    quadrature_order: int = 32,
) -> float:
    """Average the insertion probability over a uniform mutation position."""
    lower = float(edge_lower)
    upper = float(edge_upper)
    if not upper > lower:
        raise ValueError("edge_upper must exceed edge_lower")
    if quadrature_order < 2:
        raise ValueError("quadrature_order must be at least two")
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    mutation_times = 0.5 * (upper - lower) * nodes + 0.5 * (upper + lower)
    probabilities = np.asarray([
        derived_probability_exact_time(
            tree, focal_node, sample_time, mutation_time, ne
        )
        for mutation_time in mutation_times
    ])
    return float(0.5 * np.dot(weights, probabilities))
