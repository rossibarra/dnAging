#!/usr/bin/env python3
"""Direct ancient-lineage insertion probabilities on a fixed modern tree.

This prototype conditions on a local modern genealogy rather than replacing it
with a population-frequency diffusion. Backward from ancient sampling time T,
an added lineage coalesces with each extant modern lineage at rate 1/(2 Ne) in a
diploid population. It carries the derived allele exactly when its first
coalescence before mutation time is into ancestry below the mutation edge.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tskit


@dataclass(frozen=True)
class PiecewiseConstantNe:
    """Diploid effective size on contiguous half-open generation intervals."""

    left: np.ndarray
    right: np.ndarray
    size: np.ndarray

    @classmethod
    def from_tsv(cls, path, series=None, extend_last=False):
        with Path(path).open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if series is not None:
            rows = [r for r in rows if r.get("series") == series]
        if not rows:
            raise ValueError(f"no demographic epochs in {path}")
        right = np.asarray([float(r["time_right"]) for r in rows])
        if extend_last:
            right[-1] = np.inf
        obj = cls(np.asarray([float(r["time_left"]) for r in rows]), right,
                  np.asarray([float(r["effective_population_size"]) for r in rows]))
        obj.validate()
        return obj

    def validate(self):
        if not (len(self.left) == len(self.right) == len(self.size)):
            raise ValueError("demographic epoch arrays differ in length")
        if self.left[0] != 0 or np.any(self.right <= self.left):
            raise ValueError("demographic epochs must start at zero and have positive width")
        if np.any(~np.isfinite(self.left)) or np.any(~np.isfinite(self.right[:-1])):
            raise ValueError("only the final demographic boundary may be infinite")
        if np.any(~np.isfinite(self.size)) or np.any(self.size <= 0):
            raise ValueError("effective population sizes must be finite and positive")
        if len(self.left) > 1 and not np.allclose(self.left[1:], self.right[:-1]):
            raise ValueError("demographic epochs must be contiguous")

    def at(self, times):
        times = np.asarray(times, dtype=float)
        index = np.searchsorted(self.right, times, side="right")
        if np.any(index >= len(self.size)):
            raise ValueError("demography does not cover all requested times")
        return self.size[index]

    def internal_breaks(self, lower, upper):
        return self.right[(self.right > lower) & (self.right < upper)]

    def tau(self, time):
        total = 0.0
        for left, right, size in zip(self.left, self.right, self.size):
            width = max(0.0, min(float(time), right) - left)
            total += width / (2.0 * size)
            if time <= right:
                return total
        raise ValueError("demography does not cover requested time")

    def time_at_tau(self, target):
        remaining = float(target)
        for left, right, size in zip(self.left, self.right, self.size):
            capacity = (right - left) / (2.0 * size)
            if remaining <= capacity:
                return float(left + remaining * 2.0 * size)
            remaining -= capacity
        raise ValueError("demography does not extend to requested diffusion time")


def _as_demography(ne, maximum):
    if isinstance(ne, PiecewiseConstantNe):
        ne.validate()
        ne.at([maximum])
        return ne
    value = float(ne)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("Ne must be finite and positive")
    return PiecewiseConstantNe(np.array([0.0]), np.array([max(maximum, 1.0) + 1.0]),
                               np.array([value]))


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
    if not np.isfinite([T, m]).all() or T < 0 or m < 0:
        raise ValueError("sample time, mutation time and Ne must be finite and valid")
    if T >= m:
        return 0.0
    if focal_node == tskit.NULL or focal_node < 0:
        raise ValueError("focal_node must identify a node in the tree")

    intervals = _lineage_intervals(tree, focal_node)
    demography = _as_demography(ne, m)
    knots = {T, m}
    knots.update(demography.internal_breaks(T, m))
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
        local_ne = float(demography.at([midpoint])[0])
        interval_hazard = k * (upper - lower) / (2.0 * local_ne)
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


def _cumulative_insertion_functions(
    tree: tskit.Tree, focal_node: int, times, ne: float
):
    """Evaluate cumulative total hazard K, focal incidence B and integral C.

    K(x) = integral_0^x k(u)/(2Ne) du,
    B(x) = integral_0^x exp(-K(u)) d(u)/(2Ne) du, and
    C(x) = integral_0^x B(u) du.
    """
    times = np.asarray(times, dtype=float)
    if np.any(times < 0) or not np.all(np.isfinite(times)):
        raise ValueError("evaluation times must be finite and nonnegative")
    maximum = float(np.max(times)) if times.size else 0.0
    demography = _as_demography(ne, maximum)
    intervals = _lineage_intervals(tree, focal_node)
    knots = {0.0, maximum}
    for lower, upper, _ in intervals:
        if 0 < lower < maximum:
            knots.add(lower)
        if 0 < upper < maximum:
            knots.add(upper)
    knots.update(demography.internal_breaks(0, maximum))
    knots = np.asarray(sorted(knots), dtype=float)
    if len(knots) == 1:
        zeros = np.zeros_like(times)
        return zeros, zeros, zeros

    starts = knots[:-1]
    stops = knots[1:]
    k_values, d_values = [], []
    for lower, upper in zip(starts, stops):
        midpoint = 0.5 * (lower + upper)
        active = [(a, b, derived) for a, b, derived in intervals
                  if a <= midpoint < b]
        if not active:
            raise ValueError(f"modern tree has no lineage at time {midpoint:g}")
        k_values.append(len(active))
        d_values.append(sum(derived for _, _, derived in active))
    k_values = np.asarray(k_values, dtype=float)
    d_values = np.asarray(d_values, dtype=float)
    ne_values = demography.at(0.5 * (starts + stops))
    rates = k_values / (2.0 * ne_values)

    K_start = np.zeros(len(starts))
    B_start = np.zeros(len(starts))
    C_start = np.zeros(len(starts))
    for j in range(1, len(starts)):
        delta = stops[j - 1] - starts[j - 1]
        r = rates[j - 1]
        one_minus_exp = -np.expm1(-r * delta)
        amplitude = np.exp(-K_start[j - 1]) * d_values[j - 1] / k_values[j - 1]
        K_start[j] = K_start[j - 1] + r * delta
        B_start[j] = B_start[j - 1] + amplitude * one_minus_exp
        C_start[j] = (C_start[j - 1] + B_start[j - 1] * delta
                      + amplitude * (delta - one_minus_exp / r))

    index = np.searchsorted(knots, times, side="right") - 1
    index = np.clip(index, 0, len(starts) - 1)
    delta = times - starts[index]
    r = rates[index]
    one_minus_exp = -np.expm1(-r * delta)
    amplitude = np.exp(-K_start[index]) * d_values[index] / k_values[index]
    K = K_start[index] + r * delta
    B = B_start[index] + amplitude * one_minus_exp
    C = C_start[index] + B_start[index] * delta + amplitude * (
        delta - one_minus_exp / r
    )
    return K, B, C


def derived_probability_uniform_edge_grid(
    tree: tskit.Tree,
    focal_node: int,
    sample_times,
    edge_lower: float,
    edge_upper: float,
    ne: float,
):
    """Vectorized uniform-edge insertion probability over sample ages.

    This evaluates the nested coalescence/mutation-time integrals analytically
    between modern-tree events, avoiding per-site numerical quadrature.
    """
    sample_times = np.asarray(sample_times, dtype=float)
    lower = float(edge_lower)
    upper = float(edge_upper)
    if not upper > lower:
        raise ValueError("edge_upper must exceed edge_lower")
    sample_times_bounded = np.minimum(sample_times, upper)
    integration_lower = np.minimum(np.maximum(sample_times, lower), upper)
    evaluation_times = np.concatenate(
        (sample_times_bounded, integration_lower, [upper])
    )
    K, B, C = _cumulative_insertion_functions(
        tree, focal_node, evaluation_times, ne
    )
    n = len(sample_times)
    K_T = K[:n]
    B_T = B[:n]
    C_lower = C[n:2 * n]
    C_upper = C[-1]
    available = np.maximum(upper - integration_lower, 0.0)
    numerator = C_upper - C_lower - available * B_T
    probability = np.exp(K_T) * numerator / (upper - lower)
    probability[available == 0] = 0.0
    return np.clip(probability, 0.0, 1.0)
