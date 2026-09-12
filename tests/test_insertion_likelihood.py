import numpy as np
import tskit

from insertion_likelihood import (
    derived_probability_exact_time,
    derived_probability_uniform_edge,
    derived_probability_uniform_edge_grid,
)


def two_tip_tree(root_time=10.0):
    tables = tskit.TableCollection(sequence_length=1)
    left = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    right = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    root = tables.nodes.add_row(time=root_time)
    tables.edges.add_row(0, 1, root, left)
    tables.edges.add_row(0, 1, root, right)
    tables.sort()
    return tables.tree_sequence().first(), left


def test_exact_time_two_tip_tree_matches_closed_form():
    tree, focal = two_tip_tree()
    ne = 100.0
    observed = derived_probability_exact_time(
        tree, focal, sample_time=0, mutation_time=5, ne=ne
    )
    expected = 0.5 * (1.0 - np.exp(-5.0 / ne))
    assert np.isclose(observed, expected)


def test_sample_older_than_mutation_is_ancestral():
    tree, focal = two_tip_tree()
    assert derived_probability_exact_time(
        tree, focal, sample_time=6, mutation_time=5, ne=100
    ) == 0.0


def test_uniform_edge_two_tip_tree_matches_closed_form():
    tree, focal = two_tip_tree()
    ne = 100.0
    width = 10.0
    observed = derived_probability_uniform_edge(
        tree, focal, sample_time=0, edge_lower=0, edge_upper=width, ne=ne
    )
    expected = 0.5 * (1.0 - (ne / width) * (1.0 - np.exp(-width / ne)))
    assert np.isclose(observed, expected, rtol=1e-12)


def test_uniform_edge_respects_mutation_existence_boundary():
    tree, focal = two_tip_tree()
    probability = derived_probability_uniform_edge(
        tree, focal, sample_time=8, edge_lower=0, edge_upper=10, ne=100
    )
    assert 0 < probability < 0.01


def test_vectorized_uniform_edge_matches_quadrature():
    tree, focal = two_tip_tree()
    sample_times = np.array([0.0, 1.5, 5.0, 8.0, 10.0, 12.0])
    observed = derived_probability_uniform_edge_grid(
        tree, focal, sample_times, edge_lower=0, edge_upper=10, ne=100
    )
    expected = np.array([
        derived_probability_uniform_edge(
            tree, focal, t, edge_lower=0, edge_upper=10, ne=100,
            quadrature_order=128,
        )
        for t in sample_times
    ])
    assert np.allclose(observed, expected, atol=2e-6)
