import gzip
import numpy as np
import tskit

from insertion_likelihood import (
    PiecewiseConstantNe,
    derived_probability_exact_time,
    derived_probability_uniform_edge,
    derived_probability_uniform_edge_grid,
)
from msprime_insertion_validation import vcf_call_matrix


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


def test_piecewise_ne_exact_time_matches_closed_form():
    tree, focal = two_tip_tree()
    demography = PiecewiseConstantNe(
        np.array([0.0, 2.0]), np.array([2.0, 20.0]),
        np.array([100.0, 50.0]))
    observed = derived_probability_exact_time(
        tree, focal, sample_time=0, mutation_time=5, ne=demography)
    total_hazard = 2 * 2 / (2 * 100) + 2 * 3 / (2 * 50)
    assert np.isclose(observed, 0.5 * (1 - np.exp(-total_hazard)))


def test_tau_inverse_crosses_epochs():
    demography = PiecewiseConstantNe(
        np.array([0.0, 100.0]), np.array([100.0, 1000.0]),
        np.array([10.0, 100.0]))
    assert np.isclose(demography.time_at_tau(5.5), 200.0)
    assert np.isclose(demography.tau(200.0), 5.5)


def test_multi_ancient_vcf_reader(tmp_path):
    path = tmp_path / "calls.vcf.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmodern_01\tmodern_02\tancient_01\tancient_02\n")
        handle.write("1\t1\t7\tA\tG\t.\tPASS\t.\tGT\t0\t1\t1\t0\n")
    names, calls = vcf_call_matrix(path, 2)
    assert names == ["ancient_01", "ancient_02"]
    assert calls[7][0] == 1
    assert np.array_equal(calls[7][1], [1, 0])


def test_explicit_final_epoch_extrapolation(tmp_path):
    path = tmp_path / "ne.tsv"
    path.write_text("series\ttime_left\ttime_right\teffective_population_size\n"
                    "x\t0\t100\t50\n")
    finite = PiecewiseConstantNe.from_tsv(path, "x")
    extended = PiecewiseConstantNe.from_tsv(path, "x", extend_last=True)
    assert np.isclose(finite.at([99])[0], 50)
    assert np.isclose(extended.at([1000000])[0], 50)
