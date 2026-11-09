import numpy as np
import tskit

from chunked_true_topology_polegon import extract_chunk, merge_chunks


def example_ts():
    tables = tskit.TableCollection(sequence_length=10)
    samples = [tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE) for _ in range(2)]
    left_root = tables.nodes.add_row(time=2)
    right_root = tables.nodes.add_row(time=3)
    for left, right, root in ((0, 5, left_root), (5, 10, right_root)):
        for sample in samples:
            tables.edges.add_row(left, right, root, sample)
    for position, node in ((2, samples[0]), (7, samples[1])):
        site = tables.sites.add_row(position, ancestral_state="0")
        tables.mutations.add_row(site, node=node, derived_state="1")
    tables.sort()
    return tables.tree_sequence()


def test_extract_and_merge_chunks_duplicates_internal_nodes_only():
    truth = example_ts()
    first, _ = extract_chunk(truth, 0, 5)
    second, _ = extract_chunk(truth, 5, 10)
    merged = merge_chunks([(0, first), (5, second)], truth.sequence_length)
    assert merged.num_samples == truth.num_samples
    assert merged.num_nodes == truth.num_samples + 2
    assert np.array_equal(merged.sites_position, truth.sites_position)
    for position in (2, 7):
        assert set(merged.at(position).samples()) == set(truth.at(position).samples())


def test_merge_preserves_a_gap_outside_production_chunks():
    truth = example_ts()
    middle, _ = extract_chunk(truth, 2, 8)
    merged = merge_chunks([(2, middle)], truth.sequence_length)
    assert merged.sequence_length == truth.sequence_length
    assert merged.at(1).num_edges == 0
    assert merged.at(9).num_edges == 0
