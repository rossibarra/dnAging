import numpy as np
import tskit

from true_topology_polegon import polegon_tables, replace_node_times, root_stems


def example_tree():
    tables = tskit.TableCollection(10)
    a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    p = tables.nodes.add_row(time=10)
    tables.edges.add_row(0, 10, parent=p, child=a)
    tables.edges.add_row(0, 10, parent=p, child=b)
    site = tables.sites.add_row(5, ancestral_state="A")
    tables.mutations.add_row(site=site, node=a, derived_state="G")
    tables.sort()
    return tables.tree_sequence()


def test_polegon_tables_and_time_replacement():
    ts = example_tree()
    nodes, branches, mutations = polegon_tables(ts, 1e-8)
    assert nodes.shape == (3,)
    assert branches.shape == (3, 4)
    assert mutations.shape == (1, 4)
    dated = replace_node_times(ts, np.array([0, 0, 20]))
    assert dated.node(2).time == 20
    assert np.array_equal(dated.tables.edges.parent, ts.tables.edges.parent)
    assert np.all(np.isfinite(dated.tables.mutations.time))
    mutation = dated.mutation(0)
    assert dated.node(mutation.node).time <= mutation.time <= dated.node(2).time


def test_time_replacement_resorts_edges_after_parent_order_changes():
    tables = tskit.TableCollection(10)
    a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    p1 = tables.nodes.add_row(time=10)
    p2 = tables.nodes.add_row(time=20)
    tables.edges.add_row(0, 5, parent=p1, child=a)
    tables.edges.add_row(5, 10, parent=p2, child=b)
    tables.sort()
    ts = tables.tree_sequence()
    dated = replace_node_times(ts, np.array([0, 0, 30, 15]))
    assert dated.num_edges == 2
    assert np.array_equal(np.sort(dated.tables.edges.parent), np.array([2, 3]))


def test_root_stems_cover_sequence():
    stems = root_stems(example_tree())
    assert np.array_equal(stems[0], np.array([0, 10, -1, 2]))
