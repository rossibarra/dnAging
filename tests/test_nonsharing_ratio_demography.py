import numpy as np
import pytest
import tskit

demes = pytest.importorskip("demes")

from nonsharing_ratio_demography import (
    make_edges, rates_from_arg, rates_from_demes, resolve_rates,
)
from nonsharing_ratio import cumulative_hazard


def write(tmp_path, graph, name="g.yaml"):
    p = tmp_path / name
    demes.dump(graph, str(p))
    return p


def constant_graph(n=20000, units="generations", gt=None):
    b = demes.Builder(time_units=units, generation_time=gt)
    b.add_deme("pop", epochs=[dict(start_size=n)])
    return b.resolve()


def test_constant_deme_matches_one_over_two_n(tmp_path):
    edges = make_edges(10000, 50)
    e, r = rates_from_demes(write(tmp_path, constant_graph()), edges)
    assert np.allclose(r, 1 / 40000)
    assert np.isinf(e[-1]) and r.size == e.size - 1


def test_ploidy_one(tmp_path):
    e, r = rates_from_demes(write(tmp_path, constant_graph()), make_edges(1000, 10), ploidy=1)
    assert np.allclose(r, 1 / 20000)


def test_bottleneck_epochs(tmp_path):
    b = demes.Builder(time_units="generations")
    b.add_deme("pop", epochs=[dict(start_size=25000, end_time=6000),
                              dict(start_size=5000, end_time=1500),
                              dict(start_size=50000, end_time=0)])
    p = write(tmp_path, b.resolve())
    e, r = rates_from_demes(p, make_edges(8000, 80))  # bins of 100 align to boundaries
    mid = lambda t: r[int(t // 100)]
    assert np.isclose(mid(500), 1 / 100000)
    assert np.isclose(mid(3000), 1 / 10000)
    assert np.isclose(mid(7000), 1 / 50000)


def test_misaligned_bin_averages_hazard(tmp_path):
    b = demes.Builder(time_units="generations")
    b.add_deme("pop", epochs=[dict(start_size=10000, end_time=150), dict(start_size=20000, end_time=0)])
    e, r = rates_from_demes(write(tmp_path, b.resolve()), make_edges(200, 1))
    # one bin [0,200): [0,150) has N=20000 (rate 1/40000), [150,200) has N=10000 (rate 1/20000)
    expected = (150 / 40000 + 50 / 20000) / 200
    assert np.isclose(r[0], expected)


def test_exponential_growth_hazard_matches_closed_form(tmp_path):
    # N grows from 1000 (older, t=1000) to 10000 (now).
    b = demes.Builder(time_units="generations")
    b.add_deme("pop", epochs=[dict(start_size=1000, end_time=1000),
                              dict(start_size=1000, end_size=10000, end_time=0,
                                   size_function="exponential")])
    e, r = rates_from_demes(write(tmp_path, b.resolve()), make_edges(1000, 4))
    k = np.log(10) / 1000  # N(t) = 10000 exp(-k t)
    exact = (np.exp(k * 1000) - 1) / (2 * 10000 * k)  # int_0^1000 dt/(2N)
    assert np.isclose(cumulative_hazard(e[:-1], r[:-1])[-1], exact, rtol=1e-6)


def test_years_are_converted_to_generations(tmp_path):
    # 25-year generations: a size change at 2500 years is at 100 generations.
    b = demes.Builder(time_units="years", generation_time=25)
    b.add_deme("pop", epochs=[dict(start_size=10000, end_time=2500), dict(start_size=20000, end_time=0)])
    e, r = rates_from_demes(write(tmp_path, b.resolve()), make_edges(200, 2))
    assert np.allclose(r[:-1], [1 / 40000, 1 / 20000])


def test_generation_time_override_on_years_graph(tmp_path):
    b = demes.Builder(time_units="years", generation_time=25)
    b.add_deme("pop", epochs=[dict(start_size=10000, end_time=2500), dict(start_size=20000, end_time=0)])
    e, r = rates_from_demes(write(tmp_path, b.resolve()), make_edges(200, 2), generation_time=50)
    # with 50-year generations the change is at 50 generations, inside the first bin
    assert r[0] == pytest.approx((50 / 40000 + 50 / 20000) / 100)


def test_multi_deme_requires_choice(tmp_path):
    b = demes.Builder(time_units="generations")
    b.add_deme("anc", epochs=[dict(start_size=1000, end_time=500)])
    b.add_deme("a", ancestors=["anc"], epochs=[dict(start_size=2000)])
    b.add_deme("b", ancestors=["anc"], epochs=[dict(start_size=3000)])
    p = write(tmp_path, b.resolve())
    with pytest.raises(ValueError):
        rates_from_demes(p, make_edges(1000, 10))
    e, r = rates_from_demes(p, make_edges(1000, 10), deme="a")  # follows ancestor beyond 500
    assert np.isclose(r[0], 1 / 4000) and np.isclose(r[9], 1 / 2000)


def two_tip_arg(root=10000.0, n=26):
    tables = tskit.TableCollection(sequence_length=100)
    leaves = [tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0) for _ in range(n)]
    # caterpillar with coalescences evenly spaced
    parent_time = np.linspace(500, root, n - 1)
    prev = leaves[0]
    for leaf, t in zip(leaves[1:], parent_time):
        node = tables.nodes.add_row(time=t)
        tables.edges.add_row(0, 100, node, prev)
        tables.edges.add_row(0, 100, node, leaf)
        prev = node
    tables.sort()
    return tables.tree_sequence()


def test_rates_from_arg_is_positive_and_shaped():
    ts = two_tip_arg()
    e, r = rates_from_arg(ts, make_edges(12000, 24))
    assert r.size == e.size - 1 and np.isinf(e[-1])
    assert np.all(r >= 0) and r[-1] > 0


def test_resolve_dispatch(tmp_path):
    ts = two_tip_arg()
    edges = make_edges(1000, 5)
    assert resolve_rates("arg", edges, ts=ts)[1].size == 6
    with pytest.raises(ValueError):
        resolve_rates("arg", edges)
    e, r = resolve_rates(write(tmp_path, constant_graph()), edges)
    assert np.allclose(r, 1 / 40000)
