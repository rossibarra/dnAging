import numpy as np

from conftest import (FakeChunk, FakeStore, MISS, install_repo_stub, make_args,
                      make_table, pack)


def _run_one_site(monkeypatch, tmp_path, draw_ids, polarity):
    tab = make_table(n=6)
    tab["age_tau"] = np.array([0.1, 4.0])
    tab["n_panel"] = np.array([6])

    ancient = FakeChunk(["1"], [100], ["A"], ["C"], pack([[1]], [[1]]))
    panel = FakeChunk(["1"], [100], ["A"], ["C"],
                      pack([[1], [1], [0], [0], [0], [0]],
                           [[1], [1], [1], [1], [1], [1]]))
    mapping = {
        "anc.vcf": [(["ancient"], ancient)],
        "panel.vcf": [([f"p{i}" for i in range(6)], panel)],
    }
    store = FakeStore(2, {0: ([200.0] * len(draw_ids),
                              [300.0] * len(draw_ids), draw_ids)})
    install_repo_stub(monkeypatch, vcf_mapping=mapping, store=store,
                      polarity=np.asarray([polarity], dtype=np.int16), rows={100: 0})
    args = make_args(tmp_path, epsilon=0.0, extra=("--min-n", "6"))

    import posterior_sample_age_infer as inf
    return inf.run_chromosome(args, tab)


def test_site_with_missing_arg_draw_is_dropped(monkeypatch, tmp_path):
    _order, _grid, ll, stats = _run_one_site(
        monkeypatch, tmp_path, draw_ids=[0], polarity=[0, 0])

    assert np.array_equal(ll, np.zeros_like(ll))
    assert stats["sites_used"] == 0
    assert stats["sites_incomplete_draws"] == 1
    assert stats["draws_missing_interval"] == 1


def test_site_with_rejected_arg_draw_is_dropped(monkeypatch, tmp_path):
    _order, _grid, ll, stats = _run_one_site(
        monkeypatch, tmp_path, draw_ids=[0, 1], polarity=[0, MISS])

    assert np.array_equal(ll, np.zeros_like(ll))
    assert stats["sites_used"] == 0
    assert stats["sites_incomplete_draws"] == 1
    assert stats["sites_bad_polarity"] == 1
    assert stats["draws_bad_polarity"] == 1


def test_site_with_all_arg_draws_is_used(monkeypatch, tmp_path):
    _order, _grid, ll, stats = _run_one_site(
        monkeypatch, tmp_path, draw_ids=[0, 1], polarity=[0, 0])

    assert np.any(ll != 0.0)
    assert stats["sites_used"] == 1
    assert stats["sites_incomplete_draws"] == 0
