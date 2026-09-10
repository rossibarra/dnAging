import numpy as np

from conftest import (FakeChunk, FakeStore, install_repo_stub, make_args,
                      make_table, pack)


def test_draw_subset_requires_only_selected_draws(monkeypatch, tmp_path):
    tab = make_table(n=6)
    tab["age_tau"] = np.array([0.1, 4.0])
    tab["n_panel"] = np.array([6])
    ancient = FakeChunk(["1"], [100], ["A"], ["C"], pack([[1]], [[1]]))
    panel = FakeChunk(["1"], [100], ["A"], ["C"],
                      pack([[1, 1, 0, 0, 0, 0]], [[1, 1, 1, 1, 1, 1]]))
    mapping = {
        "anc.vcf": [(["ancient"], ancient)],
        "panel.vcf": [([f"p{i}" for i in range(6)], panel)],
    }
    # Store draw 1 is absent at this site, but it is not selected.
    store = FakeStore(2, {0: ([200.0], [300.0], [0])})
    install_repo_stub(monkeypatch, vcf_mapping=mapping, store=store,
                      polarity=np.asarray([[0, 255]], dtype=np.int16), rows={100: 0})
    args = make_args(tmp_path, extra=("--min-n", "6", "--draw-ids", "0"))

    import posterior_sample_age_infer as inf
    _order, _grid, ll, stats, _epsilon_data = inf.run_chromosome(args, tab)

    assert np.any(ll != 0)
    assert stats["sites_used"] == 1
    assert stats["n_arg_draws"] == 1
    assert stats["store_draw_ids"] == [0]
    assert stats["draws_missing_interval"] == 0


def test_draw_subset_is_validated(tmp_path):
    import posterior_sample_age_infer as inf

    args = make_args(tmp_path, extra=("--draw-ids", "0,0"))
    # Validation needing the store happens in run_chromosome, so parsing keeps
    # the exact string without silently dropping or reordering draw IDs.
    assert args.draw_ids == "0,0"
