import json

import numpy as np

import posterior_sample_age_infer as inf
from conftest import (FakeChunk, FakeStore, install_repo_stub, make_args,
                      make_table, pack)


def test_save_epsilon_data_preserves_calls_and_phi(monkeypatch, tmp_path):
    tab = make_table(n=6)
    tab["age_tau"] = np.array([0.1, 4.0])
    tab["n_panel"] = np.array([6])
    ancient = FakeChunk(["1"], [100], ["A"], ["C"],
                        pack([[1, 0]], [[1, 1]]))
    panel = FakeChunk(["1"], [100], ["A"], ["C"],
                      pack([[1, 1, 0, 0, 0, 0]], [[1, 1, 1, 1, 1, 1]]))
    mapping = {
        "anc.vcf": [(["alt_sample", "ref_sample"], ancient)],
        "panel.vcf": [([f"p{i}" for i in range(6)], panel)],
    }
    store = FakeStore(1, {0: ([200.0], [300.0], [0])})
    install_repo_stub(monkeypatch, vcf_mapping=mapping, store=store,
                      polarity=np.asarray([[0]], dtype=np.int16), rows={100: 0})
    args = make_args(tmp_path, epsilon=0.01,
                     extra=("--min-n", "6", "--save-epsilon-data"))

    order, grid, _ll, stats, data = inf.run_chromosome(args, tab)
    assert order == ["alt_sample", "ref_sample"]
    assert stats["sites_used"] == 1
    assert data["position"].tolist() == [100]
    assert data["phi_alt"].shape == (1, 1, len(grid))
    assert data["observed_alt"].tolist() == [[1, 0]]
    assert data["called"].tolist() == [[1, 1]]
    assert data["panel_alt_count"].tolist() == [2]
    assert data["panel_called"].tolist() == [6]

    inf.write_outputs(tmp_path / "out", order, grid, np.zeros_like(grid),
                      np.zeros((2, len(grid))), stats, args, data)
    saved = np.load(tmp_path / "out" / "epsilon_calibration_data.npz",
                    allow_pickle=False)
    assert saved["phi_alt"].shape == (1, 1, len(grid))
    assert saved["samples"].tolist() == order
    assert json.loads(str(saved["meta"]))["format_version"] == 1
