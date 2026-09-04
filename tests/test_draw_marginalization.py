import numpy as np
import pytest

import posterior_sample_age_infer as inf
from conftest import (FakeChunk, FakeStore, install_repo_stub, make_args,
                      make_table, pack)


def _run_two_sites(monkeypatch, tmp_path, *, second_site_draw_values,
                   ploidy=1):
    tab = make_table(n=6)
    tab["age_tau"] = np.array([0.1, 4.0])
    tab["n_panel"] = np.array([6])

    dosage = 1 if ploidy == 1 else 2
    called = 1 if ploidy == 1 else 2
    ancient = FakeChunk(["1", "1"], [100, 200], ["A", "A"], ["C", "C"],
                        pack([[dosage, dosage]], [[called, called]]))
    panel = FakeChunk(
        ["1", "1"], [100, 200], ["A", "A"], ["C", "C"],
        pack([[1, 1], [1, 1], [0, 0], [0, 0], [0, 0], [0, 0]],
             [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1]]))
    mapping = {
        "anc.vcf": [(["ancient"], ancient)],
        "panel.vcf": [([f"p{i}" for i in range(6)], panel)],
    }
    store = FakeStore(2, {
        0: ([201.0, 202.0], [211.0, 212.0], [0, 1]),
        1: ([301.0, 302.0], [311.0, 312.0], [0, 1]),
    })
    install_repo_stub(
        monkeypatch, vcf_mapping=mapping, store=store,
        polarity=np.asarray([[0, 0], [0, 0]], dtype=np.int16),
        rows={100: 0, 200: 1})

    first = {201.0: 0.9, 202.0: 0.1}
    second = dict(zip((301.0, 302.0), second_site_draw_values))

    def fake_phi_lookup(_tab, _d0, t_lo, _t_hi, key="table", **_kwargs):
        value = (first | second)[float(t_lo)]
        if key == "table2":
            value = value ** 2
        return np.full(len(tab["Tgrid"]), value)

    monkeypatch.setattr(inf, "phi_lookup", fake_phi_lookup)
    args = make_args(tmp_path, ploidy=ploidy, epsilon=0.0,
                     extra=("--min-n", "6"))
    return inf.run_chromosome(args, tab)


@pytest.mark.parametrize(
    "second_site_draw_values, expected",
    [((0.9, 0.1), (0.9 ** 2 + 0.1 ** 2) / 2),
     ((0.1, 0.9), 0.9 * 0.1)],
)
def test_haploid_sites_are_multiplied_within_draw_before_marginalizing(
        monkeypatch, tmp_path, second_site_draw_values, expected):
    _order, _grid, ll, stats = _run_two_sites(
        monkeypatch, tmp_path,
        second_site_draw_values=second_site_draw_values)

    assert np.exp(ll) == pytest.approx(expected)
    assert stats["sites_used"] == 2
    # Both cases have per-site draw means of 0.5; the former implementation
    # therefore returned 0.25 for each and could not distinguish their pairing.
    if second_site_draw_values == (0.9, 0.1):
        assert np.exp(ll) != pytest.approx(0.25)


def test_diploid_likelihood_also_preserves_draw_identity(monkeypatch, tmp_path):
    _order, _grid, ll, stats = _run_two_sites(
        monkeypatch, tmp_path, second_site_draw_values=(0.9, 0.1), ploidy=2)

    expected = (0.9 ** 4 + 0.1 ** 4) / 2
    assert np.exp(ll) == pytest.approx(expected)
    assert stats["sites_used"] == 2
