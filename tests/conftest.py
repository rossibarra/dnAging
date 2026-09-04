# -*- coding: utf-8 -*-
"""Shared synthetic fixtures.

`normalize_tes` is not importable and there is no test data (no VCF, no SINGER
store, no Ne TSV, no .npz table), so everything is built here and the adapter
boundary -- `_import_repo` and `_resolve_rows` -- is stubbed. The real
`read_ancient` / `read_panel_alt` / `run_chromosome` code paths run unmodified
against these fakes; only the two functions that would import the repo are
replaced.
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import posterior_sample_age_infer as inf   # noqa: E402

N_PANEL = 6            # table's n_sample: panel is 6 called alleles
MISS = inf._MISS


# ---------------------------------------------------------------------------
# packed genotype codes: low nibble = alleles called, high nibble = ALT dosage
# ---------------------------------------------------------------------------


def pack(alt, called):
    """VCF chunk code as the adapter decodes it: (alt << 4) | called."""
    return (np.asarray(alt, np.int64) << 4) | np.asarray(called, np.int64)


class FakeChunk:
    """A VcfChunk-shaped object.

    `chrom` is POISONED on purpose: `_chunk_sites` must reach the attribute via
    `_attr(chunk, "chromosomes", "chrom")`, which only touches the fallback when
    the first name is missing. A nested `getattr(c, "chromosomes", getattr(c,
    "chrom"))` evaluates the fallback eagerly and blows up here -- so every test
    that reads a chunk also pins fix 2.
    """

    def __init__(self, chrom, positions, ref, alt, codes):
        self.chromosomes = np.asarray([str(c) for c in chrom])
        self.positions = np.asarray(positions, dtype=np.int64)
        self.ref = np.asarray(ref)
        self.alt = np.asarray(alt)
        self.codes = np.asarray(codes)

    @property
    def chrom(self):
        raise AssertionError("eager fallback: 'chromosomes' exists, do not read 'chrom'")


class FakeStore:
    """SINGER interval store: n_posterior_draws + intervals(rows) -> below/above/draw_id."""

    def __init__(self, n_draws, rows):
        self.n_posterior_draws = n_draws
        self._rows = rows            # row -> (below, above, draw_id)

    def intervals(self, rows):
        r = int(np.asarray(rows, dtype=np.int64)[0])
        below, above, draw = self._rows[r]
        return SimpleNamespace(below=np.asarray(below, float),
                               above=np.asarray(above, float),
                               draw_id=np.asarray(draw, np.int64))


def make_read_vcf_chunks(mapping):
    """A read_vcf_chunks stub: {path string: [(names, chunk), ...]}."""
    def read_vcf_chunks(path, sample_filter=None, chunk_records=None,
                        multiallelic=None, progress=None):
        assert multiallelic == "skip", "the caller must ask for multiallelic skipping"
        for names, chunk in mapping[str(path)]:
            yield list(names), chunk, None, None
    return read_vcf_chunks


def install_repo_stub(monkeypatch, *, vcf_mapping, store, polarity, rows):
    """Stub the whole adapter boundary; return nothing but the patched module."""
    read_vcf_chunks = make_read_vcf_chunks(vcf_mapping)
    monkeypatch.setattr(inf, "_import_repo", lambda: (
        lambda _p: store, lambda _p, _s: (polarity, None), read_vcf_chunks, MISS))
    monkeypatch.setattr(inf, "_resolve_rows",
                        lambda _store, _chrom, pos: np.asarray(
                            [rows.get(int(p), -1) for p in pos], dtype=np.int64))
    return inf


# ---------------------------------------------------------------------------
# synthetic frequency table
# ---------------------------------------------------------------------------


def make_table(n=N_PANEL, with_table2=True):
    """A table whose lookups are exactly predictable.

    Constant along the age axis and with every Tgrid point BELOW age[0], so
    phi_lookup's log-age interpolation and branch quadrature are both exact and
    the T >= t_i mask never fires: phi[d0-1, k] == 0.1*d0 + 0.02*k. That makes
    the inference tests assertions about the LIKELIHOOD, not about interpolation.
    """
    age = np.array([100.0, 1000.0])
    Tgrid = np.array([0.0, 10.0, 20.0, 30.0])
    base = 0.1 * np.arange(1, n + 1)[:, None] + 0.02 * np.arange(len(Tgrid))[None, :]
    table = np.repeat(base[:, None, :], len(age), axis=1).astype(np.float32)
    tab = {"table": table, "d0": np.arange(1, n + 1), "age": age,
           "Tgrid": Tgrid, "n_sample": n}
    if with_table2:
        # any value in [p^2, p] is admissible; 0.5(p^2+p) is strictly inside
        tab["table2"] = (0.5 * (table.astype(np.float64) ** 2
                                + table.astype(np.float64))).astype(np.float32)
    return tab


@pytest.fixture
def tab():
    return make_table()


@pytest.fixture
def tab_no2():
    return make_table(with_table2=False)


def expected_phi(tab, d0):
    """The value make_table() guarantees phi_lookup returns for this d0."""
    return 0.1 * d0 + 0.02 * np.arange(len(tab["Tgrid"]))


# ---------------------------------------------------------------------------
# args
# ---------------------------------------------------------------------------


def make_args(tmp_path, ploidy=1, epsilon=1e-3, chrom="1", extra=()):
    argv = ["--freq-table", str(tmp_path / "tab.npz"),
            "--store", str(tmp_path / "store"),
            "--draw-polarity", str(tmp_path / "pol"),
            "--panel-vcf", "panel.vcf", "--vcf", "anc.vcf",
            "--chrom", chrom, "--output", str(tmp_path / "out"),
            "--ploidy", str(ploidy), "--epsilon", str(epsilon), "--quiet"]
    return inf.parse_args(argv + list(extra))


# ---------------------------------------------------------------------------
# Ne TSV
# ---------------------------------------------------------------------------


def write_ne(path, rows, series_col=True, header=None):
    """rows: (series, time_left, time_right, Ne) tuples written as a TSV."""
    cols = (["series"] if series_col else []) + ["time_left", "time_right",
                                                 "effective_population_size"]
    lines = ["\t".join(header if header is not None else cols)]
    for s, lo, hi, ne in rows:
        vals = ([s] if series_col else []) + [str(lo), str(hi), str(ne)]
        lines.append("\t".join(vals))
    Path(path).write_text("\n".join(lines) + "\n")
    return path
