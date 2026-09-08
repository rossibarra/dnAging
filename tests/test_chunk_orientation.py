# -*- coding: utf-8 -*-
"""Genotype-code orientation must come from the counts, not from a shape guess.

normalizeTE's VcfChunk documents `codes` as (records, samples). The adapter used
to decide with `codes.shape[0] == len(names)`, which is true for ANY chunk holding
exactly as many records as there are samples -- a final partial chunk of 26
records in a 26-sample panel, for instance. Such a chunk was transposed silently:
not an error, just wrong allele counts, which can flip which site looks
monomorphic and invert an ancient observation.
"""
import numpy as np
import pytest

import posterior_sample_age_infer as inf
from conftest import FakeChunk, pack


def _chunk(n_site, n_samp, alt):
    return FakeChunk(["1"] * n_site, list(range(100, 100 + n_site)),
                     ["A"] * n_site, ["C"] * n_site,
                     pack(alt, np.ones((n_site, n_samp), np.int64)))


def test_square_chunk_is_not_transposed():
    """3 sites x 3 samples: the ambiguous case the old shape test got wrong."""
    # site 0 carried by every sample, sites 1 and 2 by none
    alt = np.array([[1, 1, 1], [0, 0, 0], [0, 0, 0]])
    names = ["p0", "p1", "p2"]
    codes = inf._codes_by_site(_chunk(3, 3, alt), names)
    assert codes.shape == (3, 3)
    counts = (codes >> 4).sum(axis=1)
    assert counts.tolist() == [3, 0, 0]
    # the transposed reading would give [1, 1, 1] -- three singletons instead of
    # one fixed site and two monomorphic ones
    assert counts.tolist() != [1, 1, 1]


def test_non_square_chunk_in_reader_layout():
    alt = np.array([[1, 1, 0, 0], [1, 0, 0, 0]])          # 2 sites, 4 samples
    codes = inf._codes_by_site(_chunk(2, 4, alt), ["p0", "p1", "p2", "p3"])
    assert codes.shape == (2, 4)
    assert (codes >> 4).sum(axis=1).tolist() == [2, 1]


def test_transposed_input_is_accepted_when_unambiguous():
    """A (samples, sites) chunk is still handled where the shape settles it."""
    # 2 sites, 4 samples, supplied as (samples, sites): only one reading fits
    ch = FakeChunk(["1", "1"], [100, 200], ["A", "A"], ["C", "C"],
                   pack(np.array([[1, 1], [1, 0], [0, 0], [0, 0]]),      # 4x2
                        np.ones((4, 2), np.int64)))
    codes = inf._codes_by_site(ch, ["p0", "p1", "p2", "p3"])
    assert codes.shape == (2, 4)
    assert (codes >> 4).sum(axis=1).tolist() == [2, 1]


def test_shape_matching_neither_layout_is_rejected():
    ch = FakeChunk(["1", "1"], [100, 200], ["A", "A"], ["C", "C"],
                   pack(np.zeros((5, 7), np.int64), np.ones((5, 7), np.int64)))
    with pytest.raises(SystemExit, match="neither"):
        inf._codes_by_site(ch, ["p0", "p1", "p2"])


def test_chunk_codes_selects_requested_samples_from_a_square_chunk():
    """_chunk_codes returns (len(want), sites) and must pick the right samples."""
    alt = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])      # 3 sites, 3 samples
    ch = _chunk(3, 3, alt)
    got = inf._chunk_codes(ch, ["a", "b", "c"], ["c"])
    assert got.shape == (1, 3)
    assert (got >> 4).tolist() == [[0, 0, 1]]              # sample c carries site 2
