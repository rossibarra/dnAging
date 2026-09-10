from argparse import Namespace
from math import comb

import mpmath as mp
import numpy as np

from betabinom_real_data import VariableNePhi, orientation
from precompute_betabinom_variable_ne import build
from tests._reference import moments


def test_variable_ne_table_matches_independent_moment_reference(tmp_path):
    ne_file = tmp_path / "ne.tsv"
    ne_file.write_text(
        "series\ttime_left\ttime_right\teffective_population_size\n"
        "posterior_mean\t0\t1000000\t50000\n"
    )
    out = tmp_path / "table"
    build(Namespace(ne=ne_file, ne_series="posterior_mean", n_panel=4,
                    u_min=1e-5, u_max=1e-2, n_u=4, dps=120, output=out))
    table = VariableNePhi(out)
    u = float(table.u[2])
    sample_t = 100.0
    origin_t = sample_t + 2 * 50000 * u
    got = float(table.integrate(2, 3, np.array([sample_t]),
                                np.array([origin_t]), np.array([1.0])))

    mom = moments(4, u, 1 / 100000)
    den = mp.fsum(comb(1, j) * (-1) ** j * mom[2 + j] for j in range(2))
    num = mp.fsum(comb(1, j) * (-1) ** j * mom[3 + j] for j in range(2))
    np.testing.assert_allclose(got, float(num / den), rtol=2e-6)


class _Mutation:
    def __init__(self, state):
        self.derived_state = state


class _Site:
    def __init__(self, state):
        self.ancestral_state = state


def test_orientation_handles_swap_and_strand_complement():
    assert orientation(_Site("A"), _Mutation("G"), "A", "G") == (True, "exact")
    assert orientation(_Site("G"), _Mutation("A"), "A", "G") == (False, "exact")
    assert orientation(_Site("T"), _Mutation("C"), "A", "G") == (
        True, "strand_complement")
