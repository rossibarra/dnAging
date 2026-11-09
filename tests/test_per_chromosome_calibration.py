import numpy as np

from per_chromosome_calibration import combined_age_loglik


def test_combined_age_loglik_marginalizes_scales_then_sums_chromosomes():
    ll = np.zeros((2, 1, 2, 3))
    ll[0, 0, 0] = [0, 1, 2]
    ll[1, 0, 0] = [2, 1, 0]
    ll[0, 0, 1] = [1, 0, 1]
    ll[1, 0, 1] = [0, 1, 0]
    weights = np.full((2, 2, 1), -np.log(2))
    observed = combined_age_loglik(ll, weights)
    expected = (
        np.logaddexp(ll[0, 0, 0], ll[1, 0, 0]) - np.log(2)
        + np.logaddexp(ll[0, 0, 1], ll[1, 0, 1]) - np.log(2)
    )
    assert np.allclose(observed, expected)
