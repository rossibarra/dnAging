import numpy as np

from radiocarbon_inferred_stage2 import resampled_likelihood


def test_resampled_likelihood_sums_blocks_then_marginalizes_draws_by_chromosome():
    # epsilon, draw, chromosome, block, sample, age
    blocks=np.zeros((2,2,2,2,1,1))
    blocks[:,0,1]=-100
    blocks[:,1,0]=-100
    weights=np.ones((2,2),dtype=int)
    observed=resampled_likelihood(blocks,weights)
    expected=2*np.log((1+np.exp(-200))/2)
    assert observed.shape==(2,1,1)
    np.testing.assert_allclose(observed,expected)
