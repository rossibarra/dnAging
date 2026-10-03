import numpy as np
import tskit
from insertion_real_data import (block_bootstrap_maps, call_loglik, chromosome_layout,
                                 marginalize_independent_chromosomes)

def test_draws_are_marginalized_independently_by_chromosome():
    # Each chromosome has one good draw, but they are different draw indices.
    ll=np.array([[[[0.]], [[-100.]]], [[[-100.]], [[0.]]]])
    observed=marginalize_independent_chromosomes(ll)[0,0]
    expected=2*np.log((1+np.exp(-100))/2)
    assert np.isclose(observed,expected)
    # Incorrect genome-wide coupling would be about -100, not -log(4).
    assert observed > -2

def test_one_draw_reduces_to_chromosome_sum():
    ll=np.arange(6,dtype=float).reshape(1,2,1,3)
    assert np.array_equal(marginalize_independent_chromosomes(ll),ll[0].sum(0))

def test_argtest_chromosome_offsets_are_read_from_metadata():
    tables=tskit.TableCollection(sequence_length=30)
    tables.metadata_schema=tskit.MetadataSchema({"codec":"json"})
    tables.metadata={"chrom_offsets":[{"chrom":str(i),"offset":3*(i-1),"length":3}
                                                for i in range(1,11)]}
    layout=chromosome_layout(tables.tree_sequence())
    assert layout["1"] == (0.0,3.0)
    assert layout["10"] == (27.0,3.0)

def test_block_bootstrap_maps_preserves_shape_and_constant_blocks():
    grid=np.array([0.,1.,2.])
    one=np.array([0.,3.,1.])
    ll=np.tile(one,(2,3,4,2,1))
    maps=block_bootstrap_maps(ll,grid,n_bootstrap=17,seed=4)
    assert maps.shape == (17,2)
    assert np.all(maps == 1)

def test_ref_derived_orientation_is_applied_once():
    p=np.array([0.8])
    # Carrying the derived state has the same likelihood whether it is ALT or REF.
    alt_derived=call_loglik(p,np.array([1]),np.array([1]),True,0)[0,0]
    ref_derived=call_loglik(p,np.array([0]),np.array([1]),False,0)[0,0]
    alt_ancestral=call_loglik(p,np.array([0]),np.array([1]),True,0)[0,0]
    ref_ancestral=call_loglik(p,np.array([1]),np.array([1]),False,0)[0,0]
    assert np.isclose(alt_derived, np.log(.8))
    assert np.isclose(ref_derived, np.log(.8))
    assert np.isclose(alt_ancestral, np.log(.2))
    assert np.isclose(ref_ancestral, np.log(.2))
