import json
from argparse import Namespace

import numpy as np
import tskit
from insertion_real_data import (block_bootstrap_maps, call_loglik, call_loglik_grid,
                                 chromosome_layout, parse_epsilon_grid,
                                 load_draw_checkpoint, save_draw_checkpoint,
                                 marginalize_epsilon_independent_chromosomes,
                                 marginalize_independent_chromosomes, merge)

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

def test_epsilon_grid_marginalizes_draws_inside_each_chromosome():
    ll=np.zeros((2,3,2,1,4))
    ll[0,:,1]=-100
    ll[1,:,0]=-100
    observed=marginalize_epsilon_independent_chromosomes(ll)
    assert observed.shape==(3,1,4)
    expected=2*np.log((1+np.exp(-100))/2)
    np.testing.assert_allclose(observed,expected)

def test_merge_streams_epsilon_grid_with_same_result(tmp_path):
    epsilon=np.array([0.,.1]); grid=np.array([0.,10.,20.])
    draw_arrays=[]; parts=[]
    for draw_id in range(2):
        part=tmp_path/f"draw_{draw_id}"; part.mkdir(); parts.append(part)
        values=np.arange(2*2*1*3,dtype=float).reshape(2,2,1,3)+draw_id
        draw_arrays.append(values)
        np.save(part/"epsilon_grid.npy",epsilon)
        np.save(part/"ll_by_chrom_epsilon.npy",values)
        np.save(part/"grid.npy",grid)
        np.save(part/"ll_by_chrom.npy",values[0])
        np.save(part/"used_chr1.npy",np.array([10,20]))
        np.save(part/"used_chr2.npy",np.array([30]))
        (part/"samples.txt").write_text("sample1\n")
        (part/"run.json").write_text(json.dumps(
            {"draw_id":draw_id,"chromosomes":["1","2"],"epsilon":0,
             "epsilon_grid":epsilon.tolist(),"ne":"ne.tsv"}))
    output=tmp_path/"merged"
    merge(Namespace(draw_results=parts,manifest=None,draw_root=None,skip_d0=False,
                    bootstrap_replicates=5,bootstrap_seed=1,output=output))
    expected=marginalize_epsilon_independent_chromosomes(np.stack(draw_arrays))
    np.testing.assert_allclose(np.load(output/"ll_marginal_epsilon.npy"),expected)

def test_argtest_chromosome_offsets_are_read_from_metadata():
    tables=tskit.TableCollection(sequence_length=30)
    tables.metadata_schema=tskit.MetadataSchema({"codec":"json"})
    tables.metadata={"chrom_offsets":[{"chrom":str(i),"offset":3*(i-1),"length":3}
                                                for i in range(1,11)]}
    layout=chromosome_layout(tables.tree_sequence())
    assert layout["1"] == (0.0,3.0)
    assert layout["10"] == (27.0,3.0)

def test_fai_reconstructs_padded_concatenation(tmp_path):
    fai=tmp_path/"ref.fa.fai"
    fai.write_text("1\t10\t0\t0\t0\n2\t20\t0\t0\t0\nscaffold\t5\t0\t0\t0\n")
    tables=tskit.TableCollection(sequence_length=32)
    layout=chromosome_layout(tables.tree_sequence(),fai)
    assert layout=={"1":(0.0,11.0),"2":(11.0,21.0)}

def test_draw_checkpoint_round_trip_and_config_guard(tmp_path):
    checkpoint=tmp_path/"draw.npz"; config={"draw":3,"grid":[0,1]}
    arrays={name:np.full((2,),index,dtype=float) for index,name in enumerate(
        ("ll","ll_epsilon","ll_by_block","sites_by_block","ll_by_d0","sites_by_d0"))}
    counts={"1":{"targets":10,"used":2}}; used={"1":np.array([4,9])}
    save_draw_checkpoint(checkpoint,config,0,7,arrays,counts,used,["1"])
    chrom,row,observed_arrays,observed_counts,observed_used=load_draw_checkpoint(
        checkpoint,config,["1"])
    assert (chrom,row)==(0,7)
    assert observed_counts==counts
    np.testing.assert_array_equal(observed_used["1"],[4,9])
    for name in arrays:
        np.testing.assert_array_equal(observed_arrays[name],arrays[name])
    assert not (tmp_path/"draw.npz.tmp").exists()
    try:
        load_draw_checkpoint(checkpoint,{"draw":4,"grid":[0,1]},["1"])
    except ValueError:
        pass
    else:
        raise AssertionError("incompatible checkpoint configuration accepted")

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

def test_epsilon_grid_matches_repeated_scalar_likelihoods():
    p=np.array([.1,.5,.9])
    observed=np.array([1,0]); called=np.array([1,1])
    epsilon=parse_epsilon_grid("0:0.1:0.05")
    actual=call_loglik_grid(p,observed,called,True,epsilon)
    expected=np.stack([call_loglik(p,observed,called,True,e) for e in epsilon])
    np.testing.assert_allclose(actual,expected)

def test_epsilon_grid_parser_is_inclusive_and_validates_range():
    np.testing.assert_allclose(parse_epsilon_grid("0:0.15:0.05"),[0,.05,.1,.15])
    np.testing.assert_allclose(parse_epsilon_grid("0.01,0.05"),[.01,.05])
    for bad in ("0:.2:0", "0.1,0.1", "0,0.5"):
        try:
            parse_epsilon_grid(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid epsilon grid accepted: {bad}")
