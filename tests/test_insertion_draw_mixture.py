import json
from argparse import Namespace

import numpy as np
import tskit
import tszip
import insertion_real_data as insertion
from insertion_real_data import (block_bootstrap_maps, call_loglik, call_loglik_grid,
                                 call_loglik_scale_epsilon_grid,
                                 chromosome_layout, parse_epsilon_grid,
                                 parse_time_scale_grid,
                                 joint_accumulator_path, open_joint_accumulator,
                                 load_draw_checkpoint, save_draw_checkpoint,
                                 marginalize_epsilon_independent_chromosomes,
                                 marginalize_scale_epsilon_independent_chromosomes,
                                 marginalize_independent_chromosomes, merge, run_draw)

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
    time_scale=np.array([1.])
    draw_arrays=[]; parts=[]
    for draw_id in range(2):
        part=tmp_path/f"draw_{draw_id}"; part.mkdir(); parts.append(part)
        values=np.arange(2*2*1*3,dtype=float).reshape(2,2,1,3)+draw_id
        draw_arrays.append(values)
        np.save(part/"epsilon_grid.npy",epsilon)
        np.save(part/"ll_by_chrom_epsilon.npy",values)
        np.save(part/"time_scale_grid.npy",time_scale)
        np.save(part/"ll_by_chrom_scale_epsilon.npy",values[None])
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
    expected_joint=marginalize_scale_epsilon_independent_chromosomes(
        np.stack(draw_arrays)[:,None])
    np.testing.assert_allclose(
        np.load(output/"ll_marginal_scale_epsilon.npy"),expected_joint)
    np.testing.assert_allclose(
        np.load(output/"ll_marginal_scale_epsilon.npy")[0],expected)

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
        ("ll","ll_epsilon","ll_scale_epsilon","ll_by_block","sites_by_block",
         "ll_by_d0","sites_by_d0"))}
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

def test_scale_epsilon_chunk_matches_repeated_scalar_likelihoods():
    p=np.array([[.1,.5,.9],[.2,.4,.8]])
    observed=np.array([1,0]); called=np.array([1,1]); epsilon=np.array([0.,.1])
    actual=call_loglik_scale_epsilon_grid(p,observed,called,True,epsilon)
    expected=np.stack([
        np.stack([call_loglik(row,observed,called,True,e) for e in epsilon])
        for row in p])
    np.testing.assert_allclose(actual,expected)

def test_joint_accumulator_is_disk_backed_and_resumable(tmp_path):
    checkpoint=tmp_path/"draw.checkpoint.npz"
    sidecar=joint_accumulator_path(checkpoint); shape=(3,2,1,2,4)
    accumulator=open_joint_accumulator(sidecar,shape)
    assert isinstance(accumulator,np.memmap)
    accumulator[1,1,0,0]=np.arange(4); accumulator.flush(); del accumulator
    resumed=open_joint_accumulator(sidecar,shape,resume=True)
    np.testing.assert_array_equal(resumed[1,1,0,0],np.arange(4))
    assert sidecar.stat().st_size < np.prod(shape)*8+1024

def test_joint_accumulator_rejects_orphan_and_wrong_shape(tmp_path):
    sidecar=tmp_path/"joint.npy"
    accumulator=open_joint_accumulator(sidecar,(2,3)); accumulator.flush(); del accumulator
    for resume,shape in ((False,(2,3)),(True,(3,2))):
        try:
            open_joint_accumulator(sidecar,shape,resume=resume)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid joint sidecar state accepted")

def test_run_draw_promotes_disk_backed_joint_sidecar(tmp_path):
    tables=tskit.TableCollection(sequence_length=10)
    tables.metadata_schema=tskit.MetadataSchema({"codec":"json"})
    tables.metadata={"chrom_offsets":[{"chrom":"1","offset":0,"length":10}]}
    left=tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE,time=0)
    right=tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE,time=0)
    root=tables.nodes.add_row(time=10)
    tables.edges.add_row(0,10,root,left); tables.edges.add_row(0,10,root,right)
    site=tables.sites.add_row(1,"A")
    tables.mutations.add_row(site=site,node=left,derived_state="G")
    tables.sort(); tree_path=tmp_path/"draw.tsz"
    tszip.compress(tables.tree_sequence(),tree_path)
    calls=tmp_path/"calls"; calls.mkdir()
    np.savez(calls/"chr1.npz",samples=np.array(["ancient"]),
             position=np.array([1]),ref=np.array(["A"]),alt=np.array(["G"]),
             observed_alt=np.array([[1]],dtype=np.uint8),
             called=np.array([[1]],dtype=np.uint8))
    ne=tmp_path/"ne.tsv"
    ne.write_text("series\ttime_left\ttime_right\teffective_population_size\n"
                  "constant\t0\t100\t1000\n")
    output=tmp_path/"output"; checkpoint=tmp_path/"checkpoint.npz"
    run_draw(Namespace(tree=tree_path,draw_id=0,calls=calls,chrom_fai=None,
                       ne=ne,ne_series="constant",ne_extrapolate_last=False,
                       epsilon=.01,epsilon_grid="0,0.1",
                       time_scale_grid="1,1.2",time_scale_chunk_size=1,
                       age_min=0,age_max=10,age_step=5,block_size=10,
                       checkpoint=checkpoint,checkpoint_every=1,output=output))
    joint=np.load(output/"ll_by_chrom_scale_epsilon.npy")
    assert joint.shape==(2,2,1,1,3)
    assert not joint_accumulator_path(checkpoint).exists()
    with np.load(checkpoint,allow_pickle=False) as saved:
        assert "ll_scale_epsilon" not in saved.files
    metadata=json.loads((output/"run.json").read_text())
    assert metadata["joint_accumulator"]=="disk-backed numpy memmap"
    assert metadata["time_scale_chunk_size"]==1

def test_joint_checkpoint_restart_matches_uninterrupted_run(tmp_path,monkeypatch):
    tables=tskit.TableCollection(sequence_length=10)
    tables.metadata_schema=tskit.MetadataSchema({"codec":"json"})
    tables.metadata={"chrom_offsets":[{"chrom":"1","offset":0,"length":10}]}
    left=tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE,time=0)
    right=tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE,time=0)
    root=tables.nodes.add_row(time=10)
    tables.edges.add_row(0,10,root,left); tables.edges.add_row(0,10,root,right)
    for position in (1,2):
        site=tables.sites.add_row(position,"A")
        tables.mutations.add_row(site=site,node=left,derived_state="G")
    tables.sort(); tree_path=tmp_path/"draw.tsz"
    tszip.compress(tables.tree_sequence(),tree_path)
    calls=tmp_path/"calls"; calls.mkdir()
    np.savez(calls/"chr1.npz",samples=np.array(["ancient"]),
             position=np.array([1,2]),ref=np.array(["A","A"]),
             alt=np.array(["G","G"]),
             observed_alt=np.array([[1],[0]],dtype=np.uint8),
             called=np.array([[1],[1]],dtype=np.uint8))
    ne=tmp_path/"ne.tsv"
    ne.write_text("series\ttime_left\ttime_right\teffective_population_size\n"
                  "constant\t0\t100\t1000\n")

    def arguments(output,checkpoint):
        return Namespace(tree=tree_path,draw_id=0,calls=calls,chrom_fai=None,
                         ne=ne,ne_series="constant",ne_extrapolate_last=False,
                         epsilon=.01,epsilon_grid="0,0.1",
                         time_scale_grid="1,1.2",time_scale_chunk_size=1,
                         age_min=0,age_max=10,age_step=5,block_size=10,
                         checkpoint=checkpoint,checkpoint_every=1,output=output)

    resumed_output=tmp_path/"resumed"; resumed_checkpoint=tmp_path/"resumed.npz"
    original=insertion.derived_probability_uniform_edge_grid
    calls_seen=0
    def interrupt_on_second_site(*args,**kwargs):
        nonlocal calls_seen
        calls_seen+=1
        if calls_seen==3:
            raise RuntimeError("intentional checkpoint interruption")
        return original(*args,**kwargs)
    monkeypatch.setattr(insertion,"derived_probability_uniform_edge_grid",
                        interrupt_on_second_site)
    try:
        insertion.run_draw(arguments(resumed_output,resumed_checkpoint))
    except RuntimeError as error:
        assert str(error)=="intentional checkpoint interruption"
    else:
        raise AssertionError("test interruption did not occur")
    assert joint_accumulator_path(resumed_checkpoint).exists()
    monkeypatch.setattr(insertion,"derived_probability_uniform_edge_grid",original)
    insertion.run_draw(arguments(resumed_output,resumed_checkpoint))

    full_output=tmp_path/"full"; full_checkpoint=tmp_path/"full.npz"
    insertion.run_draw(arguments(full_output,full_checkpoint))
    for filename in ("ll_by_chrom.npy","ll_by_chrom_epsilon.npy",
                     "ll_by_chrom_scale_epsilon.npy","ll_by_block.npy",
                     "ll_by_chrom_d0.npy"):
        np.testing.assert_array_equal(np.load(resumed_output/filename),
                                      np.load(full_output/filename))

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

def test_time_scale_grid_parser_is_inclusive_and_positive():
    np.testing.assert_allclose(parse_time_scale_grid("0.8:1.2:0.1"),
                               [.8,.9,1.,1.1,1.2])
    np.testing.assert_allclose(parse_time_scale_grid("0.9,1.25"),[.9,1.25])
    for bad in ("0:1:0.1", "1.1,1.1", "1:0.8:0.1", "nan"):
        try:
            parse_time_scale_grid(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid time-scale grid accepted: {bad}")
