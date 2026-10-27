import csv
import json
from argparse import Namespace

import numpy as np

import radiocarbon_calibration as cal


def test_radiocarbon_profile_integrates_density_not_midpoint():
    age=np.array([0.,10.,20.])
    ll=np.array([[0.,-10.,0.],[-4.,0.,-4.]])
    rc_age=np.array([0.,10.,20.])
    rc_log_weight=cal.normalized_log_weights(rc_age,np.array([1.,0.,1.]))
    profile=cal.radiocarbon_log_profile(age,ll,rc_age,rc_log_weight)
    assert profile[0]>profile[1]


def test_cross_validation_and_unknown_age_inference(tmp_path):
    likelihood=tmp_path/"likelihood"; likelihood.mkdir()
    age=np.array([0.,10.,20.]); epsilon=np.array([0.,.1])
    samples=["dated_young","dated_old","unknown"]
    ll=np.full((2,3,3),-8.)
    ll[0,0]=[0.,-4.,-8.]
    ll[0,1]=[-8.,-4.,0.]
    ll[0,2]=[-5.,0.,-5.]
    ll[1]=ll[0]-3
    np.save(likelihood/"grid.npy",age)
    np.save(likelihood/"epsilon_grid.npy",epsilon)
    np.save(likelihood/"ll_marginal_epsilon.npy",ll)
    (likelihood/"samples.txt").write_text("\n".join(samples)+"\n")
    radiocarbon=tmp_path/"radiocarbon.tsv"
    radiocarbon.write_text(
        "sample\tage_generations\tdensity\n"
        "dated_young\t0\t1\n"
        "dated_old\t20\t1\n")
    output=tmp_path/"output"
    cal.run(Namespace(likelihood_dir=likelihood,radiocarbon=radiocarbon,
                      folds=None,n_folds=0,unknown_samples=None,
                      age_prior=None,bootstrap_replicates=100,
                      bootstrap_seed=7,output=output))
    metadata=json.loads((output/"run.json").read_text())
    assert metadata["n_dated"]==2
    assert metadata["n_unknown"]==1
    assert metadata["epsilon"]["mle"]==0
    with (output/"cross_validated_ages.tsv").open() as handle:
        rows=list(csv.DictReader(handle,delimiter="\t"))
    assert {row["sample"] for row in rows}=={"dated_young","dated_old"}
    assert {float(row["map_T"]) for row in rows}=={0.,20.}
    with (output/"unknown_ages.tsv").open() as handle:
        unknown=list(csv.DictReader(handle,delimiter="\t"))
    assert len(unknown)==1 and unknown[0]["sample"]=="unknown"
    assert float(unknown[0]["map_T"])==10.
