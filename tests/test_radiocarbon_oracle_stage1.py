from argparse import Namespace

import numpy as np

import radiocarbon_oracle_stage1 as oracle


def test_oracle_merge_builds_calibration_result(tmp_path):
    parts=tmp_path/"parts"; parts.mkdir()
    grid=np.array([0.,500.,1000.]); samples=np.array(["modern_age","dated1","dated2"])
    truth=np.array([0.,500.,1000.])
    for epsilon in (0.,.01,.05):
        ll=np.array([[-1.,-2.,-3.],[-3.,0.,-3.],[-3.,-2.,0.]])-20*epsilon
        np.savez_compressed(parts/f"epsilon_{epsilon:.3f}.npz",
            fit_epsilon=epsilon,Tgrid=grid,samples=samples,true_ages=truth,
            log_likelihood=ll,injected_epsilon=.01,injection_seed=7,
            flipped=3,allele_total=300,counts_json="{}")
    output=tmp_path/"result"
    oracle.merge(Namespace(parts=parts,rc_sd=100,n_folds=2,
                           bootstrap_replicates=20,bootstrap_seed=4,output=output))
    assert (output/"likelihood"/"ll_marginal_epsilon.npy").exists()
    assert (output/"calibration"/"cross_validated_ages.tsv").exists()
    assert (output/"oracle_run.json").exists()
