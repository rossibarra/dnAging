import csv
from pathlib import Path

import numpy as np

from direct_frequency_age_infer import (
    add_bernoulli_log_likelihood,
    normalize_log_likelihood,
    read_derived_frequency,
    read_sample_genotypes,
)


def test_reads_binary_genotype_and_collapsed_derived_frequency(tmp_path: Path):
    frequency = tmp_path / "replicate_000001.seed_101.frequencies.tsv"
    frequency.write_text(
        "replicate\tseed\tgeneration\tallele_state\tmutation_id\t"
        "origin_generation\tpopulation_frequency\n"
        "1\t101\t1\tancestral\tNA\tNA\t0.75\n"
        "1\t101\t1\tderived\t0\t1\t0.10\n"
        "1\t101\t1\tderived\t1\t1\t0.15\n"
        "1\t101\t2\tancestral\tNA\tNA\t0.25\n"
        "1\t101\t2\tderived\t1\t1\t0.75\n"
    )
    sample = tmp_path / "replicate_000001.seed_101.samples.tsv"
    sample.write_text(
        "replicate\tseed\tgeneration\thaploid_genotype\tmutation_id\n"
        "1\t101\t2\t1\t1\n"
        "1\t101\t1\t0\tNA\n"
    )

    np.testing.assert_allclose(read_derived_frequency(frequency, 2), [0.25, 0.75])
    replicate, seed, generations, genotypes = read_sample_genotypes(sample)
    assert (replicate, seed) == (1, 101)
    np.testing.assert_array_equal(generations, [1, 2])
    np.testing.assert_array_equal(genotypes, [0, 1])


def test_known_frequency_bernoulli_likelihood_and_normalization():
    log_likelihood = np.zeros((2, 3))
    frequency = np.asarray([0.0, 0.25, 1.0])
    add_bernoulli_log_likelihood(log_likelihood, frequency, np.asarray([0, 1]))

    np.testing.assert_allclose(log_likelihood[0, :2], [0.0, np.log(0.75)])
    assert np.isneginf(log_likelihood[0, 2])
    assert np.isneginf(log_likelihood[1, 0])
    np.testing.assert_allclose(log_likelihood[1, 1:], [np.log(0.25), 0.0])

    posterior = normalize_log_likelihood(log_likelihood)
    np.testing.assert_allclose(posterior.sum(axis=1), 1.0)
    assert np.argmax(posterior[0]) == 0
    assert np.argmax(posterior[1]) == 2
