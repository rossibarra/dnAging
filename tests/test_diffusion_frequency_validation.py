from pathlib import Path

import numpy as np

from diffusion_frequency_validation import (
    read_final_trajectories,
    read_sample_mutation_ids,
    sample_panel_focal,
)


def test_reads_all_final_lineages_and_their_trajectories(tmp_path: Path):
    frequency = tmp_path / "x.frequencies.tsv"
    frequency.write_text(
        "replicate\tseed\tgeneration\tallele_state\tmutation_id\t"
        "origin_generation\tpopulation_frequency\n"
        "1\t101\t2\tderived\told\t1\t0.2\n"
        "1\t101\t2\tderived\tfocal\t2\t0.1\n"
        "1\t101\t3\tderived\tfocal\t2\t0.3\n"
        "1\t101\t4\tderived\tfocal\t2\t0.4\n"
        "1\t101\t4\tderived\tother\t3\t0.2\n"
    )
    got = read_final_trajectories(frequency, np.asarray([2, 3]), 4)
    assert [x[:3] for x in got] == [("focal", 2, 0.4), ("other", 3, 0.2)]
    np.testing.assert_allclose(got[0][3], [0.1, 0.3])
    np.testing.assert_allclose(got[1][3], [0.0, 0.0])


def test_joint_panel_sampling_selects_at_most_one_polymorphic_lineage():
    trajectories = [
        ("a", 1, 0.2, np.asarray([0.1])),
        ("b", 2, 0.3, np.asarray([0.2])),
    ]
    first = sample_panel_focal(101, trajectories, 26)
    second = sample_panel_focal(101, trajectories, 26)
    assert first[0][0] in {"a", "b"}
    assert 0 < first[1] < 26
    assert first[0][0] == second[0][0] and first[1] == second[1]


def test_sample_ids_are_sorted_and_panel_draw_is_reproducible(tmp_path: Path):
    sample = tmp_path / "x.samples.tsv"
    sample.write_text(
        "replicate\tseed\tgeneration\thaploid_genotype\tmutation_id\n"
        "1\t101\t3\t1\tfocal\n1\t101\t2\t0\tNA\n"
    )
    replicate, seed, generations, ids = read_sample_mutation_ids(sample)
    assert (replicate, seed) == (1, 101)
    np.testing.assert_array_equal(generations, [2, 3])
    np.testing.assert_array_equal(ids, ["NA", "focal"])
