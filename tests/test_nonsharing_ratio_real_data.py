import gzip
import json

import msprime
import numpy as np
import pytest
import tskit

import nonsharing_ratio as nr
import nonsharing_ratio_io as nio
from nonsharing_ratio_real_data import main

N_PANEL, N_CTRL = 10, 3
ANC_AGES = [200.0, 800.0]
LENGTH = 400_000


def simulate(seed=3):
    sets = [msprime.SampleSet(N_PANEL, time=0, ploidy=1), msprime.SampleSet(N_CTRL, time=0, ploidy=1)]
    sets += [msprime.SampleSet(1, time=a, ploidy=1) for a in ANC_AGES]
    ts = msprime.sim_ancestry(samples=sets, ploidy=1, population_size=2000,
                              sequence_length=LENGTH, recombination_rate=1e-8,
                              discrete_genome=True, random_seed=seed)
    return msprime.sim_mutations(ts, rate=5e-8, model=msprime.BinaryMutationModel(),
                                 discrete_genome=True, random_seed=seed + 1)


def write_vcf(path, ts, columns, names, flip_every=0, missing_every=0, chrom="chr1"):
    """Haploid VCF with REF=0/ALT=1 ('0'/'1' alleles); flip_every swaps REF/ALT on some rows."""
    geno = ts.genotype_matrix()
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
        for i, site in enumerate(ts.sites()):
            row = geno[i, columns].astype(int)
            ref, alt = "0", "1"
            if flip_every and i % flip_every == 0:
                ref, alt, row = "1", "0", 1 - row
            gts = [str(g) for g in row]
            if missing_every and i % missing_every == 0:
                gts[0] = "."
            fh.write("\t".join([chrom, str(int(site.position)), ".", ref, alt, ".", ".", ".", "GT"] + gts) + "\n")


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("nsr")
    ts = simulate()
    panel_ids = list(range(N_PANEL))
    panel_ts = ts.simplify(panel_ids, filter_sites=False)
    arg = d / "draw0.trees"
    panel_ts.dump(str(arg))
    ctrl_cols = list(range(N_PANEL, N_PANEL + N_CTRL))
    anc_cols = list(range(N_PANEL + N_CTRL, N_PANEL + N_CTRL + len(ANC_AGES)))
    write_vcf(d / "controls.vcf.gz", ts, ctrl_cols, [f"c{i}" for i in range(N_CTRL)], flip_every=3)
    write_vcf(d / "ancient.vcf.gz", ts, anc_cols, [f"a{i}" for i in range(len(ANC_AGES))],
              flip_every=4, missing_every=5)
    (d / "g.fai").write_text(f"chr1\t{LENGTH}\t0\t0\t0\n")
    return d, ts, panel_ts


def test_read_vcf_modes(tmp_path):
    p = tmp_path / "x.vcf"
    p.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
                 "1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1/1\n"
                 "1\t20\t.\tA\tG\t.\t.\t.\tGT\t0/1\t./.\n"
                 "1\t30\t.\tA\tGT\t.\t.\t.\tGT\t0|0\t0|0\n")
    s = nio.read_vcf(p, mode="sample")
    assert s.columns == ["s1", "s2"] and s.skipped["non_snp_or_multiallelic"] == 1
    assert s.chroms["1"]["gt"].tolist() == [[-1, -1], [1, -1]]
    # phased mode rejects any record with an unphased genotype instead of guessing
    ph = nio.read_vcf(p, mode="phased")
    assert ph.skipped["bad_gt"] == 2 and not ph.chroms
    q = tmp_path / "phased.vcf"
    q.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
                 "1\t10\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|1\n"
                 "1\t20\t.\tA\tG\t.\t.\t.\tGT\t.|1\t0|0\n")
    ph = nio.read_vcf(q, mode="phased")
    assert ph.columns == ["s1_1", "s1_2", "s2_1", "s2_2"]
    assert ph.chroms["1"]["gt"].tolist() == [[0, -1], [1, 1], [1, 0], [1, 0]]


def test_polarise_handles_flipped_ref_alt(data):
    d, ts, panel_ts = data
    ctrl = nio.read_vcf(d / "controls.vcf.gz", mode="sample")
    anc, der = nio.arg_site_alleles(panel_ts)
    pos = panel_ts.tables.sites.position
    out, stats = nio.polarise(ctrl, {"chr1": (0.0, LENGTH)}, pos, anc, der)
    truth = ts.genotype_matrix()[:, N_PANEL:N_PANEL + N_CTRL].T
    usable = der != ""
    assert (out[:, usable] == truth[:, usable]).all()
    assert stats["allele_mismatch"] == 0


def test_run_draw_matches_direct_computation(data, tmp_path):
    d, ts, panel_ts = data
    out = tmp_path / "draw0.npz"
    main(["run-draw", "--tree", str(d / "draw0.trees"), "--draw-id", "0",
          "--controls-vcf", str(d / "controls.vcf.gz"), "--ancient-vcf", str(d / "ancient.vcf.gz"),
          "--controls-mode", "sample", "--chrom-fai", str(d / "g.fai"),
          "--a0", "3000", "--block-size", "100000", "--c-grid-max", "6000", "--c-grid-n", "30",
          "--output", str(out)])
    res = np.load(out)
    assert list(res["targets"]) == ["a0", "a1", "c0", "c1", "c2"]
    assert res["num"].shape == (4, 5)  # 400 kb / 100 kb blocks

    # independent computation from the true genotype matrix
    geno = ts.genotype_matrix()
    panel = nr.panel_derived_matrix(panel_ts)
    anc, der = nio.arg_site_alleles(panel_ts)
    usable = (panel >= 0).all(axis=0) & (der != "")
    targets = geno[:, N_PANEL:].T[:, usable]
    ages = nr.edge_midpoint_ages(panel_ts)
    # reproduce missing calls: ancient a0 missing at every 5th ARG site index
    site_idx = np.nonzero(usable)[0]
    called_full = np.ones((N_CTRL + len(ANC_AGES), usable.sum()), bool)
    called_full[len(ANC_AGES) - 2, :] = np.array([i % 5 != 0 for i in site_idx])
    order = list(range(N_PANEL + N_CTRL, N_PANEL + N_CTRL + 2)) + list(range(N_PANEL, N_PANEL + N_CTRL))
    tg = geno[:, order].T[:, usable]
    mask = nr.select_sites(panel[:, usable], tg, ages[usable], 3000.0, True, called_full)
    num, den = nr.nonsharing_counts(panel[:, usable], tg, mask, called_full)
    assert np.allclose(res["num"].sum(0), num[0]) and np.allclose(res["den"].sum(0), den[0])
    meta = json.loads(str(res["meta"]))
    assert meta["n_selected_sites"] == int(mask.sum())

    # merge works and returns one row per ancient sample
    merged = tmp_path / "merged.tsv"
    main(["merge", "--draw-results", str(out), "--bootstrap-replicates", "50", "--output", str(merged)])
    rows = merged.read_text().strip().split("\n")
    assert len(rows) == 3 and rows[1].startswith("a0\t")


def test_overlapping_sample_names_are_rejected(data, tmp_path):
    d, ts, panel_ts = data
    names = tmp_path / "panel.txt"
    names.write_text("c0 p1 p2\n")
    with pytest.raises(ValueError, match="share sample names"):
        main(["run-draw", "--tree", str(d / "draw0.trees"), "--draw-id", "0",
              "--controls-vcf", str(d / "controls.vcf.gz"), "--ancient-vcf", str(d / "ancient.vcf.gz"),
              "--controls-mode", "sample", "--chrom-fai", str(d / "g.fai"), "--panel-names", str(names),
              "--a0", "3000", "--output", str(tmp_path / "x.npz")])
