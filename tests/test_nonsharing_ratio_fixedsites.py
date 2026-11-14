"""Tests for the additions that make the non-sharing-ratio CLI runnable on real data:
panel-fixed sites, Ne-TSV c(t), degenerate-target flags, pseudo-haploid GT conventions."""
import json

import numpy as np
import pytest

import nonsharing_ratio as nr
import nonsharing_ratio_demography as nd
import nonsharing_ratio_io as nio
from nonsharing_ratio_real_data import add_fixed_sites, main
from test_nonsharing_ratio_real_data import (
    ANC_AGES, LENGTH, N_CTRL, N_PANEL, simulate, write_vcf,
)

HDR = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"


# ------------------------------------------------------------------ Ne TSV ---

def write_ne(path, rows, series=True):
    cols = "series\ttime_left\ttime_right\teffective_population_size" if series else \
        "time_left\ttime_right\teffective_population_size"
    lines = [cols]
    for s, a, b, n, *rep in rows:
        lines.append("\t".join(map(str, ([s] if series else []) + [a, b, n])))
    path.write_text("\n".join(lines) + "\n")


def test_ne_tsv_constant_and_epoch_average(tmp_path):
    p = tmp_path / "ne.tsv"
    write_ne(p, [("posterior_mean", 0, 100, 1000), ("posterior_mean", 100, 300, 4000)])
    edges = nd.make_edges(300, 3)          # bins 0-100, 100-200, 200-300, then open
    e, r = nd.rates_from_ne_tsv(p, edges)
    assert np.allclose(r[:3], [1 / 2000, 1 / 8000, 1 / 8000])
    assert r[-1] == pytest.approx(1 / 8000) and np.isinf(e[-1])
    # a bin straddling the epoch boundary gets the time-weighted average
    e2, r2 = nd.rates_from_ne_tsv(p, nd.make_edges(200, 1))
    assert r2[0] == pytest.approx(0.5 * (1 / 2000 + 1 / 8000))
    # ploidy=1 gives 1/N
    assert nd.rates_from_ne_tsv(p, edges, ploidy=1)[1][0] == pytest.approx(1 / 1000)
    # last epoch is extended beyond the file's final boundary
    assert nd.rates_from_ne_tsv(p, nd.make_edges(1000, 2))[1][-2] == pytest.approx(1 / 8000)


def test_ne_tsv_series_handling(tmp_path):
    p = tmp_path / "ne.tsv"
    write_ne(p, [("a", 0, 10, 100), ("b", 0, 10, 200)])
    with pytest.raises(ValueError, match="several series"):
        nd.rates_from_ne_tsv(p, nd.make_edges(10, 1))
    assert nd.rates_from_ne_tsv(p, nd.make_edges(10, 1), series="b")[1][0] == pytest.approx(1 / 400)
    with pytest.raises(ValueError, match="not in"):
        nd.rates_from_ne_tsv(p, nd.make_edges(10, 1), series="zzz")
    q = tmp_path / "rep.tsv"
    q.write_text("series\treplicate_index\ttime_left\ttime_right\teffective_population_size\n"
                 "replicate\t0\t0\t10\t100\nreplicate\t1\t0\t10\t200\n")
    with pytest.raises(ValueError, match="several replicates"):
        nd.rates_from_ne_tsv(q, nd.make_edges(10, 1), series="replicate")
    plain = tmp_path / "plain.tsv"
    write_ne(plain, [(None, 0, 10, 50)], series=False)
    assert nd.rates_from_ne_tsv(plain, nd.make_edges(10, 1))[1][0] == pytest.approx(1 / 100)


# -------------------------------------------------------------------- flags ---

def test_target_flags():
    f = np.array([0.0, 0.2, 0.05, 0.3, np.nan])
    with np.errstate(all="ignore"):
        H = np.log(f / 0.1)
    flags = nr.target_flags(f, H)
    assert list(flags) == ["f0", "", "H<=0", "", "no_data"]
    # zero denominator is no_data even if f was forced finite
    assert nr.target_flags([0.2], [0.5], den=[0])[0] == "no_data"
    # T_hat for a flagged target is reported as 0 by invert_hazard: the flag is the only signal
    num = np.array([[0.0, 5.0, 10.0, 10.0]]); den = np.array([[10.0, 10.0, 10.0, 10.0]])
    T, Hh, ff = nr.ages_from_counts(num, den, np.array([False, False, True, True]),
                                    np.array([0.0, 1000.0, np.inf]), np.array([1e-3, 1e-3]))
    assert T[0] == 0.0 and list(nr.target_flags(ff, Hh)) == ["f0", "H<=0"]


# ------------------------------------------------------ pseudo-haploid GT ---

def test_pseudohaploid_gt_conventions(tmp_path):
    """Real ancient data: diploid-coded 0/0, 1/1, ./. (and some pipelines write haploid 0/1/.)."""
    p = tmp_path / "ph.vcf"
    p.write_text(HDR + "s1\ts2\ts3\n"
                 "1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/0\t1/1\t./.\n"
                 "1\t20\t.\tA\tG\t.\t.\t.\tGT\t0\t1\t.\n"
                 "1\t30\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/.\t./1\n"
                 "1\t40\t.\tC\t.\t.\t.\t.\tGT\t0/0\t0/0\t./.\n")   # monomorphic: ALT '.'
    v = nio.read_vcf(p, mode="sample")
    assert v.chroms["1"]["gt"].tolist() == [[0, 0, -1], [1, 1, -1], [-1, -1, -1]]
    assert v.skipped["non_snp_or_multiallelic"] == 1            # the ALT='.' record
    # polarise: REF=ancestral A, ALT=derived G at all three sites
    pos = np.array([10.0, 20.0, 30.0])
    out, st = nio.polarise(v, {"1": (0, 100)}, pos, np.array(["A"] * 3, object),
                           np.array(["G"] * 3, object))
    assert out.tolist() == [[0, 0, -1], [1, 1, -1], [-1, -1, -1]]
    assert st["matched"] == 3
    # flipped REF/ALT relative to the ARG gives the same derived calls
    q = tmp_path / "ph_flip.vcf"
    q.write_text(HDR + "s1\ts2\ts3\n1\t10\t.\tG\tA\t.\t.\t.\tGT\t1/1\t0/0\t./.\n")
    vq = nio.read_vcf(q, mode="sample")
    oq, _ = nio.polarise(vq, {"1": (0, 100)}, pos[:1], np.array(["A"], object), np.array(["G"], object))
    assert oq.tolist() == [[0], [1], [-1]]
    # position shift is applied to the VCF coordinate
    o2, st2 = nio.polarise(v, {"1": (0, 100)}, pos + 1, np.array(["A"] * 3, object),
                           np.array(["G"] * 3, object), position_shift=1)
    assert st2["matched"] == 3 and (o2 == out).all()
    assert nio.polarise(v, {"1": (0, 100)}, pos + 1, np.array(["A"] * 3, object),
                        np.array(["G"] * 3, object), position_shift=0)[1]["matched"] == 0


# -------------------------------------------------- scan / polarise fixed ---

def test_scan_panel_fixed_and_polarise(tmp_path):
    p = tmp_path / "panel.vcf"
    cols = [f"p{i}" for i in range(5)]
    rows = [
        ("1", 10, "A", "G", ["1", "1", "1", "1", "1"]),    # fixed ALT
        ("1", 20, "A", "G", ["1", "1", "1", "1", "."]),    # fixed ALT with one missing
        ("1", 30, "A", "G", ["1", "0", "1", "1", "1"]),    # segregating
        ("1", 40, "A", "G", ["0", "0", "0", "0", "0"]),    # fixed REF
        ("1", 50, "A", "G", ["1", "1", "1", ".", "."]),    # too few called
        ("1", 60, "A", "GT", ["1"] * 5),                   # not a SNP
        ("2", 5, "C", "T", ["1/1", "1|1", "1/1", "1/1", "1/1"]),  # diploid-coded: slow path
    ]
    with p.open("w") as fh:
        fh.write(HDR + "\t".join(cols) + "\n")
        for c, pos, ref, alt, g in rows:
            fh.write("\t".join([c, str(pos), ".", ref, alt, ".", ".", ".", "GT"] + g) + "\n")
    fx = nio.scan_panel_fixed(p, min_called=4)
    assert fx["pos"].tolist() == [10, 20, 40, 5]
    assert fx["fixed_idx"].tolist() == [1, 1, 0, 1]
    assert fx["n_called"].tolist() == [5, 4, 5, 5] and fx["n_columns"] == 5
    only1 = nio.scan_panel_fixed(p, chroms=["1"], min_called=4)
    assert only1["pos"].tolist() == [10, 20, 40]

    layout = {"1": (0.0, 100.0), "2": (101.0, 100.0)}
    # position 10 is already an ARG site and must not be added twice
    ex = nio.polarise_fixed_sites(fx, layout, np.array([10.0]), "ref")
    assert ex["positions"].tolist() == [20.0, 106.0]                  # 5 + offset 101
    assert ex["stats"]["fixed_ancestral_dropped"] == 1 and ex["stats"]["already_in_arg"] == 1
    assert ex["derived"].tolist() == ["G", "T"] and ex["n_called"].tolist() == [4, 5]
    # table polarity: ancestral G at pos 20 means the fixed ALT (G) is ancestral -> dropped;
    # pos 5 has no entry -> dropped; pos 10 table entry irrelevant (already in ARG)
    tab = tmp_path / "pol.txt"
    tab.write_text("# chrom pos anc\n1 20 G\n1 40 G\n")
    table = nio.read_polarity_table(tab)
    ex2 = nio.polarise_fixed_sites(fx, layout, np.array([]), "table", table)
    assert ex2["stats"]["no_polarity"] == 2 and ex2["stats"]["kept"] == 1
    assert ex2["positions"].tolist() == [40.0]    # fixed REF A is derived when ancestral is G
    assert ex2["ancestral"].tolist() == ["G"] and ex2["derived"].tolist() == ["A"]


def test_add_fixed_sites_orders_and_marks():
    pos = np.array([10.0, 30.0])
    anc = np.array(["A", "A"], object); der = np.array(["G", "G"], object)
    panel = np.array([[1, 0], [0, 0], [1, 1]], dtype=np.int8)
    ages = np.array([5.0, 7.0])
    extra = {"positions": np.array([20.0]), "ancestral": np.array(["C"], object),
             "derived": np.array(["T"], object), "n_called": np.array([2], dtype=np.int32)}
    p2, a2, d2, pan2, ag2, ex, ref_count = add_fixed_sites(
        pos, anc, der, panel, ages, extra)
    assert p2.tolist() == [10.0, 20.0, 30.0] and ex.tolist() == [False, True, False]
    assert pan2[:, 1].tolist() == [1, 1, 1] and ref_count.tolist() == [2, 2, 1]
    assert np.isinf(ag2[1])
    assert d2.tolist() == ["G", "T", "G"] and a2.tolist() == ["A", "C", "A"]
    # The site is fixed among called references, is weighted by the two calls, and obeys
    # both the fixed-everywhere exclusion and include_panel_fixed=False.
    mask = nr.select_sites(pan2, np.zeros((1, 3), int), ag2, a0=100.0)
    assert mask.tolist() == [False, True, False]
    assert not nr.select_sites(pan2, np.ones((1, 3), int), ag2, a0=100.0).any()
    assert not nr.select_sites(pan2, np.zeros((1, 3), int), ag2, a0=100.0,
                               include_panel_fixed=False).any()
    num, den = nr.nonsharing_counts(pan2, np.zeros((1, 3), int), mask,
                                    reference_counts=ref_count)
    assert num[0, 0] == 2 and den[0, 0] == 2


# ----------------------------------------------- end to end with the CLI ---

@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = tmp_path_factory.mktemp("nsrfix")
    ts = simulate(seed=11)
    p = N_PANEL
    full = ts.simplify(list(range(p)), filter_sites=False)          # ARG with every site
    geno = full.genotype_matrix()
    fixed_ids = np.nonzero(geno.sum(axis=1) == p)[0]
    assert fixed_ids.size >= 5, "simulation should yield panel-fixed derived sites"
    arg_no_fixed = full.delete_sites(fixed_ids)
    full.dump(str(d / "full.trees")); arg_no_fixed.dump(str(d / "nofixed.trees"))
    ctrl = list(range(p, p + N_CTRL)); anc = list(range(p + N_CTRL, p + N_CTRL + len(ANC_AGES)))
    write_vcf(d / "controls.vcf.gz", ts, ctrl, [f"c{i}" for i in range(N_CTRL)])
    write_vcf(d / "ancient.vcf.gz", ts, anc, [f"a{i}" for i in range(len(ANC_AGES))])
    write_vcf(d / "panel.vcf.gz", ts, list(range(p)), [f"p{i}" for i in range(p)])
    (d / "g.fai").write_text(f"chr1\t{LENGTH}\t0\t0\t0\n")
    return d, fixed_ids


def run(d, tree, out, *extra):
    main(["run-draw", "--tree", str(d / tree), "--draw-id", "0",
          "--controls-vcf", str(d / "controls.vcf.gz"), "--ancient-vcf", str(d / "ancient.vcf.gz"),
          "--controls-mode", "sample", "--chrom-fai", str(d / "g.fai"),
          "--a0", "3000", "--block-size", "100000", "--c-grid-max", "6000", "--c-grid-n", "30",
          "--output", str(out), *extra])
    return np.load(out)


def test_panel_fixed_sites_recover_the_full_arg_result(world, tmp_path):
    d, fixed_ids = world
    ref = run(d, "full.trees", tmp_path / "full.npz")
    without = run(d, "nofixed.trees", tmp_path / "nofixed.npz")
    main(["prepare-fixed", "--panel-vcf", str(d / "panel.vcf.gz"), "--min-called", str(N_PANEL),
          "--output", str(tmp_path / "fixed.npz")])
    withfx = run(d, "nofixed.trees", tmp_path / "withfixed.npz",
                 "--panel-fixed", str(tmp_path / "fixed.npz"), "--fixed-polarity", "ref")
    meta = json.loads(str(withfx["meta"]))
    assert meta["n_selected_extra_fixed_sites"] == fixed_ids.size
    assert meta["panel_fixed_input"]["kept"] == fixed_ids.size
    # dropping the fixed sites changes the counts; adding them back restores the full-ARG counts
    assert not np.allclose(without["num"].sum(0), ref["num"].sum(0))
    assert np.allclose(withfx["num"], ref["num"]) and np.allclose(withfx["den"], ref["den"])
    assert (without["den"] <= withfx["den"] + 1e-9).all()


def test_ne_option_sets_rates(world, tmp_path):
    d, _ = world
    ne = tmp_path / "ne.tsv"
    write_ne(ne, [("posterior_mean", 0, 1000, 2000), ("posterior_mean", 1000, 5000, 8000)])
    res = run(d, "full.trees", tmp_path / "ne.npz", "--ne", str(ne), "--ne-series", "posterior_mean")
    assert json.loads(str(res["meta"]))["c_source"] == str(ne)
    assert res["rates"][0] == pytest.approx(1 / 4000) and res["rates"][-1] == pytest.approx(1 / 16000)
    with pytest.raises(ValueError, match="at most one"):
        run(d, "full.trees", tmp_path / "bad.npz", "--ne", str(ne), "--demography", str(ne))


def test_chroms_restriction_and_flag_output(world, tmp_path):
    d, _ = world
    res = run(d, "full.trees", tmp_path / "sub.npz", "--chroms", "chr1")
    assert json.loads(str(res["meta"]))["chroms"] == ["chr1"]
    assert "flags" in res.files
    merged = tmp_path / "m.tsv"
    main(["merge", "--draw-results", str(tmp_path / "sub.npz"), "--bootstrap-replicates", "20",
          "--output", str(merged)])
    header = merged.read_text().split("\n")[0].split("\t")
    assert header[-1] == "flag"
    with pytest.raises(ValueError, match="not in the ARG layout"):
        run(d, "full.trees", tmp_path / "sub2.npz", "--chroms", "chrZ")


# --------------------------------------- vectorised ARG helpers vs reference ---

def test_vectorised_arg_helpers_match_reference():
    import msprime
    ts = simulate(seed=5)
    panel_ts = ts.simplify(list(range(N_PANEL)), filter_sites=False)
    # add recurrent / multi-allelic sites and a back mutation so that edge cases are covered
    tables = panel_ts.dump_tables()
    sid = tables.sites.add_row(position=123_457, ancestral_state="0")
    tables.mutations.add_row(site=sid, node=0, derived_state="1")
    tables.mutations.add_row(site=sid, node=1, derived_state="2")        # 3 alleles -> -1
    sid2 = tables.sites.add_row(position=123_459, ancestral_state="0")
    tables.mutations.add_row(site=sid2, node=2, derived_state="1")
    tables.mutations.add_row(site=sid2, node=3, derived_state="1")       # recurrent, same state
    sid3 = tables.sites.add_row(position=130_001, ancestral_state="0")
    tables.mutations.add_row(site=sid3, node=int(panel_ts.at(130_001).root), derived_state="1")  # on the root
    tables.sort(); tables.build_index(); tables.compute_mutation_parents()
    panel_ts = tables.tree_sequence()
    assert np.array_equal(nr.panel_derived_matrix(panel_ts), nr.panel_derived_matrix_slow(panel_ts))
    fast, slow = nr.edge_midpoint_ages(panel_ts), nr.edge_midpoint_ages_slow(panel_ts)
    assert np.array_equal(np.isnan(fast), np.isnan(slow))
    assert np.allclose(fast[~np.isnan(fast)], slow[~np.isnan(slow)])
    assert np.isinf(fast).sum() == np.isinf(slow).sum() > 0
    # vectorised allele strings equal the per-site definition
    anc, der = nio.arg_site_alleles(panel_ts)
    for site in panel_ts.sites():
        states = {m.derived_state for m in site.mutations}
        want = next(iter(states)) if len(states) == 1 and next(iter(states)) != site.ancestral_state else ""
        assert der[site.id] == want and anc[site.id] == site.ancestral_state


def test_empty_panel_fixed_list_is_a_noop(world, tmp_path):
    """The real panel VCF lists only variable sites, so prepare-fixed yields no sites."""
    d, _ = world
    var_only = tmp_path / "panel_var.vcf"
    var_only.write_text(HDR + "p0\tp1\n1\t10\t.\tA\tG\t.\t.\t.\tGT\t0\t1\n")
    main(["prepare-fixed", "--panel-vcf", str(var_only), "--output", str(tmp_path / "none.npz")])
    ref = run(d, "full.trees", tmp_path / "a.npz")
    res = run(d, "full.trees", tmp_path / "b.npz", "--panel-fixed", str(tmp_path / "none.npz"),
               "--fixed-polarity", "ref")
    assert np.allclose(ref["num"], res["num"]) and np.allclose(ref["den"], res["den"])
    # the default polarity needs a table when fixed sites exist
    with pytest.raises(ValueError, match="polarity-table"):
        main(["prepare-fixed", "--panel-vcf", str(d / "panel.vcf.gz"), "--min-called", str(N_PANEL),
              "--output", str(tmp_path / "some.npz")])
        run(d, "nofixed.trees", tmp_path / "c.npz", "--panel-fixed", str(tmp_path / "some.npz"))


def test_chromosome_name_styles_are_interchangeable(world, tmp_path):
    """Real ARG layout uses '1'; a VCF written as 'chr1' must give identical counts."""
    d, _ = world
    (tmp_path / "num.fai").write_text(f"1\t{LENGTH}\t0\t0\t0\n")
    a = run(d, "full.trees", tmp_path / "chr.npz")                       # fai 'chr1', VCF 'chr1'
    main(["run-draw", "--tree", str(d / "full.trees"), "--draw-id", "0",
          "--controls-vcf", str(d / "controls.vcf.gz"), "--ancient-vcf", str(d / "ancient.vcf.gz"),
          "--controls-mode", "sample", "--chrom-fai", str(tmp_path / "num.fai"),
          "--a0", "3000", "--block-size", "100000", "--c-grid-max", "6000", "--c-grid-n", "30",
          "--output", str(tmp_path / "num.npz")])
    b = np.load(tmp_path / "num.npz")
    assert np.allclose(a["num"], b["num"]) and np.allclose(a["den"], b["den"]) and a["den"].sum() > 0
