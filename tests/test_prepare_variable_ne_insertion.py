import gzip

from prepare_variable_ne_insertion import numbered_draws, read_calls


def test_numbered_draws_are_numeric_and_ignore_other_files(tmp_path):
    for name in ("run.combined.99.tsz", "run.combined.9.tsz", "notes.tsz"):
        (tmp_path / name).touch()
    assert [x[0] for x in numbered_draws(tmp_path)] == [9, 99]


def test_read_calls_uses_singer_position_map_and_ancient_columns(tmp_path):
    mapping = tmp_path / "position_map.tsv"
    mapping.write_text("site_id\toriginal_pos\tsinger_pos\n7\t100\t103\n")
    vcf = tmp_path / "all.vcf.gz"
    with gzip.open(vcf, "wt") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmodern_01\tancient_01\tancient_02\n")
        out.write("1\t100\t7\tA\tG\t.\tPASS\t.\tGT\t1\t0\t1\n")
    samples, calls = read_calls(vcf, mapping)
    assert samples == ["ancient_01", "ancient_02"]
    assert list(calls) == [103]
    assert calls[103][0:2] == ("A", "G")
    assert calls[103][2].tolist() == [0, 1]
