import gzip
from pathlib import Path

from msprime_exact_time_validation import seeds_for_replicate, vcf_calls


def test_seed_streams_are_distinct_and_reproducible():
    assert seeds_for_replicate(7, 1000) == (1007, 101007, 201007)
    assert len(set(seeds_for_replicate(7, 1000))) == 3
    assert seeds_for_replicate(7, 1000) != seeds_for_replicate(8, 1000)


def test_vcf_calls_reads_26_haploid_modern_and_one_ancient(tmp_path: Path):
    path = tmp_path / "calls.vcf.gz"
    names = [f"modern_{i:02d}" for i in range(1, 27)] + ["ancient_01"]
    genotypes = ["1"] * 3 + ["0"] * 23 + ["1"]
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                 "FILTER", "INFO", "FORMAT", *names]) + "\n")
        handle.write("\t".join(["1", "10", "42", "A", "T", ".", "PASS", ".",
                                 "GT", *genotypes]) + "\n")
    assert vcf_calls(path, 26) == {42: (3, 1)}
