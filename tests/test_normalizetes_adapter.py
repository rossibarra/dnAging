from types import SimpleNamespace

import numpy as np

import posterior_sample_age_infer as inf


def test_chunk_sites_accepts_current_normalizetes_vcfchunk_fields():
    chunk = SimpleNamespace(
        chrom=np.array(["simulation", "simulation"]),
        position=np.array([10, 20]),
        ref_index=np.array([0, 2], dtype=np.uint8),
        alt_index=np.array([1, 3], dtype=np.uint8),
    )

    chrom, position, ref, alt = inf._chunk_sites(chunk)

    assert chrom.tolist() == ["simulation", "simulation"]
    assert position.tolist() == [10, 20]
    assert ref.tolist() == [0, 2]
    assert alt.tolist() == [1, 3]


def test_vcf_readers_ignore_normalizetes_header_and_summary_events(monkeypatch):
    def read_chunks(*_args, **_kwargs):
        yield ["sample"], None, None, None
        yield ["sample"], None, {"records": 0}, "digest"

    monkeypatch.setattr(inf, "_import_repo", lambda: (None, None, read_chunks, 255))

    order, calls = inf.read_ancient(["empty.vcf"], None, "simulation", None,
                                    100, True)
    panel = inf.read_panel_alt(["empty.vcf"], "simulation", 100, True, 26)

    assert order == []
    assert calls == {}
    assert panel == {}
