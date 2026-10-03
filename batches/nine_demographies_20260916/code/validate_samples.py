import gzip, json, sys
from pathlib import Path
import numpy as np
root = Path(sys.argv[1])
reference = None
seeds = []
for c in range(1,4):
    rep = root / "simulations" / f"replicate_{c:03d}"
    m = json.loads((rep / "metadata.json").read_text())
    ages = m["true_ancient_ages"]
    np.testing.assert_allclose(ages, np.arange(0,10001,500)/m["years_per_generation"])
    assert m["n_modern_haploid"] == 100 and m["n_ancient_haploid"] == 21
    assert m["sequence_length"] == 50000000
    assert not set(m["modern_sample_nodes"]) & set(m["ancient_sample_nodes"])
    with gzip.open(rep / "all_samples.vcf.gz", "rt") as h:
        names = next(line.rstrip().split("\t")[9:] for line in h if line.startswith("#CHROM"))
    assert names[100:] == [f"ancient_{i:02d}" for i in range(1,22)]
    signature = (names, ages, m["epochs"])
    if reference is None: reference = signature
    assert signature == reference, "chromosomes differ in sample identity, age, or demography"
    seeds.append(m["seeds"]["ancestry"])
assert len(set(seeds)) == 3
print("PASS: same 21 held-out individuals and ages across three independent chromosomes")
