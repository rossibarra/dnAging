import json
import os
from pathlib import Path

import msprime
import numpy as np
import tskit

SOURCE = Path("/Users/jeffreyross-ibarra/Projects/mutrates/haploid_simulations_100mb_batch 2")
NE = int(os.environ.get("SIM_NE", "10000"))
LENGTH = float(os.environ.get("SIM_LENGTH", "100000000"))
OUT = Path(os.environ.get("SIM_BATCH_OUT", "/tmp/dnAging-ne10000-infinite-100mb"))
OUT.mkdir(parents=True, exist_ok=True)
ages = np.loadtxt(SOURCE / "ancient_sample_ages.tsv", skiprows=1, dtype=int)
params = json.loads((SOURCE / "parameters.json").read_text())
seed_by_sim = {int(x["simulation"]): (int(x["ancestry_seed"]), int(x["mutation_seed"]))
               for x in params["replicates"]}

records = []
selected = {
    int(x) for x in os.environ.get("SIM_IDS", "").split(",") if x.strip()
}
for sim, age in ages:
    if selected and int(sim) not in selected:
        continue
    name = f"simulation_{sim:02d}"
    target = OUT / name
    target.mkdir(exist_ok=True)
    ancestry_seed, mutation_seed = seed_by_sim[int(sim)]
    print(f"START {name} age={age} ancestry_seed={ancestry_seed} mutation_seed={mutation_seed}",
          flush=True)
    full = msprime.sim_ancestry(
        samples=[msprime.SampleSet(26, time=0, ploidy=1),
                 msprime.SampleSet(1, time=float(age), ploidy=1)],
        population_size=NE,
        sequence_length=LENGTH,
        recombination_rate=1e-8,
        ploidy=2,
        model="hudson",
        random_seed=ancestry_seed,
    )
    full = msprime.sim_mutations(
        full, rate=1e-8, model=msprime.InfiniteSites(msprime.NUCLEOTIDES),
        discrete_genome=False, random_seed=mutation_seed,
    )
    assert full.num_sites == full.num_mutations
    modern = [u for u in full.samples() if full.node(u).time == 0]
    ancient = [u for u in full.samples() if full.node(u).time == age]
    assert len(modern) == 26 and len(ancient) == 1
    modern_poly = []
    for var in full.variants(samples=modern):
        if 0 < np.count_nonzero(var.genotypes) < 26:
            modern_poly.append(var.site.id)
    discard = np.setdiff1d(np.arange(full.num_sites), np.asarray(modern_poly), assume_unique=True)
    asc = full.delete_sites(discard)
    assert asc.num_sites == asc.num_mutations
    modern_ts = asc.simplify(modern, filter_sites=False)
    ancient_ts = asc.simplify(ancient, filter_sites=False)
    tables = modern_ts.dump_tables()
    tables.provenances.clear()
    tables.metadata_schema = tskit.MetadataSchema.null()
    tables.metadata = b""
    modern_ts = tables.tree_sequence()
    modern_ts.dump(target / f"{name}.trees")
    with (target / f"{name}_modern.vcf").open("w") as fh:
        modern_ts.write_vcf(fh, contig_id="1", individual_names=[f"modern_{i:02d}" for i in range(1, 27)])
    with (target / f"{name}_ancient.vcf").open("w") as fh:
        ancient_ts.write_vcf(fh, contig_id="1", individual_names=["ancient_01"])
    record = {
        "simulation": int(sim), "true_age": int(age),
        "ancestry_seed": ancestry_seed, "mutation_seed": mutation_seed,
        "modern_sites": int(modern_ts.num_sites), "trees": int(modern_ts.num_trees),
    }
    records.append(record)
    print(f"DONE {name}: {modern_ts.num_sites} sites, {modern_ts.num_trees} trees", flush=True)

(OUT / "run.json").write_text(json.dumps({
    "Ne": NE, "mutation_rate": 1e-8, "recombination_rate": 1e-8,
    "sequence_length": LENGTH, "modern_haplotypes": 26,
    "ancient_haplotypes": 1, "mutation_model": "InfiniteSites(NUCLEOTIDES)",
    "discrete_genome": False, "replicates": records,
}, indent=2) + "\n")
