#!/usr/bin/env python3
"""Simulate a 100 Mb region with 26 modern haploid chromosomes plus one ancient
individual, and write:

    <prefix>.trees          ARG for the 26 MODERN samples only, provenance stripped
    <prefix>_modern.vcf     26 haploid columns
    <prefix>_ancient.vcf    ONE chromosome of the ancient individual, haploid

The ancient sample takes part in the coalescent, so it shares the genealogy and the
mutations, but it is absent from the written tree sequence and its sampling time is
not recoverable from it.

Usage:  simulate_ancient.py --age X --seed S [--prefix simulation_01]
"""
import argparse, sys
import msprime, tskit, numpy as np

p = argparse.ArgumentParser()
p.add_argument("--age",    type=float, required=True, help="ancient sample age, generations")
p.add_argument("--seed",   type=int,   required=True)
p.add_argument("--prefix", default="simulation_01")
p.add_argument("--ne",     type=float, default=25_000, help="DIPLOID effective size")
p.add_argument("--length", type=float, default=100e6)
p.add_argument("--rec",    type=float, default=1e-8)
p.add_argument("--mu",     type=float, default=1e-8)
p.add_argument("--nmod",   type=int,   default=26, help="modern haploid chromosomes")
p.add_argument("--ancient-carried-only", action="store_true",
               help="emit only sites where the ancient carries the derived allele "
                    "(default: emit every ascertained site as 0 or 1)")
a = p.parse_args()

# ---- 1. ancestry -----------------------------------------------------------
# ploidy=2 at the simulation level sets the coalescent timescale to 2*Ne
# generations, which is what makes --ne a DIPLOID effective size.  The ancient is
# simulated as a diploid individual so that "one chromosome" is meaningful.
ts = msprime.sim_ancestry(
    samples=[msprime.SampleSet(a.nmod, time=0,     ploidy=1),
             msprime.SampleSet(1,      time=a.age, ploidy=2)],
    population_size=a.ne,
    sequence_length=a.length,
    recombination_rate=a.rec,
    ploidy=2,
    random_seed=a.seed,
)

# ---- 2. mutations on the FULL tree sequence --------------------------------
# Must come before any simplify, or the ancient would not share the mutations.
ts = msprime.sim_mutations(ts, rate=a.mu, random_seed=a.seed + 1)

modern  = [n for n in ts.samples() if ts.node(n).time == 0.0]
ancient = [n for n in ts.samples() if ts.node(n).time == a.age]
anc_chrom = ancient[0]                      # ONE chromosome of the ancient individual

# ---- 3. drop multiallelic sites --------------------------------------------
# discrete_genome=True (the default) lets two mutations hit one position.  Drop
# those on the FULL ts so the modern and ancient outputs stay position-consistent.
multi = [v.site.id for v in ts.variants() if np.unique(v.genotypes).size > 2]
ts = ts.delete_sites(multi)
print(f"[sim] dropped {len(multi):,} multiallelic sites; {ts.num_sites:,} remain",
      file=sys.stderr)

# ---- 4. the ascertainment set: polymorphic among the 26 modern -------------
# Real ancient samples are genotyped only at SNPs ascertained in a modern panel,
# so mutations private to the ancient lineage must not appear in its VCF.
mod_poly = {v.site.id for v in ts.variants(samples=modern)
            if np.unique(v.genotypes).size > 1}
print(f"[sim] {len(mod_poly):,} sites polymorphic in the modern panel", file=sys.stderr)

# ---- 5. ARG for the modern samples, with the ancient and its time removed ---
ts_mod = ts.simplify(samples=modern)
tables = ts_mod.dump_tables()
tables.provenances.clear()        # provenance records Ne, seeds AND the sample times
tables.metadata_schema = tskit.MetadataSchema.null()
tables.metadata = b""
tables.tree_sequence().dump(f"{a.prefix}.trees")

# ---- 6. VCFs ---------------------------------------------------------------
with open(f"{a.prefix}_modern.vcf", "w") as f:
    ts_mod.write_vcf(f, contig_id="1")

asc_only = ts.delete_sites([s.id for s in ts.sites() if s.id not in mod_poly])
ts_anc = asc_only.simplify(samples=[anc_chrom],
                           filter_sites=a.ancient_carried_only)
with open(f"{a.prefix}_ancient.vcf", "w") as f:
    ts_anc.write_vcf(f, contig_id="1")

print(f"[sim] wrote {a.prefix}.trees ({ts_mod.num_samples} samples, "
      f"{ts_mod.num_trees:,} trees, {ts_mod.num_sites:,} sites)", file=sys.stderr)
print(f"[sim] wrote {a.prefix}_modern.vcf and {a.prefix}_ancient.vcf "
      f"({ts_anc.num_sites:,} sites)", file=sys.stderr)
