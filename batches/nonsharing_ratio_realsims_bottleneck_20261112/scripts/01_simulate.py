#!/usr/bin/env python3
"""Simulate one replicate (BOTTLENECK demography from ../demography_bottleneck.yaml, haploid) and write panel / controls / ancient VCFs.

Conventions copied from validation/nonsharing_ratio_validation.py: Ne=20000 via
ploidy=1 + population_size=2*Ne (pair coalescence rate 1/(2Ne)), mu=1.25e-8, rho=1e-8,
20 Mb, BinaryMutationModel, discrete genome.  VCF POS = msprime position + 1
(1-based), CHROM = "1", REF = ancestral = "A", ALT = derived = "C".
Usage: 01_simulate.py REP(1..5) OUTDIR
"""
import gzip, json, sys
from pathlib import Path
import numpy as np, msprime, tskit

NE, MU, RHO, L = 20000, 1.25e-8, 1e-8, 20_000_000
N_PANEL, N_CTRL, N_ANC = 26, 10, 10
BASE_SEED = 20261112

rep = int(sys.argv[1]); out = Path(sys.argv[2]) / f"rep{rep}"; out.mkdir(parents=True, exist_ok=True)
seed = BASE_SEED + 100 * rep          # ancestry seed; mutations seed+1; ancient ages RNG seed+2
# stratified-uniform ancient ages: one uniform draw in each of 10 bins of width 600 over [0, 6000]
rng = np.random.default_rng(seed + 2)
ages = np.round(np.arange(N_ANC) * 600 + rng.uniform(0, 600, N_ANC), 2)

sets = [msprime.SampleSet(N_PANEL, time=0, ploidy=1), msprime.SampleSet(N_CTRL, time=0, ploidy=1)]
sets += [msprime.SampleSet(1, time=float(a), ploidy=1) for a in ages]
import demes
demog = msprime.Demography.from_demes(demes.load(str(Path(__file__).resolve().parent.parent / "demography_bottleneck.yaml")))
ts = msprime.sim_ancestry(samples=sets, ploidy=1, demography=demog, sequence_length=L,
                          recombination_rate=RHO, discrete_genome=True, random_seed=seed)
ts = msprime.sim_mutations(ts, rate=MU, model=msprime.BinaryMutationModel(),
                           discrete_genome=True, random_seed=seed + 1)
assert ts.num_samples == N_PANEL + N_CTRL + N_ANC
# sample order == order of SampleSets (verified below via node times)
st = ts.tables.nodes.time[ts.samples()]
assert np.all(st[:N_PANEL + N_CTRL] == 0) and np.allclose(st[N_PANEL + N_CTRL:], ages)

G = ts.genotype_matrix().astype(np.int8)          # (sites, 46)
pos1 = (ts.sites_position.astype(np.int64) + 1)
iP = np.arange(N_PANEL); iC = np.arange(N_PANEL, N_PANEL + N_CTRL); iA = np.arange(N_PANEL + N_CTRL, ts.num_samples)
names_P = [f"P{i+1:02d}" for i in range(N_PANEL)]
names_C = [f"C{i+1:02d}" for i in range(N_CTRL)]
names_A = [f"A{i+1:02d}" for i in range(N_ANC)]

def write_vcf(path, names, cols, sites):
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=01_simulate.py\n"
                 f"##contig=<ID=1,length={L}>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
        for s in sites:
            fh.write(f"1\t{pos1[s]}\t.\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(map(str, G[s, cols])) + "\n")

allsites = np.arange(ts.num_sites)
nP = G[:, iP].sum(1); nM = G[:, np.r_[iP, iC]].sum(1); nA = G[:, iA].sum(1)
seg = allsites[(nP > 0) & (nP < N_PANEL)]
keep_anc = allsites[nM > 0]
dropped = allsites[nM == 0]
write_vcf(out / "panel_full.vcf.gz", names_P, iP, allsites)           # pre-dropping, for prepare-fixed
write_vcf(out / "panel_segregating.vcf.gz", names_P, iP, seg)         # SINGER input (segregating in panel)
write_vcf(out / "controls.vcf.gz", names_C, iC, allsites)
write_vcf(out / "ancient.vcf.gz", names_A, iA, keep_anc)
(out / "panel_names.txt").write_text("\n".join(names_P) + "\n")
with open(out / "truth.tsv", "w") as fh:
    fh.write("rep\tsample\ttrue_age\n")
    for n, a in zip(names_A, ages): fh.write(f"{rep}\t{n}\t{a}\n")

# True panel ARG for sanity check (a): simplify to panel, shift to the VCF 1-based convention, A/C alleles.
sub = ts.simplify(list(iP))
t = sub.dump_tables()
t.sites.set_columns(position=t.sites.position + 1,
                    ancestral_state=np.frombuffer(b"A" * t.sites.num_rows, dtype=np.int8),
                    ancestral_state_offset=np.arange(t.sites.num_rows + 1, dtype=np.uint32),
                    metadata=t.sites.metadata, metadata_offset=t.sites.metadata_offset)
ds = np.array([b"C" if x == "1" else b"A" for x in
               tskit.unpack_strings(t.mutations.derived_state, t.mutations.derived_state_offset)])
t.mutations.set_columns(site=t.mutations.site, node=t.mutations.node, time=np.full(t.mutations.num_rows, tskit.UNKNOWN_TIME),
                        derived_state=np.frombuffer(b"".join(ds.tolist()), dtype=np.int8),
                        derived_state_offset=np.arange(len(ds) + 1, dtype=np.uint32),
                        parent=t.mutations.parent,
                        metadata=t.mutations.metadata, metadata_offset=t.mutations.metadata_offset)
t.sequence_length = L + 1
t.sort()
tsp = t.tree_sequence()
tsp.dump(out / "true_panel.trees")

info = dict(demography='bottleneck (0.3Ne at 500 gens ago, exponential recovery to Ne by 300)', rep=rep, seed_ancestry=seed, seed_mutations=seed + 1, seed_ages=seed + 2, ne=NE, mu=MU, rho=RHO, length=L,
            n_panel=N_PANEL, n_controls=N_CTRL, n_ancient=N_ANC, ancient_ages=ages.tolist(),
            n_sites_all=int(ts.num_sites), n_sites_panel_segregating=int(seg.size),
            n_sites_panel_fixed_derived=int(((nP == N_PANEL)).sum()),
            ancient_vcf_sites_kept=int(keep_anc.size),
            ancient_sites_removed_no_modern_derived=int(dropped.size),
            removed_sites_with_ancient_derived=int((nA[dropped] > 0).sum()),
            removed_sites_with_no_derived_anywhere=int((nA[dropped] == 0).sum()),
            true_panel_ts_sites=int(tsp.num_sites), true_panel_ts_trees=int(tsp.num_trees),
            mean_pairwise_div_panel=float(ts.simplify(list(iP)).diversity(mode="site")))
(out / "sim_info.json").write_text(json.dumps(info, indent=1))
print(json.dumps(info, indent=1))
