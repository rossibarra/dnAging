#!/usr/bin/env python3
"""300 msprime simulations across randomised parameters, run in parallel.

Per simulation, drawn independently:
    Ne        uniform on [10,000, 100,000]          diploid effective size
    mu        LOG-uniform on [1e-9, 1e-8]           per base per generation
    rec       mu * one of {0.5, 1, 2}
    n_modern  uniform integer on [10, 40]           haploid chromosomes at time 0
    age       uniform on [100, 10,000] generations  one ancient haploid chromosome

Writes per simulation:
    sim_NNN.trees        ARG of the MODERN samples only; the ancient sample and its
                         sampling time are absent, and provenance is stripped so Ne,
                         the seeds and the sample times are not recoverable
    sim_NNN_modern.vcf   n_modern haploid columns
    sim_NNN_ancient.vcf  one haploid column, every biallelic site (no ascertainment)

and one parameters.tsv for the batch.
"""
import argparse, os, sys, time, json
import multiprocessing as mp
import numpy as np, msprime, tskit

SEQLEN = 10e6
MASTER = 20260907

def draw(i):
    rng = np.random.default_rng(MASTER + i)
    mu  = float(10 ** rng.uniform(-9, -8))
    return dict(
        sim       = f"sim_{i:03d}",
        seed      = int(rng.integers(1, 2**31 - 1)),
        Ne        = int(rng.integers(10_000, 100_001)),
        mu        = mu,
        rec_ratio = float(rng.choice([0.5, 1.0, 2.0])),
        n_modern  = int(rng.integers(10, 41)),
        age       = float(rng.uniform(100, 10_000)),
        seq_len   = SEQLEN,
    )

def one(p, outdir, tsv_only=False):
    t0 = time.time()
    rec = p["mu"] * p["rec_ratio"]
    ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(p["n_modern"], time=0,        ploidy=1),
                 msprime.SampleSet(1,             time=p["age"], ploidy=1)],
        population_size=p["Ne"], sequence_length=p["seq_len"],
        recombination_rate=rec, ploidy=2, random_seed=p["seed"])
    ts = msprime.sim_mutations(ts, rate=p["mu"], random_seed=p["seed"] + 1)

    modern  = [n for n in ts.samples() if ts.node(n).time == 0.0]
    ancient = [n for n in ts.samples() if ts.node(n).time == p["age"]][0]

    # Keep only sites with exactly ONE mutation.  That guarantees the allele list
    # is [ancestral, derived] so every VCF genotype is 0 or 1; filtering on the
    # number of observed allele indices is not enough, because a site can carry a
    # third allele in its list that no sample happens to show, which still emits
    # GT=2.  Applied to the FULL ts so both VCFs stay position-consistent.
    n_all = ts.num_sites
    multi = [s_.id for s_ in ts.sites() if len(s_.mutations) != 1]
    # a site at position 0 makes VCF POS=0, which tskit refuses to write
    multi += [s_.id for s_ in ts.sites() if s_.position < 1.0]
    ts = ts.delete_sites(sorted(set(multi)))

    mod_poly = sum(1 for v in ts.variants(samples=modern)
                   if np.unique(v.genotypes).size > 1)

    # ARG for the modern samples only, provenance and metadata cleared
    ts_mod = ts.simplify(samples=modern)
    tb = ts_mod.dump_tables()
    tb.provenances.clear()
    tb.metadata_schema = tskit.MetadataSchema.null()
    tb.metadata = b""
    if not tsv_only:
        tb.tree_sequence().dump(os.path.join(outdir, f"{p['sim']}.trees"))
        with open(os.path.join(outdir, f"{p['sim']}_modern.vcf"), "w") as f:
            ts_mod.write_vcf(f, contig_id="1")

    # Every biallelic site, ancient scored 0 or 1: no ascertainment applied.
    # Written from the FULL ts selecting the ancient individual -- simplifying to a
    # single sample leaves it isolated with no edges, which makes tskit report every
    # genotype as missing.
    anc_ind = ts.node(ancient).individual
    if not tsv_only:
        with open(os.path.join(outdir, f"{p['sim']}_ancient.vcf"), "w") as f:
            ts.write_vcf(f, contig_id="1", individuals=[anc_ind])

    anc_row = list(ts.samples()).index(ancient)
    anc_carried = sum(1 for v in ts.variants() if v.genotypes[anc_row] != 0)
    p = dict(p, rec=rec, n_trees=ts_mod.num_trees, n_sites_raw=n_all,
             n_sites_multiallelic=len(multi), n_sites=ts.num_sites,
             n_sites_modern_poly=mod_poly, ancient_carried=anc_carried,
             wall_sec=round(time.time() - t0, 1))
    print(f"  {p['sim']}  Ne={p['Ne']:>6,}  mu={p['mu']:.2e}  r/mu={p['rec_ratio']}  "
          f"n={p['n_modern']:>2}  age={p['age']:>7.0f}  sites={ts.num_sites:>7,}  "
          f"{p['wall_sec']:>6.1f}s", flush=True)
    return p

def worker(args):
    try:
        return one(*args)
    except Exception as e:                       # never let one run kill the batch
        print(f"  !! {args[0]['sim']} FAILED: {type(e).__name__}: {e}", flush=True)
        return dict(args[0], error=f"{type(e).__name__}: {e}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--procs", type=int, default=12)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--tsv-only", action="store_true",
                    help="recompute the stats and rewrite parameters.tsv without "
                         "touching the .trees/.vcf files")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    todo = [(draw(i), a.outdir, a.tsv_only)
            for i in range(a.start, a.start + a.n)]
    print(f"[batch] {len(todo)} simulations, {a.procs} processes, "
          f"{SEQLEN/1e6:g} Mb each -> {a.outdir}", flush=True)
    t0 = time.time()
    with mp.Pool(a.procs, maxtasksperchild=4) as pool:
        rows = pool.map(worker, todo, chunksize=1)
    cols = ["sim","seed","Ne","mu","rec","rec_ratio","n_modern","age","seq_len",
            "n_trees","n_sites_raw","n_sites_multiallelic","n_sites",
            "n_sites_modern_poly","ancient_carried","wall_sec","error"]
    # merge with anything already present, keyed on sim, so a partial rerun
    # updates rows rather than replacing the whole table
    path = os.path.join(a.outdir, "parameters.tsv")
    merged = {}
    if os.path.exists(path):
        with open(path) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                v = line.rstrip("\n").split("\t")
                merged[v[0]] = dict(zip(hdr, v))
    for r in rows:
        merged[r["sim"]] = {c: str(r.get(c, "")) for c in cols}
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for k in sorted(merged):
            f.write("\t".join(merged[k].get(c, "") for c in cols) + "\n")
    ok = sum(1 for r in rows if "error" not in r)
    print(f"[batch] {ok}/{len(rows)} succeeded in {(time.time()-t0)/60:.1f} min; "
          f"wrote parameters.tsv", flush=True)
