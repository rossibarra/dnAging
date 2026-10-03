#!/usr/bin/env python3
"""Ancient-lineage insertion on concatenated, chromosome-independent ARG draws."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import tszip
from betabinom_real_data import orientation, summarize
from insertion_likelihood import PiecewiseConstantNe, derived_probability_uniform_edge_grid

def logmeanexp(values, axis=0):
    peak = np.max(values, axis=axis, keepdims=True)
    return np.squeeze(peak, axis=axis) + np.log(np.mean(np.exp(values-peak), axis=axis))

def call_loglik(p, observed_alt, called, is_alt, epsilon):
    """Log likelihood when p is the probability of the ARG-derived state."""
    q=np.clip(epsilon+(1-2*epsilon)*p,1e-300,1-1e-15)
    observed_alt=np.asarray(observed_alt); called=np.asarray(called)
    derived=observed_alt if is_alt else called-observed_alt
    return (derived[:,None]*np.log(q)[None,:] +
            (called-derived)[:,None]*np.log1p(-q)[None,:])

def marginalize_independent_chromosomes(draw_ll):
    """draw_ll has shape (draw, chromosome, sample, age)."""
    if draw_ll.ndim != 4 or not np.all(np.isfinite(draw_ll)):
        raise ValueError("draw likelihoods must be finite with four dimensions")
    return np.sum(logmeanexp(draw_ll, axis=0), axis=0)

def block_bootstrap_maps(draw_block_ll, grid, n_bootstrap=1000, seed=20260914):
    """Bootstrap blocks within chromosomes and return MAPs (bootstrap, sample)."""
    if draw_block_ll.ndim != 5:
        raise ValueError("block likelihoods require draw, chromosome, block, sample, age")
    rng=np.random.default_rng(seed); _,nchrom,nblock,_,_=draw_block_ll.shape
    maps=[]
    # Work in modest batches so 50 draws x 3 chromosomes remains memory-safe.
    for start in range(0,n_bootstrap,50):
        k=min(50,n_bootstrap-start)
        weights=np.zeros((k,nchrom,nblock),dtype=np.int16)
        for boot in range(k):
            for chrom in range(nchrom):
                weights[boot,chrom]=np.bincount(rng.integers(0,nblock,nblock),minlength=nblock)
        sampled=np.einsum("kcb,dcbst->kdcst",weights,draw_block_ll,optimize=True)
        peak=np.max(sampled,axis=1,keepdims=True)
        ll=np.sum(np.squeeze(peak,axis=1)+np.log(np.mean(np.exp(sampled-peak),axis=1)),axis=1)
        maps.append(grid[np.argmax(ll,axis=-1)])
    return np.concatenate(maps,axis=0)

def chromosome_layout(ts):
    """Read ARGtest's authoritative concatenation offsets from TS metadata."""
    metadata = ts.metadata if isinstance(ts.metadata, dict) else {}
    records = metadata.get("chrom_offsets")
    if not records:
        raise ValueError("tree sequence lacks ARGtest chrom_offsets metadata")
    layout = {str(row["chrom"]): (float(row["offset"]), float(row["length"]))
              for row in records}
    return layout

def call_files(path):
    found={}
    for p in path.glob("chr*.npz"):
        chrom=p.stem[3:]
        found[chrom]=p
    if not found:
        raise ValueError(f"no chr*.npz call files in {path}")
    return found

def run_draw(args):
    grid = np.arange(args.age_min, args.age_max + args.age_step/2, args.age_step)
    ne = PiecewiseConstantNe.from_tsv(args.ne, args.ne_series, args.ne_extrapolate_last)
    calls_by_chrom = {}; samples = None
    for chrom, path in sorted(call_files(args.calls).items()):
        data = np.load(path, allow_pickle=False)
        names = list(data["samples"].astype(str))
        if samples is None: samples = names
        elif samples != names: raise SystemExit(f"sample order differs in chr{chrom}")
        calls_by_chrom[chrom] = data
    ts = tszip.decompress(args.tree)
    layout = chromosome_layout(ts)
    # Accept the two conventions produced by ARGtest configurations:
    # chrom_offsets may say either "1" or "chr1", while chr1.npz naturally
    # yields the suffix "1".
    remapped={}
    for chrom,data in calls_by_chrom.items():
        key=chrom if chrom in layout else f"chr{chrom}"
        if key not in layout:
            raise ValueError(f"chrom_offsets metadata lacks chromosome {chrom}")
        if key in remapped:
            raise ValueError(f"duplicate call files resolve to chromosome {key}")
        remapped[key]=data
    calls_by_chrom=remapped
    positions = np.asarray(ts.tables.sites.position)
    chromosomes=sorted(calls_by_chrom)
    missing=[c for c in chromosomes if c not in layout]
    if missing: raise ValueError(f"chrom_offsets metadata lacks chromosomes {missing}")
    ll = np.zeros((len(chromosomes), len(samples), len(grid))); used = {}; counts = {}
    block_counts=[int(np.ceil(layout[c][1]/args.block_size)) for c in chromosomes]
    if len(set(block_counts)) != 1:
        raise ValueError("current block bootstrap requires equal block counts per chromosome")
    ll_by_block=np.zeros((len(chromosomes),block_counts[0],len(samples),len(grid)))
    sites_by_block=np.zeros((len(chromosomes),block_counts[0]),dtype=np.int64)
    ll_by_d0 = np.zeros((len(chromosomes), ts.num_samples + 1,
                         len(samples), len(grid)))
    sites_by_d0 = np.zeros((len(chromosomes), ts.num_samples + 1), dtype=np.int64)
    tree = ts.first()
    for chrom_index, chrom in enumerate(chromosomes):
        calls = calls_by_chrom[chrom]; kept=[]
        stat={"targets":len(calls["position"]),"used":0,"absent":0,"multiple":0,
              "root":0,"allele_mismatch":0,"strand_complement":0,"invalid":0}
        for row, local in enumerate(np.asarray(calls["position"], dtype=np.int64)):
            offset, length = layout[chrom]
            if not 0 <= local < length:
                raise ValueError(f"chr{chrom} position {local} outside length {length:g}")
            global_pos = float(local) + offset
            site_id = int(np.searchsorted(positions, global_pos))
            if site_id >= len(positions) or positions[site_id] != global_pos:
                stat["absent"] += 1; continue
            site=ts.site(site_id)
            if len(site.mutations) != 1: stat["multiple"] += 1; continue
            mutation=site.mutations[0]; tree.seek(global_pos); parent=tree.parent(mutation.node)
            if parent == -1: stat["root"] += 1; continue
            is_alt, match = orientation(site, mutation, str(calls["ref"][row]), str(calls["alt"][row]))
            if is_alt is None: stat["allele_mismatch"] += 1; continue
            if match == "strand_complement": stat["strand_complement"] += 1
            try:
                p = derived_probability_uniform_edge_grid(
                    tree, mutation.node, grid, tree.time(mutation.node), tree.time(parent), ne)
            except ValueError:
                stat["invalid"] += 1; continue
            observed=np.asarray(calls["observed_alt"][row]); called=np.asarray(calls["called"][row])
            contribution=call_loglik(p,observed,called,is_alt,args.epsilon)
            ll[chrom_index] += contribution
            block=min(int(local//args.block_size),block_counts[chrom_index]-1)
            ll_by_block[chrom_index,block] += contribution
            sites_by_block[chrom_index,block] += 1
            if "panel_alt_count" in calls and "panel_called" in calls:
                # Frequency strata use the observed, fixed VCF ALT count.  In
                # these simulations REF is known ancestral; inferred ARG draws
                # may disagree about polarity, but must not move a SNP between
                # frequency bins as a consequence.
                d0 = int(calls["panel_alt_count"][row])
            else:
                d0 = int(tree.num_samples(mutation.node))
            ll_by_d0[chrom_index, d0] += contribution
            sites_by_d0[chrom_index, d0] += 1
            kept.append(int(local)); stat["used"] += 1
        used[chrom]=np.asarray(kept,dtype=np.int64); counts[str(chrom)]=stat
        print(f"chr{chrom}: {stat['used']}/{stat['targets']}", flush=True)
    args.output.mkdir(parents=True, exist_ok=False)
    np.save(args.output/"grid.npy",grid); np.save(args.output/"ll_by_chrom.npy",ll)
    np.save(args.output/"ll_by_chrom_d0.npy", ll_by_d0)
    np.save(args.output/"sites_by_chrom_d0.npy", sites_by_d0)
    np.save(args.output/"ll_by_block.npy",ll_by_block)
    np.save(args.output/"sites_by_block.npy",sites_by_block)
    (args.output/"samples.txt").write_text("\n".join(samples)+"\n")
    for chrom in chromosomes: np.save(args.output/f"used_chr{chrom}.npy",used[chrom])
    (args.output/"run.json").write_text(json.dumps({"model":"direct ancient-lineage insertion",
        "tree":str(args.tree.resolve()),"draw_id":args.draw_id,"epsilon":args.epsilon,
        "chromosomes":chromosomes,"ne":str(args.ne.resolve()),"block_size":args.block_size,
        "counts_by_chromosome":counts},indent=2)+"\n")

def merge(args):
    parts=list(args.draw_results or [])
    if args.manifest:
        with args.manifest.open() as handle:
            header=handle.readline().rstrip("\n").split("\t")
            draw_col=header.index("draw_id")
            ids=[line.rstrip("\n").split("\t")[draw_col] for line in handle if line.strip()]
        parts=[args.draw_root/f"draw_{draw_id}" for draw_id in ids]
    samples=(parts[0]/"samples.txt").read_text().split(); grid=np.load(parts[0]/"grid.npy")
    arrays=[]; arrays_d0=[]; site_counts_d0=[]; block_arrays=[]; block_counts=[]; reference=None; metas=[]
    for part in parts:
        if (part/"samples.txt").read_text().split()!=samples: raise SystemExit(f"sample order differs in {part}")
        if not np.array_equal(np.load(part/"grid.npy"),grid): raise SystemExit(f"grid differs in {part}")
        meta=json.loads((part/"run.json").read_text())
        chromosomes=meta["chromosomes"]
        if metas and chromosomes != metas[0]["chromosomes"]: raise SystemExit(f"chromosomes differ in {part}")
        current=[np.load(part/f"used_chr{c}.npy") for c in chromosomes]
        if reference is None: reference=current
        elif any(not np.array_equal(x,y) for x,y in zip(reference,current)):
            raise SystemExit(f"retained sites differ in {part}; provide a common-draw site mask")
        arrays.append(np.load(part/"ll_by_chrom.npy")); metas.append(meta)
        if (part/"ll_by_block.npy").exists():
            block_arrays.append(np.load(part/"ll_by_block.npy")); block_counts.append(np.load(part/"sites_by_block.npy"))
        d0_path=part/"ll_by_chrom_d0.npy"
        if d0_path.exists() and not args.skip_d0:
            arrays_d0.append(np.load(d0_path))
            site_counts_d0.append(np.load(part/"sites_by_chrom_d0.npy"))
    draw_ll=np.stack(arrays); posterior_ll=marginalize_independent_chromosomes(draw_ll)
    args.output.mkdir(parents=True,exist_ok=False); np.save(args.output/"grid.npy",grid); np.save(args.output/"ll_marginal.npy",posterior_ll)
    if arrays_d0:
        if len(arrays_d0) != len(arrays): raise SystemExit("only some draws contain d0-stratified likelihoods")
        by_d0=np.stack(arrays_d0)  # draw, chromosome, d0, sample, age
        np.save(args.output/"ll_marginal_by_d0.npy",
                np.sum(logmeanexp(by_d0,axis=0),axis=0))
        if any(not np.array_equal(site_counts_d0[0],x) for x in site_counts_d0[1:]):
            raise SystemExit("derived-count site allocation differs across draws")
        np.save(args.output/"sites_by_d0.npy",np.sum(site_counts_d0[0],axis=0))
    if block_arrays:
        if len(block_arrays)!=len(arrays): raise SystemExit("only some draws contain block likelihoods")
        if any(not np.array_equal(block_counts[0],x) for x in block_counts[1:]):
            raise SystemExit("block site allocation differs across draws")
        draw_blocks=np.stack(block_arrays)
        np.save(args.output/"ll_by_block_and_draw.npy",draw_blocks)
        np.save(args.output/"sites_by_block.npy",block_counts[0])
        boot=block_bootstrap_maps(draw_blocks,grid,args.bootstrap_replicates,args.bootstrap_seed)
        np.save(args.output/"block_bootstrap_maps.npy",boot)
        lo,hi=np.quantile(boot,[.025,.975],axis=0)
        with (args.output/"block_bootstrap_ci.tsv").open("w") as handle:
            handle.write("sample\tbootstrap_map_ci95_lower\tbootstrap_map_ci95_upper\n")
            for sample,a,b in zip(samples,lo,hi): handle.write(f"{sample}\t{a:.6g}\t{b:.6g}\n")
    (args.output/"samples.txt").write_text("\n".join(samples)+"\n")
    with (args.output/"ages_table.tsv").open("w") as handle:
        handle.write("sample\tmap_T\tmean_T\tmedian_T\tci95_lower_T\tci95_upper_T\n")
        for sample, lp in zip(samples,posterior_ll):
            handle.write(sample+"\t"+"\t".join(f"{x:.6g}" for x in summarize(grid,lp))+"\n")
    (args.output/"run.json").write_text(json.dumps({"model":"direct ancient-lineage insertion",
        "n_arg_draws":len(parts),"draw_ids":[m["draw_id"] for m in metas],"epsilon":metas[0]["epsilon"],
        "marginalisation":"logmeanexp draws within chromosome, then sum chromosomes",
        "ne":metas[0]["ne"],"chromosomes":metas[0]["chromosomes"],
        "sites_by_chromosome":{str(c):len(x) for c,x in zip(metas[0]["chromosomes"],reference)}},indent=2)+"\n")

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest="command",required=True)
    r=sub.add_parser("run-draw"); r.add_argument("--tree",type=Path,required=True); r.add_argument("--draw-id",type=int,required=True); r.add_argument("--calls",type=Path,required=True); r.add_argument("--ne",type=Path,required=True); r.add_argument("--ne-series"); r.add_argument("--ne-extrapolate-last",action="store_true"); r.add_argument("--epsilon",type=float,default=.01); r.add_argument("--age-min",type=float,default=0); r.add_argument("--age-max",type=float,default=30000); r.add_argument("--age-step",type=float,default=100); r.add_argument("--block-size",type=float,default=5_000_000); r.add_argument("--output",type=Path,required=True); r.set_defaults(func=run_draw)
    m=sub.add_parser("merge"); m.add_argument("--draw-results",type=Path,nargs="+"); m.add_argument("--manifest",type=Path); m.add_argument("--draw-root",type=Path); m.add_argument("--skip-d0",action="store_true",help="do not merge optional frequency-stratified likelihoods"); m.add_argument("--bootstrap-replicates",type=int,default=1000); m.add_argument("--bootstrap-seed",type=int,default=20260914); m.add_argument("--output",type=Path,required=True); m.set_defaults(func=merge)
    a=p.parse_args(); a.func(a)
if __name__=="__main__": main()
