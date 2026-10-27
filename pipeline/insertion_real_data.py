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

def parse_epsilon_grid(spec):
    """Parse START:STOP:STEP (inclusive) or a comma-separated epsilon grid."""
    if spec is None:
        return None
    if ":" in spec:
        fields=spec.split(":")
        if len(fields)!=3:
            raise ValueError("epsilon grid must be START:STOP:STEP or comma-separated")
        start,stop,step=map(float,fields)
        if step<=0 or stop<start:
            raise ValueError("epsilon grid requires STEP > 0 and STOP >= START")
        values=np.arange(start,stop+step/2,step,dtype=float)
    else:
        values=np.asarray([float(x) for x in spec.split(",") if x.strip()],dtype=float)
    if len(values)==0 or np.any(~np.isfinite(values)):
        raise ValueError("epsilon grid is empty or non-finite")
    if np.any((values<0)|(values>=.5)):
        raise ValueError("epsilon values must satisfy 0 <= epsilon < 0.5")
    if len(values)>1 and np.any(np.diff(values)<=0):
        raise ValueError("epsilon grid must be strictly increasing")
    return values

def call_loglik_grid(p, observed_alt, called, is_alt, epsilon_grid):
    """Per-epsilon call log likelihood, shaped (epsilon, sample, age)."""
    eps=np.asarray(epsilon_grid,dtype=float)[:,None]
    q=np.clip(eps+(1-2*eps)*np.asarray(p)[None,:],1e-300,1-1e-15)
    observed_alt=np.asarray(observed_alt); called=np.asarray(called)
    derived=observed_alt if is_alt else called-observed_alt
    return (derived[None,:,None]*np.log(q)[:,None,:] +
            (called-derived)[None,:,None]*np.log1p(-q)[:,None,:])

def marginalize_independent_chromosomes(draw_ll):
    """draw_ll has shape (draw, chromosome, sample, age)."""
    if draw_ll.ndim != 4 or not np.all(np.isfinite(draw_ll)):
        raise ValueError("draw likelihoods must be finite with four dimensions")
    return np.sum(logmeanexp(draw_ll, axis=0), axis=0)

def marginalize_epsilon_independent_chromosomes(draw_ll):
    """draw_ll has shape (draw, epsilon, chromosome, sample, age)."""
    if draw_ll.ndim!=5 or not np.all(np.isfinite(draw_ll)):
        raise ValueError("epsilon-grid draw likelihoods must be finite with five dimensions")
    return np.sum(logmeanexp(draw_ll,axis=0),axis=1)

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

def chromosome_layout(ts, chrom_fai=None):
    """Read concatenation offsets from TS metadata or reconstruct them from a FAI."""
    metadata = ts.metadata if isinstance(ts.metadata, dict) else {}
    records = metadata.get("chrom_offsets")
    if records:
        return {str(row["chrom"]): (float(row["offset"]), float(row["length"]))
                for row in records}
    if chrom_fai is None:
        raise ValueError("tree sequence lacks ARGtest chrom_offsets metadata; "
                         "provide --chrom-fai to reconstruct concatenation offsets")
    rows=[]
    with Path(chrom_fai).open() as handle:
        for line in handle:
            fields=line.split()
            if len(fields)>=2:
                rows.append((fields[0],float(fields[1])))
    # tskit concatenation commonly adds one coordinate between each input.
    # Require an exact prefix match to the TS length so the convention cannot
    # be selected silently or approximately.
    for padding in (0.0,1.0):
        offset=0.0; layout={}
        for chrom,length in rows:
            segment_length=length+padding
            layout[str(chrom)]=(offset,segment_length)
            offset += segment_length
            if offset==float(ts.sequence_length):
                return layout
            if offset>float(ts.sequence_length):
                break
    raise ValueError(f"FAI chromosome prefix does not reconstruct tree sequence length "
                     f"{ts.sequence_length:g}")

def call_files(path):
    found={}
    for p in path.glob("chr*.npz"):
        chrom=p.stem[3:]
        found[chrom]=p
    if not found:
        raise ValueError(f"no chr*.npz call files in {path}")
    return found

def save_draw_checkpoint(path, config, next_chrom_index, next_row, arrays, counts,
                         used, chromosomes):
    """Atomically save enough state to resume the target-site loop."""
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    payload={"config_json":np.asarray(json.dumps(config,sort_keys=True)),
             "next_chrom_index":np.asarray(next_chrom_index,dtype=np.int64),
             "next_row":np.asarray(next_row,dtype=np.int64),
             "counts_json":np.asarray(json.dumps(counts,sort_keys=True)),
             **arrays}
    for index,chrom in enumerate(chromosomes):
        payload[f"used_{index}"]=np.asarray(used.get(chrom,[]),dtype=np.int64)
    temporary=path.with_name(path.name+".tmp")
    with temporary.open("wb") as handle:
        np.savez(handle,**payload)
    temporary.replace(path)

def load_draw_checkpoint(path, expected_config, chromosomes):
    """Load and validate a checkpoint, rejecting incompatible run settings."""
    with np.load(path,allow_pickle=False) as saved:
        observed=json.loads(str(saved["config_json"].item()))
        if observed!=expected_config:
            raise ValueError(f"checkpoint configuration differs from this run: {path}")
        arrays={name:np.asarray(saved[name]) for name in
                ("ll","ll_epsilon","ll_by_block","sites_by_block",
                 "ll_by_d0","sites_by_d0") if name in saved}
        used={chrom:np.asarray(saved[f"used_{index}"],dtype=np.int64)
              for index,chrom in enumerate(chromosomes)}
        return (int(saved["next_chrom_index"]),int(saved["next_row"]),arrays,
                json.loads(str(saved["counts_json"].item())),used)

def run_draw(args):
    grid = np.arange(args.age_min, args.age_max + args.age_step/2, args.age_step)
    epsilon_grid = parse_epsilon_grid(args.epsilon_grid)
    ne = PiecewiseConstantNe.from_tsv(args.ne, args.ne_series, args.ne_extrapolate_last)
    calls_by_chrom = {}; samples = None
    for chrom, path in sorted(call_files(args.calls).items()):
        data = np.load(path, allow_pickle=False)
        names = list(data["samples"].astype(str))
        if samples is None: samples = names
        elif samples != names: raise SystemExit(f"sample order differs in chr{chrom}")
        calls_by_chrom[chrom] = data
    ts = tszip.decompress(args.tree)
    layout = chromosome_layout(ts,args.chrom_fai)
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
    ll_epsilon = (None if epsilon_grid is None else
                  np.zeros((len(epsilon_grid),len(chromosomes),len(samples),len(grid))))
    block_counts=[int(np.ceil(layout[c][1]/args.block_size)) for c in chromosomes]
    if len(set(block_counts)) != 1:
        raise ValueError("current block bootstrap requires equal block counts per chromosome")
    ll_by_block=np.zeros((len(chromosomes),block_counts[0],len(samples),len(grid)))
    sites_by_block=np.zeros((len(chromosomes),block_counts[0]),dtype=np.int64)
    ll_by_d0 = np.zeros((len(chromosomes), ts.num_samples + 1,
                         len(samples), len(grid)))
    sites_by_d0 = np.zeros((len(chromosomes), ts.num_samples + 1), dtype=np.int64)
    checkpoint_path=(args.checkpoint if args.checkpoint is not None else
                     Path(str(args.output)+".checkpoint.npz"))
    if args.checkpoint_every<=0:
        raise ValueError("--checkpoint-every must be positive")
    config={"format_version":1,"tree":str(args.tree.resolve()),"draw_id":args.draw_id,
            "calls":str(args.calls.resolve()),"chromosomes":chromosomes,"samples":samples,
            "targets":[len(calls_by_chrom[c]["position"]) for c in chromosomes],
            "grid":grid.tolist(),"epsilon":args.epsilon,
            "epsilon_grid":None if epsilon_grid is None else epsilon_grid.tolist(),
            "ne":str(args.ne.resolve()),"ne_series":args.ne_series,
            "ne_extrapolate_last":args.ne_extrapolate_last,
            "block_size":args.block_size}
    start_chrom_index=0; start_row=0
    if checkpoint_path.exists():
        start_chrom_index,start_row,arrays,counts,used=load_draw_checkpoint(
            checkpoint_path,config,chromosomes)
        ll=arrays["ll"]
        saved_epsilon=arrays.get("ll_epsilon")
        if (ll_epsilon is None)!=(saved_epsilon is None):
            raise ValueError("checkpoint epsilon-grid state differs from this run")
        ll_epsilon=saved_epsilon
        ll_by_block=arrays["ll_by_block"]; sites_by_block=arrays["sites_by_block"]
        ll_by_d0=arrays["ll_by_d0"]; sites_by_d0=arrays["sites_by_d0"]
        print(f"resuming checkpoint {checkpoint_path}: chromosome index "
              f"{start_chrom_index}, row {start_row}",flush=True)
    checkpoint_arrays=lambda: {
        "ll":ll,"ll_by_block":ll_by_block,"sites_by_block":sites_by_block,
        "ll_by_d0":ll_by_d0,"sites_by_d0":sites_by_d0,
        **({} if ll_epsilon is None else {"ll_epsilon":ll_epsilon})}
    tree = ts.first()
    for chrom_index, chrom in enumerate(chromosomes):
        if chrom_index<start_chrom_index:
            continue
        calls = calls_by_chrom[chrom]
        if chrom_index==start_chrom_index:
            row_start=start_row; kept=used.get(chrom,np.array([],dtype=np.int64)).tolist()
            stat=counts.get(str(chrom),{"targets":len(calls["position"]),"used":0,
                "absent":0,"multiple":0,"root":0,"allele_mismatch":0,
                "strand_complement":0,"invalid":0})
        else:
            row_start=0; kept=[]
            stat={"targets":len(calls["position"]),"used":0,"absent":0,"multiple":0,
                  "root":0,"allele_mismatch":0,"strand_complement":0,"invalid":0}
        local_positions=np.asarray(calls["position"],dtype=np.int64)
        if not 0<=row_start<=len(local_positions):
            raise ValueError(f"checkpoint row {row_start} outside chr{chrom} targets")
        for row in range(row_start,len(local_positions)):
            if row>row_start and row%args.checkpoint_every==0:
                used[chrom]=np.asarray(kept,dtype=np.int64); counts[str(chrom)]=stat
                save_draw_checkpoint(checkpoint_path,config,chrom_index,row,
                                     checkpoint_arrays(),counts,used,chromosomes)
                print(f"checkpoint chr{chrom}: {row}/{len(local_positions)} targets, "
                      f"{stat['used']} used",flush=True)
            local=local_positions[row]
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
            if ll_epsilon is not None:
                ll_epsilon[:,chrom_index] += call_loglik_grid(
                    p,observed,called,is_alt,epsilon_grid)
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
        save_draw_checkpoint(checkpoint_path,config,chrom_index+1,0,
                             checkpoint_arrays(),counts,used,chromosomes)
        print(f"chr{chrom}: {stat['used']}/{stat['targets']}", flush=True)
    args.output.mkdir(parents=True, exist_ok=False)
    np.save(args.output/"grid.npy",grid); np.save(args.output/"ll_by_chrom.npy",ll)
    np.save(args.output/"ll_by_chrom_d0.npy", ll_by_d0)
    np.save(args.output/"sites_by_chrom_d0.npy", sites_by_d0)
    np.save(args.output/"ll_by_block.npy",ll_by_block)
    np.save(args.output/"sites_by_block.npy",sites_by_block)
    if ll_epsilon is not None:
        np.save(args.output/"epsilon_grid.npy",epsilon_grid)
        np.save(args.output/"ll_by_chrom_epsilon.npy",ll_epsilon)
    (args.output/"samples.txt").write_text("\n".join(samples)+"\n")
    for chrom in chromosomes: np.save(args.output/f"used_chr{chrom}.npy",used[chrom])
    (args.output/"run.json").write_text(json.dumps({"model":"direct ancient-lineage insertion",
        "tree":str(args.tree.resolve()),"draw_id":args.draw_id,"epsilon":args.epsilon,
        "epsilon_grid":None if epsilon_grid is None else epsilon_grid.tolist(),
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
    arrays=[]; epsilon_paths=[]; arrays_d0=[]; site_counts_d0=[]; block_arrays=[]; block_counts=[]; reference=None; metas=[]
    epsilon_grid=None
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
        epsilon_path=part/"ll_by_chrom_epsilon.npy"
        if epsilon_path.exists():
            current_grid=np.load(part/"epsilon_grid.npy")
            if epsilon_grid is None: epsilon_grid=current_grid
            elif not np.array_equal(current_grid,epsilon_grid):
                raise SystemExit(f"epsilon grid differs in {part}")
            epsilon_paths.append(epsilon_path)
        if (part/"ll_by_block.npy").exists():
            block_arrays.append(np.load(part/"ll_by_block.npy")); block_counts.append(np.load(part/"sites_by_block.npy"))
        d0_path=part/"ll_by_chrom_d0.npy"
        if d0_path.exists() and not args.skip_d0:
            arrays_d0.append(np.load(d0_path))
            site_counts_d0.append(np.load(part/"sites_by_chrom_d0.npy"))
    draw_ll=np.stack(arrays); posterior_ll=marginalize_independent_chromosomes(draw_ll)
    args.output.mkdir(parents=True,exist_ok=False); np.save(args.output/"grid.npy",grid); np.save(args.output/"ll_marginal.npy",posterior_ll)
    if epsilon_paths:
        if len(epsilon_paths)!=len(arrays):
            raise SystemExit("only some draws contain epsilon-grid likelihoods")
        expected=(len(epsilon_grid),len(chromosomes),len(samples),len(grid))
        merged_epsilon=np.zeros((len(epsilon_grid),len(samples),len(grid)))
        # Stream one chromosome across draw files. This preserves independent
        # per-chromosome draw marginalisation without stacking a potentially
        # tens-of-gigabytes draw x epsilon x chromosome array in RAM.
        for chrom_index in range(len(chromosomes)):
            running_max=None; running_sum=None
            for epsilon_path in epsilon_paths:
                saved=np.load(epsilon_path,mmap_mode="r")
                if saved.shape!=expected:
                    raise SystemExit(f"epsilon likelihood shape in {epsilon_path} is "
                                     f"{saved.shape}; expected {expected}")
                current=np.asarray(saved[:,chrom_index],dtype=float)
                if running_max is None:
                    running_max=current.copy(); running_sum=np.ones_like(current)
                else:
                    new_max=np.maximum(running_max,current)
                    running_sum=(running_sum*np.exp(running_max-new_max)
                                 +np.exp(current-new_max))
                    running_max=new_max
            merged_epsilon += running_max+np.log(running_sum/len(epsilon_paths))
        np.save(args.output/"epsilon_grid.npy",epsilon_grid)
        np.save(args.output/"ll_marginal_epsilon.npy",merged_epsilon)
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
        "epsilon_grid":None if epsilon_grid is None else epsilon_grid.tolist(),
        "marginalisation":"logmeanexp draws within chromosome, then sum chromosomes",
        "ne":metas[0]["ne"],"chromosomes":metas[0]["chromosomes"],
        "sites_by_chromosome":{str(c):len(x) for c,x in zip(metas[0]["chromosomes"],reference)}},indent=2)+"\n")

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest="command",required=True)
    r=sub.add_parser("run-draw"); r.add_argument("--tree",type=Path,required=True); r.add_argument("--draw-id",type=int,required=True); r.add_argument("--calls",type=Path,required=True); r.add_argument("--chrom-fai",type=Path,help="FAI used to reconstruct offsets when the concatenated ARG lacks chrom_offsets metadata"); r.add_argument("--ne",type=Path,required=True); r.add_argument("--ne-series"); r.add_argument("--ne-extrapolate-last",action="store_true"); r.add_argument("--epsilon",type=float,default=.01); r.add_argument("--epsilon-grid",help="also evaluate START:STOP:STEP or comma-separated epsilon values"); r.add_argument("--age-min",type=float,default=0); r.add_argument("--age-max",type=float,default=30000); r.add_argument("--age-step",type=float,default=100); r.add_argument("--block-size",type=float,default=5_000_000); r.add_argument("--checkpoint",type=Path,help="resumable checkpoint NPZ [OUTPUT.checkpoint.npz]"); r.add_argument("--checkpoint-every",type=int,default=25000,help="checkpoint after this many target rows [25000]"); r.add_argument("--output",type=Path,required=True); r.set_defaults(func=run_draw)
    m=sub.add_parser("merge"); m.add_argument("--draw-results",type=Path,nargs="+"); m.add_argument("--manifest",type=Path); m.add_argument("--draw-root",type=Path); m.add_argument("--skip-d0",action="store_true",help="do not merge optional frequency-stratified likelihoods"); m.add_argument("--bootstrap-replicates",type=int,default=1000); m.add_argument("--bootstrap-seed",type=int,default=20260914); m.add_argument("--output",type=Path,required=True); m.set_defaults(func=merge)
    a=p.parse_args(); a.func(a)
if __name__=="__main__": main()
