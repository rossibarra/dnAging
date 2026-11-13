"""Summarise one stand-in draw + its run-draw NPZ, and rewrite the NPZ with de-identified ancient names."""
import json, sys, numpy as np, tszip
b = "/quobyte/jrigrp/beil/logan_collab/dnAging/batches/nonsharing_ratio_heldout_20261111"
draw = sys.argv[1]
ts = tszip.decompress(f"{b}/singer_out/chr.10.combined.snp.te.sorted/trees/chr.10.combined.snp.te.sorted.{draw}.tsz")
pos = ts.tables.sites.position
print("ts: seqlen", ts.sequence_length, "samples", ts.num_samples, "individuals", ts.num_individuals,
      "sites", ts.num_sites, "mutations", ts.num_mutations, "trees", ts.num_trees,
      "time_units", ts.time_units, "root_time_max", float(ts.max_root_time))
print("site pos min/max", pos.min(), pos.max(), "metadata keys", list(ts.metadata.keys()) if isinstance(ts.metadata, dict) else None)
print("n sites with ancestral_state len>1:", sum(len(s.ancestral_state) > 1 for s in ts.sites()))
d = np.load(f"{b}/verify/run_draw_{draw}.npz", allow_pickle=False)
meta = json.loads(str(d["meta"]))
for k, v in meta.items():
    if k not in ("tree",): print(k, v)
tg = [str(x) for x in d["targets"]]; ic = d["is_control"]
print("targets", len(tg), "controls", [t for t, c in zip(tg, ic) if c])
T = d["T_hat"]; print("T_hat ancients n/min/median/max", len(T), np.nanmin(T), np.nanmedian(T), np.nanmax(T), "n_ctrl_T", int(ic.sum()))
print("num/den totals", float(d["num"].sum()), float(d["den"].sum()), "blocks", d["block_labels"].shape)
anon = [f"anc{i+1:03d}" if not c else t for i, (t, c) in enumerate(zip(tg, ic))]
out = {k: d[k] for k in d.files}; out["targets"] = np.array(anon)
np.savez(f"{b}/verify/run_draw_{draw}.deidentified.npz", **out)
