import tszip, numpy as np, sys
ts = tszip.decompress("/quobyte/jrigrp/beil/te_evo/singer_analysis/argtest/results/combined/tsz_files/run.combined.75.tsz")
print("seqlen", ts.sequence_length, "samples", ts.num_samples, "sites", ts.num_sites, "trees", ts.num_trees)
md = ts.metadata if isinstance(ts.metadata, dict) else {}
print("md keys", list(md.keys())[:10]); print((md.get("chrom_offsets") or [])[:12])
p = ts.tables.sites.position; print(p[:5], p[-3:])
print(ts.tables.sites[0].ancestral_state, [ts.site(i).ancestral_state for i in range(3)])
print("time units", ts.time_units, "max time", ts.max_root_time)
