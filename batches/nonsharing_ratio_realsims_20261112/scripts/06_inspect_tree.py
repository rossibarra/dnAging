import sys, tszip, numpy as np
ts = tszip.decompress(sys.argv[1]); nt = ts.tables.nodes.time
print("samples", ts.num_samples, "seqlen", ts.sequence_length, "sites", ts.num_sites, "trees", ts.num_trees, "max node time", nt.max())
roots = [t.time(t.root) for t in ts.trees() if t.num_roots == 1]
print("root time mean over trees (n=%d): %.0f  median %.0f (expected panel TMRCA ~ 2*N*(1-1/26) = %.0f)" % (len(roots[::50]), np.mean(roots[::50]), np.median(roots[::50]), 2*40000*(1-1/26)/1))
print("provenance:", [p.record[:300] for p in ts.provenances()][-1:])
print("pairwise diversity:", ts.diversity(mode="site"), "pi/(4mu)=Ne:", ts.diversity(mode="site")/(4*1.25e-8))
print("metadata keys", list(ts.metadata) if isinstance(ts.metadata, dict) else ts.metadata)
