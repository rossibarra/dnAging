import sys, numpy as np, tszip
sys.path.insert(0,"pipeline")
import nonsharing_ratio as nr
ts=tszip.decompress("batches/nonsharing_ratio_heldout_20261111/singer_out/chr.10.combined.snp.te.sorted/trees/chr.10.combined.snp.te.sorted.75.tsz")
mid=nr.edge_midpoint_ages(ts)
# child-node time per site (oldest mutation), same convention
t=ts.tables; times=np.asarray(t.nodes.time)
ch=np.full(ts.num_sites,np.nan)
for m_node,m_site in zip(np.asarray(t.mutations.node),np.asarray(t.mutations.site)):
    v=times[m_node]
    if np.isnan(ch[m_site]) or v>ch[m_site]: ch[m_site]=v
ok=~np.isnan(mid); print("sites with mutation",ok.sum())
for A0 in (12000,30000,100000,300000,1000000):
    print(A0,"midpoint>A0: %.3f  child>A0: %.3f"%((mid[ok]>A0).mean(),(ch[ok]>A0).mean()))
print("child-time quantiles 10/50/90/99:",np.percentile(ch[ok],[10,50,90,99]).round(0))
print("midpoint quantiles 10/50/90/99:",np.percentile(mid[ok][np.isfinite(mid[ok])],[10,50,90,99]).round(0))
