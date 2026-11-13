"""Ancestral alleles for chr10 from a real whole-genome draw built on all 26 haplotypes (before any
held-out lines were dropped).  Writes 'chrom pos ancestral' for single-base ancestral states.
Position convention: VCF pos = ARG site position - chromosome offset (shift 0, verified)."""
import sys
import numpy as np
sys.path.insert(0, "pipeline")
import nonsharing_ratio_real_data as rd
tree, fai, out = sys.argv[1:4]
ts = rd.load_tree(tree)
layout = rd.layout_for(ts, fai)
print({k: v for k, v in layout.items() if k in ("10", "chr10")}, flush=True)
key = "10" if "10" in layout else "chr10"
off, length = layout[key][0], layout[key][1]
t = ts.tables.sites
pos = np.asarray(t.position)
sel = np.nonzero((pos >= off) & (pos < off + length + 1))[0]
st = np.asarray(t.ancestral_state); so = np.asarray(t.ancestral_state_offset)
n1 = (so[sel + 1] - so[sel]) == 1
sel = sel[n1]
anc = st[so[sel]].tobytes().decode() if False else bytes(st[so[sel]].astype(np.uint8)).decode()
with open(out, "w") as fh:
    for i, a in zip(sel, anc):
        fh.write(f"10\t{int(round(pos[i] - off))}\t{a}\n")
print("sites in chr10 range:", int(((pos >= off) & (pos < off + length + 1)).sum()), "written (single-base):", sel.size)
