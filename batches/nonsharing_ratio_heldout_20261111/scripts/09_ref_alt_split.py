"""Does the ancient-vs-control non-sharing offset depend on whether the derived allele is the VCF REF or ALT?
Same site selection as run-draw (child-node ages, A0 given); prints aggregates only (no sample-level output).
usage: 09_ref_alt_split.py TREE CONTROLS_VCF ANCIENT_VCF PANEL_NAMES CHROM_FAI A0"""
import sys
import numpy as np
sys.path.insert(0, "pipeline")
import nonsharing_ratio as nr
import nonsharing_ratio_io as nio
import nonsharing_ratio_real_data as rd
tree, cvcf, avcf, pnames, fai, a0 = sys.argv[1:7]; a0 = float(a0)
ts = rd.load_tree(tree); layout = rd.layout_for(ts, fai)
controls = nio.read_vcf(cvcf, mode="sample", chroms=rd.vcf_names(layout))
ancient = nio.read_vcf(avcf, mode="sample", chroms=rd.vcf_names(layout))
positions = np.asarray(ts.tables.sites.position)
anc_state, der_state = nio.arg_site_alleles(ts)
panel = nr.panel_derived_matrix(ts); ages = nr.edge_child_ages(ts)
ctrl_d, _ = nio.polarise(controls, rd.match_vcf_chromosomes(layout, controls), positions, anc_state, der_state, 0.0)
anc_d, _ = nio.polarise(ancient, rd.match_vcf_chromosomes(layout, ancient), positions, anc_state, der_state, 0.0)
# VCF REF of the controls at each ARG site: class 'alt_derived' = ARG ancestral state equals REF
ref_at = np.full(positions.size, "", dtype=object)
for chrom, (offset, length) in rd.match_vcf_chromosomes(layout, controls).items() if False else rd.match_vcf_chromosomes(layout, controls).items():
    pass
m = rd.match_vcf_chromosomes(layout, controls)
for chrom, (offset, length) in m.items():
    rec = controls.chroms[chrom]
    tgt = rec["pos"].astype(float) + offset
    o = np.argsort(tgt); tgt, refs = tgt[o], np.asarray(rec["ref"])[o]
    i = np.minimum(np.searchsorted(tgt, positions), tgt.size - 1)
    hit = tgt[i] == positions
    ref_at[hit] = refs[i[hit]]
usable = (panel >= 0).all(axis=0) & (der_state != "") & (ref_at != "")
targets = np.vstack([anc_d, ctrl_d])[:, usable]
called = targets >= 0; targets = np.where(called, targets, 0)
panel, ages = panel[:, usable], ages[usable]
cls_alt_derived = (anc_state[usable] == ref_at[usable])      # ancestral == REF  (derived = ALT)
cls_ref_derived = (der_state[usable] == ref_at[usable])      # derived == REF
mask = nr.select_sites(panel, targets, ages, a0, include_panel_fixed=True, called=called)
n_anc = anc_d.shape[0]
print(f"A0={a0:.0f} child ages; usable {int(usable.sum())}, selected {int(mask.sum())}")
print(f"selected: derived=ALT (anc==REF) {int((mask & cls_alt_derived).sum())}, derived=REF {int((mask & cls_ref_derived).sum())}, neither {int((mask & ~cls_alt_derived & ~cls_ref_derived).sum())}")
res = {}
for name, cls in (("derived=ALT", cls_alt_derived), ("derived=REF", cls_ref_derived), ("all", np.ones_like(cls_alt_derived))):
    mk = mask & cls
    num, den = nr.nonsharing_counts(panel, targets, mk, called)
    f = num[0] / den[0]
    fa, fm = f[:n_anc], f[n_anc:]
    res[name] = (np.median(fa), fm.mean())
    print(f"{name:12s}: controls f {fm.mean():.4f}  ancient median f {np.median(fa):.4f} (IQR {np.percentile(fa,25):.4f}-{np.percentile(fa,75):.4f})  ratio {np.median(fa)/fm.mean():.3f}  H {np.log(np.median(fa)/fm.mean()):.3f}")
