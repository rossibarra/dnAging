#!/usr/bin/env python3
"""Choose k held-out modern lines reproducibly (stdlib only).

Candidates = every panel column except 'REF' (the reference-genome column that
argprep adds; SINGER needs it in the panel).  random.Random(SEED).sample over the
sorted candidate list.  Writes inputs/{heldout,panel}_samples.txt and selection.json.
"""
import json, random, subprocess, sys
VCF = "/quobyte/jrigrp/beil/te_evo/singer_analysis/singer-snakemake/for_te_chr10/chr.10.combined.snp.te.sorted.vcf.gz"
SEED, K = 20261111, 6
out = sys.argv[1]
names = subprocess.run(["bcftools", "query", "-l", VCF], check=True, capture_output=True, text=True).stdout.split()
cands = sorted(n for n in names if n != "REF")
held = sorted(random.Random(SEED).sample(cands, K))
panel = [n for n in names if n not in held]          # keeps VCF order, REF first
open(f"{out}/heldout_samples.txt", "w").write("\n".join(held) + "\n")
open(f"{out}/panel_samples.txt", "w").write("\n".join(panel) + "\n")
json.dump({"seed": SEED, "k": K, "n_columns_total": len(names), "n_candidates": len(cands),
           "rng": "random.Random(seed).sample(sorted(candidates), k)", "source_vcf": VCF,
           "heldout": held, "n_panel_columns": len(panel)}, open(f"{out}/selection.json", "w"), indent=1)
print(len(names), len(cands), held, len(panel))
