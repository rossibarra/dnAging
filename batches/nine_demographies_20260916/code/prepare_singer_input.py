#!/usr/bin/env python3
"""Build SINGER input from a variable-Ne replicate's modern samples.

SINGER must see only the modern panel: it is the ARG the inference treats as
known, and the ancient haplotypes are the thing being dated.  This drops the
ancient columns, leaving 26 haploid samples, which satisfies the workflow's
requirement of an even number of haplotypes.

Two details that would otherwise corrupt the run:

* **Duplicate positions.**  The simulations use a continuous genome, and the VCF
  writer truncates positions to integers, so distinct sites collide at a rate of
  roughly n^2 / 2L -- measured at 0.35% across the first three replicates.  Each
  duplicate is nudged forward to the next free integer rather than dropped: no
  site is lost, ordering is preserved, and a shift of a few base pairs at
  r = 1e-8 is far below the resolution of any tree the ARG can distinguish.
  The original position and the site id are both recorded in `position_map.tsv`
  so the inferred ARG can be joined back to simulation truth either way.

* **Polarisation.**  These VCFs are written with the ancestral state in REF, so
  the config sets `polarised: true`.  That is the same convention the earlier
  `simarg` runs used.

The ancient calls carry 1% genotyping error, but that never reaches SINGER --
error enters only through the ancient columns, which are removed here.  The ARG
is therefore estimated from error-free modern data, matching the real analysis
where the panel is modern sequence and the error is a property of aDNA.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


def prepare(args):
    rep_dir = args.simulations / f"replicate_{args.replicate:03d}"
    metadata = json.loads((rep_dir / "metadata.json").read_text())
    n_modern = int(metadata["n_modern_haploid"])

    # --input-dir puts several replicates side by side as chromosomes of one
    # dataset; without it each replicate gets its own per-replicate input tree.
    out_dir = (args.input_dir if args.input_dir is not None
               else args.output_root / f"replicate_{args.replicate:03d}" / "input")
    out_dir.mkdir(parents=True, exist_ok=True)
    vcf_out = out_dir / f"{args.chrom}.vcf.gz"

    used = set()
    shifted = 0
    rows = []
    names = None
    with gzip.open(rep_dir / "all_samples.vcf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                names = fields[9:]
                modern = [i for i, s in enumerate(names) if s.startswith("modern_")]
                if len(modern) != n_modern:
                    raise SystemExit(
                        f"expected {n_modern} modern columns, found {len(modern)}")
                continue
            pos = int(fields[1])
            while pos in used:               # nudge, do not drop
                pos += 1
                shifted += 1
            used.add(pos)
            rows.append((pos, int(fields[1]), fields[2], fields[3], fields[4],
                         [fields[9 + i] for i in modern]))

    if names is None:
        raise SystemExit("no #CHROM header found")
    rows.sort(key=lambda r: r[0])
    sample_names = [names[i] for i in modern]

    with gzip.open(vcf_out, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write(f"##contig=<ID={args.chrom},length={int(metadata['sequence_length'])}>\n")
        handle.write("#" + "\t".join(
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
             "FORMAT", *sample_names]) + "\n")
        for pos, _orig, site_id, ref, alt, gts in rows:
            handle.write("\t".join(
                [args.chrom, str(pos), site_id, ref, alt, ".", "PASS", ".",
                 "GT", *gts]) + "\n")

    with (out_dir / f"{args.chrom}.meta.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "population"])
        for name in sample_names:
            writer.writerow([name, "modern"])

    with (out_dir / f"{args.chrom}.position_map.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["site_id", "original_pos", "singer_pos"])
        for pos, orig, site_id, _ref, _alt, _gts in rows:
            writer.writerow([site_id, orig, pos])

    summary = {
        "replicate": args.replicate,
        "chrom": args.chrom,
        "modern_samples": len(sample_names),
        "sites": len(rows),
        "positions_shifted": shifted,
        "sequence_length": metadata["sequence_length"],
        "true_ancient_ages": metadata["true_ancient_ages"],
        "epochs": metadata["epochs"],
        "source_vcf": str(rep_dir / "all_samples.vcf.gz"),
        "note": "ancient columns removed; SINGER sees error-free modern panel only",
    }
    (out_dir / f"{args.chrom}.prepare_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--simulations", type=Path,
                   default=Path("msprime_variable_ne_error/simulations"))
    p.add_argument("--output-root", type=Path,
                   default=Path("msprime_variable_ne_error/singer"))
    p.add_argument("--replicate", type=int, required=True)
    p.add_argument("--chrom", default="chr1")
    p.add_argument("--input-dir", type=Path, default=None,
                   help="write directly into this SINGER input directory, so "
                        "multiple replicates can act as chromosomes of one run")
    prepare(p.parse_args())
