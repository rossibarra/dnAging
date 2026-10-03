#!/usr/bin/env python3
"""Prepare the already-filtered real-data sites and 20 ancient calls.

The completed diffusion run saved precisely the 10-draw-complete biallelic sites
and effective pseudo-haploid calls used in inference.  This adapter retains that
site ascertainment, subsets the requested cohort, and adds ancient-VCF REF/ALT so
each beta-binomial ARG draw can orient the observation to its mutation state.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epsilon-data-root", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--ancient-vcf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    wanted_samples = args.samples.read_text().split()
    targets = {}
    staged = {}
    for chrom in map(str, range(1, 11)):
        src = np.load(args.epsilon_data_root / chrom / "epsilon_calibration_data.npz",
                      allow_pickle=True)
        pos = np.asarray(src["position"], dtype=np.int64)
        targets[chrom] = set(pos.tolist())
        staged[chrom] = {
            "position": pos,
            "records": {},
        }

    sample_columns = None
    with args.ancient_vcf.open() as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                names = line.rstrip("\n").split("\t")[9:]
                missing = [sample for sample in wanted_samples if sample not in names]
                if missing:
                    raise SystemExit(f"samples absent from ancient VCF: {missing}")
                sample_columns = [9 + names.index(sample) for sample in wanted_samples]
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            if chrom not in targets:
                continue
            pos = int(fields[1])
            if pos in targets[chrom]:
                ref, alt = fields[3].upper(), fields[4].upper()
                if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
                    if sample_columns is None:
                        raise SystemExit("VCF data encountered before #CHROM header")
                    observed_alt = []
                    called = []
                    for column in sample_columns:
                        gt = fields[column].split(":", 1)[0]
                        alleles = [int(value) for value in re.split(r"[/|]", gt)
                                   if value in {"0", "1"}]
                        called.append(int(bool(alleles)))
                        observed_alt.append(int(any(value == 1 for value in alleles)))
                    staged[chrom]["records"][pos] = (
                        ref, alt, observed_alt, called)

    args.output.mkdir(parents=True, exist_ok=False)
    counts = {}
    for chrom, data in staged.items():
        pos = data["position"]
        absent = [int(p) for p in pos if int(p) not in data["records"]]
        if absent:
            raise SystemExit(f"chr{chrom}: {len(absent)} target alleles absent from VCF")
        ref = np.array([data["records"][int(p)][0] for p in pos], dtype="U1")
        alt = np.array([data["records"][int(p)][1] for p in pos], dtype="U1")
        observed_alt = np.array(
            [data["records"][int(p)][2] for p in pos], dtype=np.uint8)
        called = np.array(
            [data["records"][int(p)][3] for p in pos], dtype=np.uint8)
        np.savez(args.output / f"chr{chrom}.npz", position=pos, ref=ref, alt=alt,
                 observed_alt=observed_alt, called=called,
                 samples=np.asarray(wanted_samples))
        counts[chrom] = len(pos)
    (args.output / "metadata.json").write_text(json.dumps({
        "ancient_vcf": str(args.ancient_vcf.resolve()),
        "epsilon_data_root": str(args.epsilon_data_root.resolve()),
        "samples_file": str(args.samples.resolve()),
        "n_samples": len(wanted_samples),
        "sites_by_chromosome": counts,
        "site_policy": "10-draw-complete biallelic positions from prior diffusion run",
        "call_policy": "direct VCF GT; homozygous representation collapsed to one pseudo-haploid call",
    }, indent=2) + "\n")
    print(f"wrote {args.output}; {sum(counts.values())} sites", flush=True)


if __name__ == "__main__":
    main()
