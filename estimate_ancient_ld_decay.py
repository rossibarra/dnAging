#!/usr/bin/env python3
"""Estimate distance-binned pseudo-haploid LD from ancient VCF columns.

Sites are deterministically subsampled while reading the VCF, then nearby pairs
are sampled uniformly from log-spaced distance bins.  LD is squared Pearson
correlation among individuals called at both sites.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BINS = np.array([1, 10, 25, 50, 100, 250, 500, 1_000, 2_500, 5_000,
                 10_000, 25_000, 50_000, 100_000, 250_000, 500_000, 1_000_000])


def keep_site(chrom: str, pos: int, fraction: float) -> bool:
    value = int.from_bytes(hashlib.blake2b(f"{chrom}:{pos}".encode(), digest_size=8).digest(), "little")
    return value < int(fraction * 2**64)


def parse_gt(field: str) -> int:
    gt = field.split(":", 1)[0]
    if gt == "0":
        return 0
    if gt == "1":
        return 1
    if gt in {"0/0", "0|0"}:
        return 0
    if gt in {"0/1", "1/0", "0|1", "1|0"}:
        return 1
    if gt in {"1/1", "1|1"}:
        return 2
    return -1


def allele_code(ref: str, alt: str) -> int | None:
    """Canonical unordered biallelic SNP code, strand invariant."""
    base = {"A": 0, "C": 1, "G": 2, "T": 3}
    if ref not in base or alt not in base or ref == alt:
        return None
    pair = tuple(sorted((base[ref], base[alt])))
    comp = tuple(sorted((3 - base[ref], 3 - base[alt])))
    return min(pair[0] * 4 + pair[1], comp[0] * 4 + comp[1])


def r2_pair(a: np.ndarray, b: np.ndarray, min_overlap: int) -> tuple[float, int] | None:
    ok = (a >= 0) & (b >= 0)
    n = int(ok.sum())
    if n < min_overlap:
        return None
    x, y = a[ok].astype(float), b[ok].astype(float)
    vx, vy = x.var(), y.var()
    if vx == 0 or vy == 0:
        return None
    r = np.mean((x - x.mean()) * (y - y.mean())) / np.sqrt(vx * vy)
    return float(r * r), n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf", type=Path)
    ap.add_argument("--panel-vcf", type=Path, required=True,
                    help="Modern VCF defining the SNP ascertainment set")
    ap.add_argument("--sample-set", choices=("ancient", "modern"), default="ancient")
    ap.add_argument("--out-prefix", type=Path, required=True)
    ap.add_argument("--site-fraction", type=float, default=0.10)
    ap.add_argument("--pairs-per-bin-chrom", type=int, default=10_000)
    ap.add_argument("--min-overlap", type=int, default=20)
    ap.add_argument("--seed", type=int, default=20260912)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    panel_sites = set()
    panel_records = 0
    with args.panel_vcf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 5)
            panel_records += 1
            code = allele_code(fields[3], fields[4])
            if code is not None:
                panel_sites.add((fields[0], int(fields[1]), code))
    positions: dict[str, list[int]] = defaultdict(list)
    genotypes: dict[str, list[list[int]]] = defaultdict(list)
    sample_names = []
    selected_idx = []
    records = retained = 0
    with args.vcf.open() as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_names = line.rstrip().split("\t")[9:]
                if args.sample_set == "ancient":
                    selected_idx = [i for i, s in enumerate(sample_names) if s.startswith("a")]
                else:
                    selected_idx = [i for i, s in enumerate(sample_names) if not s.startswith("a")]
                continue
            fields = line.rstrip().split("\t")
            records += 1
            chrom, pos = fields[0], int(fields[1])
            code = allele_code(fields[3], fields[4])
            if code is None or (chrom, pos, code) not in panel_sites:
                continue
            if not keep_site(chrom, pos, args.site_fraction):
                continue
            g = [parse_gt(fields[9 + i]) for i in selected_idx]
            called = np.count_nonzero(np.asarray(g) >= 0)
            if called < args.min_overlap or len(set(x for x in g if x >= 0)) < 2:
                continue
            positions[chrom].append(pos)
            genotypes[chrom].append(g)
            retained += 1

    for chrom in positions:
        order = np.argsort(positions[chrom])
        positions[chrom] = np.asarray(positions[chrom], dtype=np.int64)[order]
        genotypes[chrom] = np.asarray(genotypes[chrom], dtype=np.int8)[order]

    rows = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        vals, overlaps = [], []
        for chrom in sorted(positions):
            p, g = positions[chrom], genotypes[chrom]
            if len(p) < 2:
                continue
            attempts = 0
            target = args.pairs_per_bin_chrom
            while len(vals) < target * (sorted(positions).index(chrom) + 1) and attempts < target * 20:
                attempts += 1
                i = int(rng.integers(0, len(p)))
                left = np.searchsorted(p, p[i] + lo, side="left")
                right = np.searchsorted(p, p[i] + hi, side="left")
                if right <= left:
                    continue
                j = int(rng.integers(left, right))
                z = r2_pair(g[i], g[j], args.min_overlap)
                if z is not None:
                    vals.append(z[0]); overlaps.append(z[1])
        a = np.asarray(vals)
        rows.append({"distance_lo": int(lo), "distance_hi": int(hi), "distance_mid": float(np.sqrt(lo * hi)),
                     "pairs": len(a), "mean_r2": float(a.mean()) if len(a) else None,
                     "median_r2": float(np.median(a)) if len(a) else None,
                     "mean_overlap": float(np.mean(overlaps)) if overlaps else None})

    # Unlinked baseline: pairs sampled from different chromosomes.
    unlinked = []
    chroms = [c for c in positions if len(positions[c])]
    for _ in range(100_000):
        c1, c2 = rng.choice(chroms, 2, replace=False)
        g1, g2 = genotypes[c1], genotypes[c2]
        z = r2_pair(g1[rng.integers(len(g1))], g2[rng.integers(len(g2))], args.min_overlap)
        if z is not None:
            unlinked.append(z[0])

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv = args.out_prefix.with_suffix(".tsv")
    with tsv.open("w") as out:
        out.write("distance_lo\tdistance_hi\tdistance_mid\tpairs\tmean_r2\tmedian_r2\tmean_overlap\n")
        for r in rows:
            out.write("\t".join(str(r[k]) for k in r) + "\n")
    summary = {"vcf_records": records, "panel_vcf": str(args.panel_vcf),
               "panel_records": panel_records, "panel_biallelic_snp_keys": len(panel_sites),
               "retained_sites": retained, "sample_set": args.sample_set,
               "samples": len(selected_idx),
               "sample_names": [sample_names[i] for i in selected_idx], "site_fraction": args.site_fraction,
               "min_overlap": args.min_overlap, "unlinked_pairs": len(unlinked),
               "unlinked_mean_r2": float(np.mean(unlinked)), "unlinked_median_r2": float(np.median(unlinked))}
    args.out_prefix.with_name(args.out_prefix.name + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    x = np.array([r["distance_mid"] for r in rows])
    y = np.array([r["mean_r2"] for r in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    label = f"{args.sample_set.capitalize()} samples"
    ax.plot(x, y, marker="o", ms=4, lw=1.5, label=label)
    ax.axhline(np.mean(unlinked), color="0.4", ls="--", lw=1.2, label="Interchromosomal baseline")
    ax.set_xscale("log"); ax.set_xlabel("Physical distance between SNPs (bp)"); ax.set_ylabel(r"Mean $r^2$")
    kind = "Pseudo-haploid" if args.sample_set == "ancient" else "Diploid dosage"
    ax.set_title(f"{kind} LD decay ({len(selected_idx)} {args.sample_set} samples)")
    ax.grid(alpha=.2); ax.legend(frameon=False); fig.tight_layout()
    fig.savefig(args.out_prefix.with_suffix(".png"), dpi=220)


if __name__ == "__main__":
    main()
