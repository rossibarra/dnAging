"""VCF reading and ARG-site polarisation for the non-sharing-ratio estimator.

VCFs are read with the stdlib (plain or gzip), biallelic SNPs only.  Genotypes
are stored as allele indices (0 = REF, 1 = ALT, -1 = missing).  `polarise`
then recodes them onto an ARG's sites as derived (1) / ancestral (0) /
missing (-1) by comparing allele *strings* with the ARG site's ancestral state
and derived state, so REF/ALT orientation in the VCF does not matter.

Coordinates: a VCF record at `chrom:pos` maps to ARG site position
`offset(chrom) + pos`, the same convention as `insertion_real_data.run_draw`.
`position_shift` lets a caller correct for 0- vs 1-based ARG input.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class VcfData:
    """Biallelic SNP genotypes for one VCF, grouped by chromosome."""
    columns: list                      # one name per genotype column (sample or haplotype)
    chroms: dict = field(default_factory=dict)   # chrom -> dict(pos, ref, alt, gt)
    skipped: dict = field(default_factory=lambda: {"non_snp_or_multiallelic": 0, "bad_gt": 0})


def _open(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def _split_gt(token, mode):
    """Allele indices (list) for one sample's GT field, or None for a bad field."""
    gt = token.split(":", 1)[0]
    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)
    try:
        alleles = [-1 if a == "." else int(a) for a in parts]
    except ValueError:
        return None
    if any(a > 1 for a in alleles):
        return None
    if len(alleles) == 1:
        return alleles
    if len(alleles) != 2:
        return None
    if mode == "phased":
        if sep != "|":
            return None
        return alleles
    # "sample": one column per sample; hom calls kept, hets and half-missing -> missing
    return [alleles[0] if alleles[0] == alleles[1] else -1]


def read_vcf(path, mode="sample", chroms=None):
    """Read a VCF into `VcfData`.

    mode = "sample": one column per sample.  Haploid calls are used as is; diploid
           homozygous calls give the allele; heterozygous or half-missing calls
           become missing (appropriate for ancient data without phasing).
    mode = "phased": one column per haplotype (`NAME_1`, `NAME_2`); requires phased
           diploid calls (haploid calls give a single column `NAME`).
    """
    if mode not in ("sample", "phased"):
        raise ValueError("mode must be 'sample' or 'phased'")
    wanted = None if chroms is None else {str(c) for c in chroms}
    data = None
    names = None
    buckets = {}
    with _open(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                names = fields[9:]
                data = VcfData(columns=list(names) if mode == "sample" else [])
                continue
            if names is None:
                raise ValueError(f"{path}: no #CHROM header before records")
            chrom = fields[0]
            if wanted is not None and chrom not in wanted:
                continue
            ref, alt = fields[3], fields[4]
            if len(ref) != 1 or len(alt) != 1 or ref == alt or alt in (".", "*"):
                data.skipped["non_snp_or_multiallelic"] += 1
                continue
            cols = []
            ok = True
            for tok in fields[9:]:
                alleles = _split_gt(tok, mode)
                if alleles is None:
                    ok = False
                    break
                cols.append(alleles)
            if not ok:
                data.skipped["bad_gt"] += 1
                continue
            if not data.columns:
                if mode == "sample" or all(len(c) == 1 for c in cols):
                    data.columns = list(names)
                else:
                    data.columns = [f"{n}_{i + 1}" if len(c) == 2 else n
                                    for n, c in zip(names, cols) for i in range(len(c))]
            flat = [a for c in cols for a in c]
            if len(flat) != len(data.columns):
                raise ValueError(f"{path}: ploidy differs between records at {chrom}:{fields[1]}")
            b = buckets.setdefault(chrom, {"pos": [], "ref": [], "alt": [], "gt": []})
            b["pos"].append(int(fields[1]))
            b["ref"].append(ref)
            b["alt"].append(alt)
            b["gt"].append(flat)
    if data is None:
        raise ValueError(f"{path}: no #CHROM header")
    for chrom, b in buckets.items():
        data.chroms[chrom] = {
            "pos": np.asarray(b["pos"], dtype=np.int64),
            "ref": np.asarray(b["ref"]),
            "alt": np.asarray(b["alt"]),
            "gt": np.asarray(b["gt"], dtype=np.int8).T,   # (columns, sites)
        }
    return data


def polarise(vcf, chrom_key_to_layout, positions, ancestral, derived,
             position_shift=0.0):
    """Derived/ancestral/missing matrix of a VCF's columns on ARG sites.

    vcf                 : VcfData
    chrom_key_to_layout : {vcf_chrom_name: (offset, length)} for chromosomes to use
    positions           : (n_sites,) ARG site positions (global, concatenated)
    ancestral, derived  : (n_sites,) allele strings per ARG site ('' if unusable)

    Returns an int8 array (n_columns, n_sites): 1 derived, 0 ancestral, -1 missing
    or the VCF has no usable record at that site, plus a dict of site counts.
    """
    positions = np.asarray(positions, dtype=float)
    n_cols = len(vcf.columns)
    out = np.full((n_cols, positions.size), -1, dtype=np.int8)
    stats = {"arg_sites": int(positions.size), "matched": 0, "allele_mismatch": 0}
    for chrom, (offset, length) in chrom_key_to_layout.items():
        rec = vcf.chroms.get(chrom)
        if rec is None:
            continue
        target = rec["pos"].astype(float) + offset + position_shift
        order = np.searchsorted(positions, target)
        order_clipped = np.minimum(order, positions.size - 1)
        hit = (order < positions.size) & (positions[order_clipped] == target)
        for j in np.nonzero(hit)[0]:
            s = order_clipped[j]
            anc, der = ancestral[s], derived[s]
            if not anc or not der:
                continue
            alleles = {0: rec["ref"][j], 1: rec["alt"][j]}
            if {alleles[0], alleles[1]} != {anc, der}:
                stats["allele_mismatch"] += 1
                continue
            stats["matched"] += 1
            gt = rec["gt"][:, j]
            col = np.full(n_cols, -1, dtype=np.int8)
            for idx in (0, 1):
                if alleles[idx] == der:
                    col[gt == idx] = 1
                elif alleles[idx] == anc:
                    col[gt == idx] = 0
            out[:, s] = col
    return out, stats


def arg_site_alleles(ts):
    """Ancestral and derived allele strings per site ('' if not a single derived state).

    Vectorised over the site/mutation tables (a real genome-wide draw has ~25 M sites).
    A site gets a derived state only if every mutation at it has the same derived state
    and that state differs from the site's ancestral state.
    """
    import tskit

    sites, muts = ts.tables.sites, ts.tables.mutations
    n = sites.num_rows
    anc = np.array(tskit.unpack_strings(sites.ancestral_state, sites.ancestral_state_offset),
                   dtype=object) if n else np.array([], dtype=object)
    der = np.full(n, "", dtype=object)
    if muts.num_rows == 0:
        return anc, der
    mder = np.array(tskit.unpack_strings(muts.derived_state, muts.derived_state_offset), dtype=object)
    uniq, code = np.unique(np.concatenate([mder, anc]), return_inverse=True)
    mcode, acode = code[: mder.size], code[mder.size:]
    msite = np.asarray(muts.site)
    order = np.argsort(msite, kind="stable")
    msite, mcode = msite[order], mcode[order]
    starts = np.concatenate([[0], np.nonzero(np.diff(msite))[0] + 1])
    site_ids = msite[starts]
    lo = np.minimum.reduceat(mcode, starts)
    hi = np.maximum.reduceat(mcode, starts)
    ok = (lo == hi) & (lo != acode[site_ids])
    der[site_ids[ok]] = uniq[lo[ok]]
    return anc, der


# ------------------------------------------------------- panel-fixed sites ---

def scan_panel_fixed(path, chroms=None, min_called=1):
    """Stream a haploid panel VCF and return sites fixed among the called panel haplotypes.

    Only biallelic SNPs where every called haplotype carries the same allele and at
    least `min_called` haplotypes are called are returned.  The VCF is streamed (no
    genotype matrix is kept), so a genome-wide panel VCF is fine.  Diploid-coded
    panel genotypes ('0|0', '1/1') are accepted when homozygous; anything else at a
    site makes that site non-fixed.

    Returns a dict of numpy arrays: chrom (str), pos (int), ref, alt, fixed_idx
    (0 = all REF, 1 = all ALT), n_called, plus 'n_columns' and 'n_records' counters.
    """
    wanted = None if chroms is None else {str(c) for c in chroms}
    out = {k: [] for k in ("chrom", "pos", "ref", "alt", "fixed_idx", "n_called")}
    n_records = 0
    n_columns = None
    with _open(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                n_columns = len(line.rstrip("\n").split("\t")) - 9
                continue
            if n_columns is None:
                raise ValueError(f"{path}: no #CHROM header before records")
            fields = line.rstrip("\n").split("\t")
            if wanted is not None and fields[0] not in wanted:
                continue
            ref, alt = fields[3], fields[4]
            if len(ref) != 1 or len(alt) != 1 or ref == alt or alt in (".", "*"):
                continue
            n_records += 1
            toks = fields[9:]
            n1, n0, nm = toks.count("1"), toks.count("0"), toks.count(".")
            if n0 + n1 + nm != len(toks):          # not plain haploid GT: slow path
                n0 = n1 = 0
                for tok in toks:
                    al = _split_gt(tok, "sample")
                    if al is None or al[0] < 0:
                        continue
                    n0 += al[0] == 0
                    n1 += al[0] == 1
            n_called = n0 + n1
            if n_called < max(1, min_called) or (n0 and n1):
                continue
            out["chrom"].append(fields[0])
            out["pos"].append(int(fields[1]))
            out["ref"].append(ref)
            out["alt"].append(alt)
            out["fixed_idx"].append(1 if n1 else 0)
            out["n_called"].append(n_called)
    res = {k: np.asarray(v) for k, v in out.items()}
    res["pos"] = res["pos"].astype(np.int64)
    res["fixed_idx"] = res["fixed_idx"].astype(np.int8)
    res["n_called"] = res["n_called"].astype(np.int32)
    res["n_columns"] = n_columns
    res["n_records"] = n_records
    return res


def read_polarity_table(path):
    """TSV/space-separated `chrom pos ancestral_allele` -> {(chrom, pos): allele string}."""
    table = {}
    with _open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 3 or f[0].lower() in ("chrom", "chr"):
                continue
            table[(f[0], int(f[1]))] = f[2].upper()
    return table


def polarise_fixed_sites(fixed, chrom_key_to_layout, arg_positions, polarity="ref",
                         polarity_table=None, position_shift=0.0):
    """Turn `scan_panel_fixed` output into extra sites for the estimator.

    Keeps only sites where the panel's fixed allele is *derived* (the only ones that
    contribute (site, reference) pairs), that are not already ARG sites, and whose
    polarity is known:

      polarity = "ref"   : ancestral = VCF REF.  This is the convention of the ARG
                           sites only if the ARG's ancestral state is the reference
                           allele; verify on real data before using it.
      polarity = "table" : ancestral allele from `polarity_table` ((chrom, pos) ->
                           allele; sites missing from the table are dropped).

    Returns dict(positions, ancestral, derived, n_called, chrom) sorted by position,
    plus stats.
    """
    arg_positions = np.asarray(arg_positions, dtype=float)
    stats = {"fixed_in_panel": int(fixed["pos"].size), "fixed_ancestral_dropped": 0,
             "already_in_arg": 0, "no_polarity": 0, "kept": 0, "polarity": polarity,
             "allele_mismatch_with_table": 0}
    keep = {"pos": [], "anc": [], "der": [], "n": [], "chrom": []}
    for i in range(fixed["pos"].size):
        chrom = str(fixed["chrom"][i])
        if chrom not in chrom_key_to_layout:
            continue
        offset = chrom_key_to_layout[chrom][0]
        gpos = float(fixed["pos"][i]) + offset + position_shift
        ref, alt, k = str(fixed["ref"][i]), str(fixed["alt"][i]), int(fixed["fixed_idx"][i])
        fixed_allele, other = (alt, ref) if k == 1 else (ref, alt)
        if polarity == "ref":
            anc, der = ref, alt
        elif polarity == "table":
            a = (polarity_table or {}).get((chrom, int(fixed["pos"][i])))
            if a is None:
                stats["no_polarity"] += 1
                continue
            if a not in (ref, alt):
                stats["allele_mismatch_with_table"] += 1
                continue
            anc, der = a, (alt if a == ref else ref)
        else:
            raise ValueError("polarity must be 'ref' or 'table'")
        if fixed_allele != der:
            stats["fixed_ancestral_dropped"] += 1
            continue
        keep["pos"].append(gpos); keep["anc"].append(anc); keep["der"].append(der)
        keep["n"].append(int(fixed["n_called"][i])); keep["chrom"].append(chrom)
    pos = np.asarray(keep["pos"], dtype=float)
    # drop positions already present among the ARG sites
    if pos.size and arg_positions.size:
        idx = np.minimum(np.searchsorted(arg_positions, pos), arg_positions.size - 1)
        dup = arg_positions[idx] == pos
    else:
        dup = np.zeros(pos.size, dtype=bool)
    stats["already_in_arg"] = int(dup.sum())
    order = np.argsort(pos[~dup], kind="stable")
    sel = lambda key: np.asarray(keep[key], dtype=object)[~dup][order]
    stats["kept"] = int((~dup).sum())
    return {"positions": pos[~dup][order], "ancestral": sel("anc"), "derived": sel("der"),
            "n_called": np.asarray(keep["n"], dtype=np.int32)[~dup][order],
            "chrom": sel("chrom"), "stats": stats}
