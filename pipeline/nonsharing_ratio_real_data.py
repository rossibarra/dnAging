#!/usr/bin/env python3
"""Non-sharing-ratio dating of ancient samples from real data.

Inputs
  --tree          one SINGER posterior draw (tszip or .trees), built on the panel only
  --controls-vcf  held-out modern genomes (NOT in the ARG panel)
  --ancient-vcf   ancient samples to date
  c(t)            by default from the ARG's own lineage counts; --demography FILE.yaml
                  (demes format) or --ne FILE.tsv (the repo's Ne trajectory TSV) use
                  1 / (ploidy * N(t)) instead

Subcommands, mirroring insertion_real_data.py:
  prepare-fixed  once: stream the panel VCF and list sites fixed among panel haplotypes
  run-draw       one ARG draw -> NPZ of per-block (site, reference) pair counts per target
  merge          NPZs from all draws -> per-sample T_hat, spread over draws, block-bootstrap CI

Method and design rules: see pipeline/nonsharing_ratio.py.  Panel-fixed sites that are
absent from the ARG's site list are recovered with `prepare-fixed` + `run-draw --panel-fixed`
(see docs/nonsharing_ratio/RUNNING.md for the polarity limitation).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import nonsharing_ratio as nr
import nonsharing_ratio_demography as nd
import nonsharing_ratio_io as nio


def load_tree(path):
    path = Path(path)
    if path.suffix == ".trees":
        import tskit
        return tskit.load(str(path))
    import tszip
    return tszip.decompress(str(path))


def layout_for(ts, chrom_fai):
    from insertion_real_data import chromosome_layout
    return chromosome_layout(ts, chrom_fai)


def match_vcf_chromosomes(layout, vcf):
    """Map each ARG chromosome key to the VCF's name for it ('1' vs 'chr1')."""
    out = {}
    for key in layout:
        for cand in (key, key[3:] if key.startswith("chr") else "chr" + key):
            if cand in vcf.chroms:
                out[cand] = layout[key]
                break
    return out


def block_layout(layout, block_size):
    """Ordered blocks: (chrom, index) labels, and per-chromosome first block id."""
    labels, first = [], {}
    for chrom in sorted(layout, key=lambda c: layout[c][0]):
        first[chrom] = len(labels)
        for b in range(int(np.ceil(layout[chrom][1] / block_size))):
            labels.append(f"{chrom}:{b}")
    return labels, first


def site_blocks(positions, layout, block_size, first):
    chroms = sorted(layout, key=lambda c: layout[c][0])
    offsets = np.array([layout[c][0] for c in chroms])
    idx = np.searchsorted(offsets, positions, side="right") - 1
    local = positions - offsets[idx]
    base = np.array([first[chroms[i]] for i in idx])
    return base + (local // block_size).astype(int)


def vcf_names(layout):
    """All spellings of the ARG chromosome keys a VCF might use ('1' / 'chr1')."""
    names = set()
    for key in layout:
        names.update({key, key[3:] if key.startswith('chr') else 'chr' + key})
    return names


def check_disjoint(names_a, names_b, label):
    clash = sorted(set(names_a) & set(names_b))
    if clash:
        raise ValueError(f"{label} share sample names: {clash[:5]}")


def restrict_to_chromosomes(ts, layout, chroms):
    """Keep only the listed chromosomes (ARG keys); returns (ts, layout subset)."""
    unknown = [c for c in chroms if c not in layout]
    if unknown:
        raise ValueError(f"chromosomes {unknown} not in the ARG layout {sorted(layout)}")
    intervals = [(layout[c][0], layout[c][0] + layout[c][1]) for c in chroms]
    sub = ts.keep_intervals(intervals, simplify=False)
    return sub, {c: layout[c] for c in chroms}


def add_fixed_sites(positions, anc_state, der_state, panel, ages, fixed_sites):
    """Append panel-fixed sites (derived in every called panel haplotype, age +inf).

    Extra sites are given a panel column with n_called ones (the called haplotypes carry
    the derived allele; uncalled ones count as non-references), so the estimator's
    'segregating and older than A0' rule keeps them through age = +inf.  Returns the
    merged arrays sorted by position and a boolean `is_extra`.
    """
    n_panel = panel.shape[0]
    k = fixed_sites["positions"].size
    if k == 0:
        return positions, anc_state, der_state, panel, ages, np.zeros(positions.size, bool)
    extra_panel = (np.arange(n_panel)[:, None] < np.minimum(fixed_sites["n_called"], n_panel)[None, :]
                   ).astype(np.int8)
    pos = np.concatenate([positions, fixed_sites["positions"]])
    order = np.argsort(pos, kind="stable")
    anc = np.concatenate([anc_state, fixed_sites["ancestral"]])[order]
    der = np.concatenate([der_state, fixed_sites["derived"]])[order]
    pan = np.concatenate([panel, extra_panel], axis=1)[:, order]
    ag = np.concatenate([ages, np.full(k, np.inf)])[order]
    is_extra = np.concatenate([np.zeros(positions.size, bool), np.ones(k, bool)])[order]
    return pos[order], anc, der, pan, ag, is_extra


def load_rates(args, ts, edges):
    if args.demography and args.ne:
        raise ValueError("give at most one of --demography and --ne")
    if args.demography:
        return nd.rates_from_demes(args.demography, edges, deme=args.deme,
                                   ploidy=args.ploidy, generation_time=args.generation_time), str(args.demography)
    if args.ne:
        return nd.rates_from_ne_tsv(args.ne, edges, series=args.ne_series, ploidy=args.ploidy), str(args.ne)
    return nd.rates_from_arg(ts, edges), "arg"


def prepare_fixed(args):
    import time
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    t = time.time()
    fixed = nio.scan_panel_fixed(args.panel_vcf, chroms=args.chroms,
                                 min_called=args.min_called)
    np.savez(args.output, **{k: v for k, v in fixed.items() if k not in ("n_columns", "n_records")},
             meta=json.dumps({"panel_vcf": str(args.panel_vcf), "n_columns": fixed["n_columns"],
                              "n_records": fixed["n_records"], "min_called": args.min_called,
                              "seconds": round(time.time() - t, 1)}))
    print(f"scanned {fixed['n_records']} biallelic SNP records, {fixed['pos'].size} fixed among "
          f">= {args.min_called} called of {fixed['n_columns']} panel columns "
          f"({time.time() - t:.0f}s)", flush=True)


def run_draw(args):
    import time
    t0 = time.time()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    ts = load_tree(args.tree)
    layout = layout_for(ts, args.chrom_fai)
    if args.chroms:
        ts, layout = restrict_to_chromosomes(ts, layout, args.chroms)
    controls = nio.read_vcf(args.controls_vcf, mode=args.controls_mode, chroms=vcf_names(layout))
    ancient = nio.read_vcf(args.ancient_vcf, mode=args.ancient_mode, chroms=vcf_names(layout))
    check_disjoint(controls.columns, ancient.columns, "controls and ancient VCFs")
    if args.panel_names:
        panel_names = [x.strip() for x in Path(args.panel_names).read_text().split() if x.strip()]
        check_disjoint(panel_names, controls.columns, "ARG panel and controls")
        check_disjoint(panel_names, ancient.columns, "ARG panel and ancient samples")
    print(f"read VCFs ({time.time() - t0:.0f}s)", flush=True)

    positions = np.asarray(ts.tables.sites.position)
    anc_state, der_state = nio.arg_site_alleles(ts)
    panel = nr.panel_derived_matrix(ts)
    ages = nr.edge_child_ages(ts) if args.age_mode == "child" else nr.edge_midpoint_ages(ts)
    print(f"ARG tables ({time.time() - t0:.0f}s): {positions.size} sites", flush=True)

    fixed_stats = None
    is_extra = np.zeros(positions.size, dtype=bool)
    if args.panel_fixed:
        raw = np.load(args.panel_fixed, allow_pickle=False)
        fixed = {k: raw[k] for k in ("chrom", "pos", "ref", "alt", "fixed_idx", "n_called")}
        table = nio.read_polarity_table(args.polarity_table) if args.polarity_table else None
        if args.fixed_polarity == "table" and table is None:
            raise ValueError("--fixed-polarity table needs --polarity-table")
        # VCF chromosome names of the panel VCF -> ARG keys
        fixed_layout = {}
        for key in layout:
            for cand in (key, key[3:] if key.startswith("chr") else "chr" + key):
                fixed_layout[cand] = layout[key]
        extra = nio.polarise_fixed_sites(fixed, fixed_layout, positions, args.fixed_polarity,
                                         table, args.position_shift)
        fixed_stats = extra["stats"]
        positions, anc_state, der_state, panel, ages, is_extra = add_fixed_sites(
            positions, anc_state, der_state, panel, ages, extra)
        print(f"panel-fixed sites added: {fixed_stats}", flush=True)

    ctrl_d, ctrl_stats = nio.polarise(controls, match_vcf_chromosomes(layout, controls),
                                      positions, anc_state, der_state, args.position_shift)
    anc_d, anc_stats = nio.polarise(ancient, match_vcf_chromosomes(layout, ancient),
                                    positions, anc_state, der_state, args.position_shift)
    print(f"polarised ({time.time() - t0:.0f}s)", flush=True)

    usable = (panel >= 0).all(axis=0) & (der_state != "")
    targets = np.vstack([anc_d, ctrl_d])[:, usable]
    called = targets >= 0
    targets = np.where(called, targets, 0)
    panel, ages, positions, is_extra = panel[:, usable], ages[usable], positions[usable], is_extra[usable]

    edges = nd.make_edges(args.c_grid_max, args.c_grid_n)
    (edges, rates), c_source = load_rates(args, ts, edges)

    mask = nr.select_sites(panel, targets, ages, args.a0,
                           include_panel_fixed=not args.no_panel_fixed, called=called)
    labels, first = block_layout(layout, args.block_size)
    blocks = site_blocks(positions, layout, args.block_size, first)
    num, den = nr.nonsharing_counts(panel, targets, mask, called, blocks, len(labels))

    n_anc = anc_d.shape[0]
    names = list(ancient.columns) + list(controls.columns)
    anc_names = list(ancient.columns)
    if args.deidentify:   # S01.. / C01.. in VCF column order; nothing else leaves the VCF
        anc_names = [f"S{i + 1:02d}" for i in range(n_anc)]
        names = anc_names + [f"C{i + 1:02d}" for i in range(len(controls.columns))]
    is_control = np.arange(len(names)) >= n_anc
    T, H, f = nr.ages_from_counts(num, den, is_control, edges, rates)
    flags = nr.target_flags(f, H, den.sum(axis=0)[~is_control])
    meta = {"tree": str(args.tree), "draw_id": args.draw_id, "a0": args.a0, "age_mode": args.age_mode,
            "c_source": c_source, "chroms": list(layout),
            "position_shift": args.position_shift,
            "n_arg_sites": int(usable.size), "n_usable_sites": int(usable.sum()),
            "n_selected_sites": int(mask.sum()),
            "n_selected_extra_fixed_sites": int((mask & is_extra).sum()),
            "n_panel": int(panel.shape[0]),
            "controls_polarise": ctrl_stats, "ancient_polarise": anc_stats,
            "controls_skipped": controls.skipped, "ancient_skipped": ancient.skipped,
            "include_panel_fixed": not args.no_panel_fixed,
            "panel_fixed_input": fixed_stats, "seconds": round(time.time() - t0, 1)}
    np.savez(args.output, targets=np.array(names), is_control=is_control, num=num, den=den,
             edges=edges, rates=rates, block_labels=np.array(labels), T_hat=T,
             flags=np.array(flags, dtype=str), meta=json.dumps(meta))
    for n, t, h, fl in zip(anc_names, T, H, flags):
        print(f"{n}\tT_hat={t:.0f}\tH={h:.5f}" + (f"\tFLAG={fl}" if fl else ""))
    if any(flags):
        print(f"WARNING: {int(np.count_nonzero(flags))} ancient target(s) flagged "
              f"({sorted(set(f for f in flags if f))}); T_hat is not a valid estimate there",
              flush=True)
    print(f"selected {int(mask.sum())} of {int(usable.sum())} sites "
          f"({int((mask & is_extra).sum())} added panel-fixed); {time.time() - t0:.0f}s", flush=True)


def merge(args):
    draws = [np.load(p, allow_pickle=False) for p in args.draw_results]
    if not draws:
        raise ValueError("no draw results given")
    ref = draws[0]
    for p, d in zip(args.draw_results, draws):
        for key in ("targets", "is_control", "block_labels"):
            if not np.array_equal(ref[key], d[key]):
                raise ValueError(f"{p}: {key} differs from the first draw")
    is_control = ref["is_control"]
    names = [str(x) for x in ref["targets"][~is_control]]
    n_blocks = ref["num"].shape[0]

    def combined(idx):
        per_draw = np.array([nr.ages_from_counts(d["num"][idx], d["den"][idx], is_control,
                                                 d["edges"], d["rates"])[0] for d in draws])
        return per_draw

    base = combined(slice(None))
    rng = np.random.default_rng(args.bootstrap_seed)
    boot = np.empty((args.bootstrap_replicates, base.shape[1]))
    for b in range(args.bootstrap_replicates):
        boot[b] = combined(rng.integers(0, n_blocks, n_blocks)).mean(axis=0)
    lo, hi = np.nanpercentile(boot, [2.5, 97.5], axis=0)
    # degenerate-case flags: per draw on all blocks, and over the pooled counts
    per_draw_flags = []
    for d in draws:
        T_, H_, f_ = nr.ages_from_counts(d["num"], d["den"], is_control, d["edges"], d["rates"])
        per_draw_flags.append(nr.target_flags(f_, H_, np.asarray(d["den"]).sum(axis=0)[~is_control]))
    per_draw_flags = np.array(per_draw_flags)
    with args.output.open("w") as fh:
        fh.write("sample\tT_hat_mean\tT_hat_median\tdraw_sd\tci_lo\tci_hi\tn_draws\tflag\n")
        for j, n in enumerate(names):
            kinds = sorted({x for x in per_draw_flags[:, j] if x})
            flag = ",".join(f"{k}:{int((per_draw_flags[:, j] == k).sum())}/{len(draws)}" for k in kinds)
            fh.write(f"{n}\t{base[:, j].mean():.1f}\t{np.median(base[:, j]):.1f}\t"
                     f"{base[:, j].std(ddof=1) if len(draws) > 1 else float('nan'):.1f}\t"
                     f"{lo[j]:.1f}\t{hi[j]:.1f}\t{len(draws)}\t{flag}\n")
    n_flagged = int((per_draw_flags != "").any(axis=0).sum())
    if n_flagged:
        print(f"WARNING: {n_flagged} of {len(names)} ancient targets are flagged in at least one "
              "draw (f0 = f is zero, H<=0 = no excess non-sharing, no_data); see the flag column")
    print(f"wrote {args.output}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run-draw", help="process one ARG draw")
    r.add_argument("--tree", type=Path, required=True)
    r.add_argument("--draw-id", type=int, required=True)
    r.add_argument("--controls-vcf", type=Path, required=True)
    r.add_argument("--ancient-vcf", type=Path, required=True)
    r.add_argument("--controls-mode", choices=["sample", "phased"], default="phased",
                   help="phased: one column per haplotype of phased diploid controls [phased]")
    r.add_argument("--ancient-mode", choices=["sample", "phased"], default="sample",
                   help="sample: hets/half-missing become missing [sample]")
    r.add_argument("--panel-names", type=Path, help="text file of the ARG panel's sample names; "
                   "used only to verify no control/ancient sample is in the panel")
    r.add_argument("--chrom-fai", type=Path,
                   help="FAI to reconstruct chromosome offsets if the ARG lacks chrom_offsets metadata")
    r.add_argument("--position-shift", type=float, default=0.0,
                   help="added to VCF position + chromosome offset to match ARG site positions. 0 is right "
                        "for the real maize ARG draws and VCFs (checked on the real data) [0]")
    r.add_argument("--a0", type=float, required=True,
                   help="keep segregating sites dated above this age (generations)")
    r.add_argument("--age-mode", choices=["midpoint", "child"], default="midpoint",
                   help="site age compared with --a0: midpoint of the mutation edge (default) or the "
                        "time of its child node (lower bound of the edge)")
    r.add_argument("--no-panel-fixed", action="store_true",
                   help="do not include panel-fixed derived sites (design rule 2 says to include them)")
    r.add_argument("--chroms", nargs="+", help="restrict to these ARG chromosome keys (smoke tests); "
                   "the ARG-derived c(t) is then computed on that subset too")
    r.add_argument("--panel-fixed", type=Path,
                   help="NPZ from `prepare-fixed`: panel-fixed sites absent from the ARG (design rule 2)")
    r.add_argument("--fixed-polarity", choices=["ref", "table"], default="table",
                   help="polarity of --panel-fixed sites: table = ancestral allele from --polarity-table "
                        "(e.g. an outgroup); ref = ancestral is the VCF REF allele, which is NOT valid for the "
                        "maize ARG (its ancestral state equals the panel VCF REF at only ~69%% of sites) "
                        "[table]")
    r.add_argument("--polarity-table", type=Path,
                   help="text file 'chrom pos ancestral_allele' (e.g. from an outgroup) for --fixed-polarity table")
    r.add_argument("--demography", type=Path, help="demes YAML; default c(t) from the ARG lineage counts")
    r.add_argument("--ne", type=Path, help="Ne trajectory TSV (time_left/time_right/effective_population_size); "
                   "c = 1/(ploidy*N)")
    r.add_argument("--ne-series", help="series name in the Ne TSV (e.g. posterior_mean)")
    r.add_argument("--deme", help="deme name when the demes file has several")
    r.add_argument("--ploidy", type=int, default=2, help="ploidy for c = 1/(ploidy*N) from demes [2]")
    r.add_argument("--generation-time", type=float, help="override the demes generation time (years)")
    r.add_argument("--c-grid-max", type=float, default=20000.0)
    r.add_argument("--c-grid-n", type=int, default=100)
    r.add_argument("--block-size", type=float, default=5_000_000)
    r.add_argument("--deidentify", action="store_true",
                   help="write ancient samples as S01.. and controls as C01.. (VCF column order) "
                        "in the NPZ and stdout")
    r.add_argument("--output", type=Path, required=True)
    r.set_defaults(func=run_draw)
    pf = sub.add_parser("prepare-fixed", help="list panel-fixed sites from the panel VCF (once)")
    pf.add_argument("--panel-vcf", type=Path, required=True)
    pf.add_argument("--chroms", nargs="+", help="VCF chromosome names to scan [all]")
    pf.add_argument("--min-called", type=int, default=21,
                   help="minimum panel haplotypes called at a site for it to count as fixed [21 = 80%% of 26]")
    pf.add_argument("--output", type=Path, required=True)
    pf.set_defaults(func=prepare_fixed)
    m = sub.add_parser("merge", help="combine draws and bootstrap")
    m.add_argument("--draw-results", type=Path, nargs="+", required=True)
    m.add_argument("--bootstrap-replicates", type=int, default=1000)
    m.add_argument("--bootstrap-seed", type=int, default=20261111)
    m.add_argument("--output", type=Path, required=True)
    m.set_defaults(func=merge)
    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
