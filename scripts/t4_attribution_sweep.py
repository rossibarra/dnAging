"""T4: separate what actually moved the bias.

bias_ideas.md reported -836 / RMSE 890 at Ne=10,000. T1's uniform arm reports
+248 / RMSE 479. Between them landed three independent sets of changes, so the
improvement is real but unattributed. This runs one config per change, on
identical sites, ARGs and table, so the only thing varying is the code path:

  A  original phi_lookup            16-node trapezoid, pre-session bugs present
  B  + my two phi_lookup bug fixes  16-node trapezoid  (isolates bc5268d)
  C  + T1 analytic integration      uniform            (isolates the quadrature)
  D  + T1 den weighting             weighted           (isolates marginalisation)

Site grouping is copied verbatim from scripts/compare_t1_marginalisation.py so
the numbers are directly comparable to T1's. eps = 0 throughout, matching that
harness and the perfect-data scope.
"""
from __future__ import annotations
import argparse, csv, importlib.util, json, sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np, tskit

T4 = Path(__file__).resolve().parent
CONFIGS = {"A_original": ("infer_a_original.py", "uniform"),
           "B_bugfixed": ("infer_b_bugfixed.py", "uniform"),
           "C_analytic": ("infer_c_current.py", "uniform"),
           "D_weighted": ("infer_c_current.py", "weighted")}

def load_mod(fname):
    spec = importlib.util.spec_from_file_location("m_" + fname.replace(".", "_"), T4 / fname)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

def ancient_alt_ids(path):
    alt = set(); rows = 0
    with open(path) as h:
        for line in h:
            if line.startswith("#"): continue
            f = line.rstrip().split("\t")
            gt = f[9].split(":")[f[8].split(":").index("GT")]
            if gt == "1": alt.add(int(f[2]))
            elif gt not in {"0", "."}: raise ValueError(f"bad GT {gt!r}")
            rows += 1
    return alt, rows

def site_groups(ts, alt_ids, cutoff):
    groups = defaultdict(lambda: [0, 0])
    for site in ts.sites():
        if len(site.mutations) != 1: continue
        mut = site.mutations[0]; tree = ts.at(site.position)
        parent = tree.parent(mut.node)
        if parent == tskit.NULL: continue
        d0 = tree.num_samples(mut.node)
        if not 1 <= d0 < ts.num_samples: continue
        lo = float(ts.node(mut.node).time)
        hi = min(float(ts.node(parent).time), cutoff)
        if lo >= cutoff or hi <= lo: continue
        groups[(int(d0), lo, hi)][1 if site.id in alt_ids else 0] += 1
    return groups

def run_one(sim_dir, table_path, ne, truth, config):
    fname, method = CONFIGS[config]
    mod = load_mod(fname)
    stem = Path(sim_dir).name
    ts = tskit.load(Path(sim_dir) / f"{stem}.trees")
    alt_ids, n_vcf = ancient_alt_ids(Path(sim_dir) / f"{stem}_ancient.vcf")
    if n_vcf != ts.num_sites:
        raise ValueError(f"{stem}: {n_vcf} VCF rows vs {ts.num_sites} sites")
    tab = mod.load_table(table_path)
    groups = site_groups(ts, alt_ids, 6.0 * ne)
    ll = np.zeros_like(tab["Tgrid"], dtype=float)
    nan_groups = 0
    for (d0, lo, hi), (n_ref, n_alt) in groups.items():
        kw = {"n_called": ts.num_samples}
        if method != "uniform" or "marginalise" in mod.phi_lookup.__code__.co_varnames:
            kw["marginalise"] = method
        p = mod.phi_lookup(tab, d0, lo, hi, **kw)
        if p is None or not np.all(np.isfinite(p)):
            nan_groups += 1
            continue                      # legacy behaviour: unusable lookup dropped
        p = np.clip(p, 1e-300, 1.0 - 1e-15)
        ll += n_alt * np.log(p) + n_ref * np.log1p(-p)
    summary, _ = mod.summarize(tab["Tgrid"], ll)
    return {"config": config, "ne": ne, "simulation": stem, "true_T": truth,
            "nan_groups": nan_groups, "n_groups": len(groups),
            "sites": sum(sum(v) for v in groups.values()), **summary}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ne", nargs="+", type=int, default=[10000])
    ap.add_argument("--configs", nargs="+", default=list(CONFIGS))
    ap.add_argument("--root", default="/tmp/dnAging-ne{ne}-infinite-10mb")
    ap.add_argument("--table", default="/tmp/dnAging-ne{ne}-t1-format5-table.npz")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    jobs = []
    for ne in a.ne:
        root = Path(a.root.format(ne=ne)); table = a.table.format(ne=ne)
        man = json.loads((root / "run.json").read_text())
        truth = {str(r["simulation"]): float(r["true_age"]) for r in man["replicates"]}
        for sd in sorted(root.glob("simulation_*")):
            k = [sd.name, sd.name.rsplit("_", 1)[-1], str(int(sd.name.rsplit("_", 1)[-1]))]
            tv = next((truth[x] for x in k if x in truth), None)
            if tv is None: raise KeyError(sd.name)
            for cfg in a.configs:
                jobs.append((str(sd), table, ne, tv, cfg))
    print(f"[t4] {len(jobs)} jobs, {a.workers} workers", flush=True)
    out = []
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(run_one, *j): j for j in jobs}
        for fu in as_completed(futs):
            j = futs[fu]
            try:
                r = fu.result(); out.append(r)
                print(f"  {r['config']:<11} ne={r['ne']:<7} {r['simulation']} "
                      f"mean_T={r.get('mean_T',float('nan')):>8.0f} "
                      f"true={r['true_T']:>7.0f} nan_grp={r['nan_groups']}", flush=True)
            except Exception as e:
                print(f"  !! {j[4]} {Path(j[0]).name}: {type(e).__name__}: {e}", flush=True)
    if not out: raise SystemExit("no results")
    out.sort(key=lambda r: (r["config"], r["ne"], r["simulation"]))
    keys = sorted({k for r in out for k in r})
    with open(a.out, "w", newline="") as h:
        w = csv.DictWriter(h, keys, delimiter="\t"); w.writeheader(); w.writerows(out)
    print(f"[t4] wrote {a.out}", flush=True)

if __name__ == "__main__":
    main()
