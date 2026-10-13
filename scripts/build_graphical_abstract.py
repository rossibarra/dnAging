#!/usr/bin/env python3
"""Build the graphical abstract for the ancient-lineage insertion method.

Panels (left to right), each with its equation from MATH.md section 8:
  1  a real msprime modern tree with three mutation-bearing edges and an
     ancient lineage inserted at age T (eq. 13);
  2  the per-site derived probability q_i(T) for those three edges, computed
     with pipeline/insertion_likelihood.py (eq. 14-15);
  3  the posterior for one ancient sample from the 300-replicate constant-Ne
     benchmark, inferred with the true modern ARG (eq. 3a);
  4  MAP versus true age for the constant-Ne and variable-Ne benchmarks.

--layout slide (default): 16:9, one row; writes docs/graphical_abstract_slide.{png,pdf,svg}.
--layout readme: 2 x 2 grid; writes docs/graphical_abstract_readme.png.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import msprime
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
from insertion_likelihood import derived_probability_uniform_edge_grid  # noqa: E402

INK = "#17202B"
INK2 = "#3B4552"
MUTED = "#5A6472"
GRID = "#E1E5EA"
PANEL = "#F4F6F8"
TREE = "#8C96A3"
BLUE = "#2563A8"
ORANGE = "#D9731F"
# sequential orange ramp for the three example edges (young -> old)
EDGE_COLORS = ["#F0A766", "#D9731F", "#8F4710"]

BATCH300 = ROOT / "msprime_batch300_insertion"
VARNE = ROOT / "msprime_variable_ne_error/insertion_results/age_demography_summary.tsv"

TREE_NE = 10_000
TREE_SEED = 10
# Hand-picked illustration of the known-ARG case: true age near 5,000 gen and
# a 95% interval that contains it (46% of the 300 replicates do).
EXAMPLE_REPLICATE = 27

AGE_TICKS = [0, 2500, 5000, 7500, 10000]
AGE_LABELS = ["0", "2.5k", "5k", "7.5k", "10k"]


def read_tsv(path):
    return np.genfromtxt(path, delimiter="\t", names=True, dtype=None, encoding=None)


def simulate_tree():
    ts = msprime.sim_ancestry(
        samples=4, ploidy=2, population_size=TREE_NE, random_seed=TREE_SEED,
        sequence_length=1,
    )
    return ts.first()


def node_x(tree):
    x = {leaf: float(i) for i, leaf in enumerate(tree.leaves())}
    for node in tree.nodes(order="postorder"):
        if node not in x:
            x[node] = float(np.mean([x[c] for c in tree.children(node)]))
    return x


def pick_edges(tree):
    """Youngest internal edge for each of three distinct clade sizes."""
    t = tree.time
    internal = [(c, tree.parent(c)) for c in tree.nodes()
                if tree.parent(c) != -1 and not tree.is_leaf(c)]
    by_size = {}
    for edge in sorted(internal, key=lambda e: t(e[0])):
        by_size.setdefault(tree.num_samples(edge[0]), edge)
    sizes = sorted(by_size)
    picks = [sizes[0], sizes[len(sizes) // 2], sizes[-1]]
    return [by_size[k] for k in picks]


def kilo_ticks(ax, axis, top):
    step = 10_000 if top > 30_000 else 5_000
    ticks = np.arange(0, top + 1, step)
    labels = [f"{v / 1000:g}k" if v else "0" for v in ticks]
    if axis == "x":
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
    else:
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)


def style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=11)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.xaxis.label.set_size(12)
    ax.yaxis.label.set_size(12)
    ax.set_facecolor("white")


def draw_tree(ax, tree, edges, sample_T, insert):
    t = tree.time
    x = node_x(tree)
    root = tree.root
    marked = {child: i for i, (child, _) in enumerate(edges)}
    for child in tree.nodes():
        parent = tree.parent(child)
        if parent == -1:
            continue
        color, lw, z = TREE, 2.4, 2
        if child in marked:
            color, lw, z = EDGE_COLORS[marked[child]], 4.5, 3
        ax.plot([x[child], x[child], x[parent]], [t(child), t(parent), t(parent)],
                color=color, lw=lw, solid_capstyle="round", zorder=z)
    ax.plot([x[root], x[root]], [t(root), t(root) * 1.06], color=TREE, lw=2.4)
    for i, (child, parent) in enumerate(edges):
        ym = 0.5 * (t(child) + t(parent))
        ax.plot(x[child], ym, marker="*", ms=17, color=EDGE_COLORS[i],
                mec="white", mew=1.2, zorder=5)
        ax.text(x[child] + 0.2, ym + 0.03 * t(root), "abc"[i], fontsize=15,
                style="italic", color=INK, family="serif", fontweight="bold",
                va="center", zorder=6)
    leaves = list(tree.leaves())
    ax.scatter([x[l] for l in leaves], [0] * len(leaves), s=34, color=TREE,
               zorder=4)

    # ancient lineage: sampled at T, first joins the tree at time u
    child, parent = insert
    u = 0.5 * (max(sample_T, t(child)) + t(parent))
    xa = max(x.values()) + 0.9
    ax.plot([xa, xa, x[child]], [sample_T, u, u], color=BLUE, lw=2.8,
            ls=(0, (4, 3)), zorder=6)
    ax.plot(xa, sample_T, "o", ms=12, color=BLUE, mec="white", mew=1.5, zorder=7)
    ax.plot(x[child], u, "o", ms=7, color=BLUE, mec="white", mew=1.2, zorder=7)
    ax.text(xa + 0.25, sample_T, "ancient\nsample", fontsize=12, color=BLUE,
            va="center")
    ax.text(xa + 0.25, u, r"$u$", fontsize=16, color=BLUE, va="center")
    ax.axhline(sample_T, color=BLUE, lw=0.9, ls=":", alpha=0.6, zorder=1)
    ax.text(-0.5, sample_T, r"$T$", fontsize=16, color=BLUE, ha="right",
            va="center")

    ax.set_xlim(-0.6, xa + 2.0)
    ax.set_ylim(-0.03 * t(root), t(root) * 1.08)
    ax.set_xticks([])
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.set_yticks([])
    ax.annotate("", xy=(-1.25, t(root) * 1.05), xytext=(-1.25, 0),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.4),
                annotation_clip=False)
    ax.text(-1.45, t(root) * 0.5, "time into the past", rotation=90,
            fontsize=12, color=INK2, ha="right", va="center")
    ax.set_facecolor("none")


def draw_q(ax, tree, edges, sample_T):
    t = tree.time
    tmax = t(tree.root)
    T = np.linspace(0, tmax, 600)
    for i, (child, parent) in enumerate(edges):
        q = derived_probability_uniform_edge_grid(
            tree, child, T, t(child), t(parent), TREE_NE)
        ax.plot(T, q, color=EDGE_COLORS[i], lw=2.6)
        ax.text(0.012 * tmax, q[0] + 0.025, "abc"[i], fontsize=15,
                style="italic", family="serif", fontweight="bold", color=INK,
                va="bottom")
    ax.axvline(sample_T, color=BLUE, lw=0.9, ls=":", alpha=0.7)
    ax.text(sample_T + 0.015 * tmax, 0.8, r"$T$", fontsize=15, color=BLUE,
            va="center")
    ax.set_xlim(0, tmax)
    ax.set_ylim(0, 0.85)
    kilo_ticks(ax, "x", tmax)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xlabel(r"sample age $T$ (generations)")
    ax.set_ylabel(r"$q_i(T)$")
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    style_axes(ax)


def read_example_posterior(replicate=EXAMPLE_REPLICATE):
    """Grid, density and run summary for the panel 3 example.

    Reads a small tracked TSV exported from the (gitignored) per-replicate
    inference outputs; its first line holds the summary as key=value pairs.
    """
    path = BATCH300 / f"insertion_results/example_posterior_replicate_{replicate:03d}.tsv"
    with open(path) as handle:
        header = handle.readline().lstrip("# ").split()
    meta = {k: float(v) for k, v in (item.split("=") for item in header)}
    data = np.loadtxt(path, skiprows=2)
    return data[:, 0], data[:, 1], meta


def draw_posterior(ax, replicate=EXAMPLE_REPLICATE):
    grid, dens, meta = read_example_posterior(replicate)
    dens = dens / dens.max()
    lo, hi = meta["ci95_lower_T"], meta["ci95_upper_T"]
    inside = (grid >= lo) & (grid <= hi)
    ax.fill_between(grid, dens, where=inside, color=BLUE, alpha=0.18, lw=0)
    ax.plot(grid, dens, color=BLUE, lw=2.6)
    true_T = meta["true_T"]
    ax.plot([true_T] * 2, [0, 1.13], color=INK, lw=1.5)
    ax.plot(true_T, 1.13, marker="v", color=INK, ms=8)
    ax.text(true_T + 120, 1.13, "true age", fontsize=12, color=INK2, va="center")
    ax.annotate("", xy=(lo, -0.06), xytext=(hi, -0.06),
                arrowprops=dict(arrowstyle="|-|", color=BLUE, lw=1.4,
                                mutation_scale=4), annotation_clip=False)
    ax.text(hi + 80, -0.06, "95%", fontsize=11, color=BLUE, va="center")
    # zoomed on the posterior so its width is visible
    ax.set_xlim(true_T - 1700, true_T + 1700)
    ticks = np.arange(4000, 7000, 1000)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{v / 1000:g}k" for v in ticks])
    ax.set_ylim(-0.17, 1.25)
    ax.set_yticks([])
    ax.spines["left"].set_bounds(0, 1.25)
    ax.set_xlabel(r"sample age $T$ (generations)")
    ax.set_ylabel(r"$P(T \mid \mathrm{data})$")
    style_axes(ax)
    return int(meta["sites_used"])


def draw_validation(ax):
    c = read_tsv(BATCH300 / "insertion_results/age_summary.tsv")
    v = read_tsv(VARNE)
    series = [(c, BLUE, "o", r"constant $N_e$"),
              (v, ORANGE, "^", r"variable $N_e(t)$")]
    stats = []
    for d, color, marker, label in series:
        ax.scatter(d["true_T"], d["map_T"], s=14, color=color, marker=marker,
                   alpha=0.45, lw=0, label=label)
        err = d["map_T"] - d["true_T"]
        stats.append((float(err.mean()),
                      float(np.polyfit(d["true_T"], d["map_T"], 1)[0])))
    ax.plot([0, 10_500], [0, 10_500], color=INK, lw=1.3, ls="--", zorder=0,
            label="1:1")
    ax.set_xlim(0, 10_500)
    ax.set_ylim(0, 13_500)
    ax.set_xticks(AGE_TICKS)
    ax.set_xticklabels(AGE_LABELS)
    ax.set_yticks(AGE_TICKS)
    ax.set_yticklabels(AGE_LABELS)
    ax.set_xlabel("true age (generations)")
    ax.set_ylabel("estimated age (MAP)")
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    leg = ax.legend(loc="upper left", fontsize=10, frameon=False,
                    handletextpad=0.3, markerscale=1.8, borderaxespad=0.1)
    for h in leg.legend_handles:
        h.set_alpha(0.9)
    style_axes(ax)
    return stats


HEADS = [("1", "Insert the lineage"),
         ("2", "One curve per site"),
         ("3", "Combine sites"),
         ("4", "Simulation check")]


def panel_texts(n_sites, stats):
    """(equation, equation fontsize, key) for each panel."""
    (cb, cs), (vb, vs) = stats
    return [
        (r"$S_T(u)=\exp\!\left[-\int_T^{\,u}\frac{\kappa(s)}{2N_e(s)}\,ds\right]$", 15,
         "$\\kappa$: tree lineages at time $s$\n"
         "$S_T(u)$: not yet joined by time $u$"),
        (r"$q_i(T)=\int_T^{\,z_i} S_T(u)\,\frac{\delta_i(u)}{2N_e(u)}\,du$", 15,
         "$\\delta_i$: lineages below mutation $i$\n"
         "$z_i$: mutation time (averaged over edge)"),
        (r"$L(T)=\prod_i q_i(T)^{x_i}\,\left[1-q_i(T)\right]^{1-x_i}$", 14,
         "$x_i = 1$: derived allele observed\n"
         f"one simulated sample, ~{n_sites / 1000:.0f}k sites,\n"
         "true modern ARG"),
        (r"$\hat T=\arg\max_T\,P(T\mid\mathrm{data})$", 14,
         "true ARG; mean $\\hat T - T$, slope\n"
         f"constant $N_e$: {cb:+.0f} gen, {cs:.2f}\n"
         f"variable $N_e(t)$: {vb:+.0f} gen, {vs:.2f}"),
    ]


def panel_frame(fig, x, y, w, h):
    fig.patches.append(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012",
        transform=fig.transFigure, facecolor=PANEL, edgecolor=GRID, lw=1,
        zorder=-10))


def panel_heading(fig, x, y, num, title):
    fig.text(x, y, num, fontsize=20, fontweight="bold", color=BLUE, family="serif")
    fig.text(x + 0.024, y, title, fontsize=16, fontweight="bold", color=INK,
             family="serif")


def draw_panels(axes, tree, edges, sample_T):
    draw_tree(axes[0], tree, edges, sample_T, edges[1])
    draw_q(axes[1], tree, edges, sample_T)
    n_sites = draw_posterior(axes[2])
    stats = draw_validation(axes[3])
    return n_sites, stats


def save(fig, out_stem, suffixes, dpi=220):
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in suffixes:
        path = out_stem.with_suffix(suffix)
        fig.savefig(path, facecolor="white", dpi=dpi if suffix == ".png" else None)
        print(path)


def build_slide(out_stem: Path):
    """16:9 slide: four panels in one row."""
    fig = plt.figure(figsize=(13.333, 7.5), facecolor="white")
    fig.text(0.03, 0.925, "Dating ancient genomes by threading them into a modern ARG",
             fontsize=26, fontweight="bold", color=INK, family="serif")

    tree = simulate_tree()
    edges = pick_edges(tree)
    sample_T = 0.35 * tree.time(edges[1][1])

    frames = [(0.020, 0.235), (0.268, 0.235), (0.516, 0.235), (0.764, 0.216)]
    y0, h = 0.04, 0.825
    for fx, fw in frames:
        panel_frame(fig, fx, y0, fw, h)
    for (ax0, aw0), (ax1, _) in zip(frames[:-1], frames[1:]):
        fig.patches.append(FancyArrowPatch(
            (ax0 + aw0 + 0.001, y0 + 0.62 * h), (ax1 - 0.001, y0 + 0.62 * h),
            transform=fig.transFigure, arrowstyle="-|>", mutation_scale=16,
            color=TREE, lw=2))
    for (fx, _), (num, title) in zip(frames, HEADS):
        panel_heading(fig, fx + 0.012, y0 + h - 0.05, num, title)

    plot_y, plot_h = 0.315, 0.455
    axes = [fig.add_axes([0.060, plot_y, 0.185, plot_h]),
            fig.add_axes([0.318, plot_y, 0.170, plot_h]),
            fig.add_axes([0.560, plot_y, 0.170, plot_h]),
            fig.add_axes([0.810, plot_y, 0.155, plot_h])]
    n_sites, stats = draw_panels(axes, tree, edges, sample_T)

    key = dict(fontsize=11, color=MUTED, va="top", linespacing=1.5)
    for (fx, _), (eq, size, text) in zip(frames, panel_texts(n_sites, stats)):
        fig.text(fx + 0.014, 0.215, eq, fontsize=size, color=INK, va="top")
        fig.text(fx + 0.014, 0.14, text, **key)
    save(fig, out_stem, (".png", ".pdf", ".svg"))


def build_readme(out_stem: Path):
    """README figure: 2 x 2 grid so text stays legible at ~900 px wide."""
    fig = plt.figure(figsize=(12, 10.6), facecolor="white")
    fig.text(0.025, 0.958, "Dating ancient genomes by threading them into a modern ARG",
             fontsize=22, fontweight="bold", color=INK, family="serif")

    tree = simulate_tree()
    edges = pick_edges(tree)
    sample_T = 0.35 * tree.time(edges[1][1])

    cw, ch = 0.47, 0.43
    cells = [(0.02, 0.495), (0.51, 0.495), (0.02, 0.05), (0.51, 0.05)]
    for (cx, cy), (num, title) in zip(cells, HEADS):
        panel_frame(fig, cx, cy, cw, ch)
        panel_heading(fig, cx + 0.015, cy + ch - 0.042, num, title)

    def plot_box(cx, cy, left):
        return [cx + left, cy + 0.215, cw - left - 0.04, ch - 0.29]

    axes = [fig.add_axes(plot_box(*cells[0], 0.07)),
            fig.add_axes(plot_box(*cells[1], 0.085)),
            fig.add_axes(plot_box(*cells[2], 0.085)),
            fig.add_axes(plot_box(*cells[3], 0.10))]
    n_sites, stats = draw_panels(axes, tree, edges, sample_T)

    key = dict(fontsize=11, color=MUTED, va="top", linespacing=1.5)
    for (cx, cy), (eq, size, text) in zip(cells, panel_texts(n_sites, stats)):
        fig.text(cx + 0.02, cy + 0.128, eq, fontsize=size, color=INK, va="center")
        fig.text(cx + 0.02, cy + 0.090, text, **key)

    fig.text(0.025, 0.012,
             "All panels: simulations with the true modern ARG. Open problems: nominal 95% "
             "intervals cover ~46% of true ages under linkage,\nand inferred ARGs give "
             "demography-dependent bias (issue #5).",
             fontsize=10.5, color=MUTED, linespacing=1.4)
    save(fig, out_stem, (".png",), dpi=160)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--layout", choices=("slide", "readme"), default="slide")
    p.add_argument("--out", type=Path, help="output path without suffix")
    a = p.parse_args()
    plt.rcParams.update({"font.family": "sans-serif", "svg.fonttype": "none",
                         "mathtext.fontset": "dejavuserif"})
    if a.layout == "slide":
        build_slide(a.out or ROOT / "docs/graphical_abstract_slide")
    else:
        build_readme(a.out or ROOT / "docs/graphical_abstract_readme")


if __name__ == "__main__":
    main()
