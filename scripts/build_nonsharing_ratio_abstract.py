#!/usr/bin/env python3
"""Graphical abstract for the non-sharing-ratio dating estimator.

Panels a-c are schematic; panel d plots simulation results (results/nonsharing_ratio/sim_inferred_vs_true.csv,
made by scripts/collect_nonsharing_ratio_sim_table.py).  Panels:
  a  an ancient lineage sampled at age T cannot coalesce during (0, T);
  b  panel haplotype carrying an old mutation (age > A0 >> T): the ancient
     target fails to share it slightly more often than a modern control;
  c  the ratio f_anc / mean f_mod = exp(int_0^T c) and the solved T.
  d  simulated ancient samples (SINGER+POLEGON ARGs, c(t) from the ARG): inferred vs true age,
     constant size and a bottleneck, A0 = 20,000.

Colours follow scripts/build_graphical_abstract.py (ancient = blue, extra
non-sharing = orange); the two are distinguishable under the common
colour-vision deficiencies, and non-sharing is also coded by marker shape.

Writes docs/nonsharing_ratio/graphical_abstract.{png,svg}.
Run through hpc_run:  ~/.claude/bin/hpc_run 'python scripts/build_nonsharing_ratio_abstract.py'
"""
from __future__ import annotations

from pathlib import Path

import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "nonsharing_ratio"

INK = "#17202B"
INK2 = "#3B4552"
MUTED = "#5A6472"
GRID = "#E1E5EA"
PANEL = "#F4F6F8"
TREE = "#8C96A3"
BLUE = "#2563A8"
BLUE_LT = "#DCE8F5"
ORANGE = "#D9731F"
ORANGE_DK = "#8F4710"

plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["font.family"] = "DejaVu Sans"


def panel_axes(fig, x, y, w, h):
    ax = fig.add_axes([x, y, w, h])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.05, 0.05), 9.9, 9.9, boxstyle="round,pad=0,rounding_size=0.25",
                                fc=PANEL, ec=GRID, lw=1.2, zorder=0))
    return ax


def heading(fig, x, y, letter, title):
    fig.text(x, y, letter, fontsize=20, fontweight="bold", color=BLUE, family="serif")
    fig.text(x + 0.03, y + 0.002, title, fontsize=15, fontweight="bold", color=INK, va="baseline")


def line(ax, xs, ys, color=INK2, lw=2.6, **kw):
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="round", **kw)


# ---------------------------------------------------------------- panel a
def panel_a(ax):
    y0, yT, ytop = 2.5, 5.2, 9.3
    ax.add_patch(Rectangle((0.9, y0), 8.8, yT - y0, fc=BLUE_LT, ec="none", alpha=0.7, zorder=1))
    xr, xm, xa = 3.0, 4.8, 7.6
    ycm = 3.9  # example modern r-m coalescence inside (0, T)
    # modern lineages
    line(ax, [xr, xr], [y0, ycm])
    line(ax, [xm, xm], [y0, ycm])
    line(ax, [xr, xm], [ycm, ycm])
    xrm = 0.5 * (xr + xm)
    # ancient lineage appears at T
    line(ax, [xa, xa], [yT, 6.9], color=BLUE)
    ax.plot([xa], [yT], "o", ms=11, color=BLUE, zorder=5)
    line(ax, [xrm, xrm], [ycm, 6.9])
    line(ax, [xrm, xa], [6.9, 6.9], color=INK2)
    xroot = 0.5 * (xrm + xa)
    line(ax, [xroot, xroot], [6.9, ytop])
    ax.plot([xr, xm], [y0, y0], "o", ms=11, color=INK2, zorder=5)
    # ghost of the ancient lineage inside (0, T)
    ax.plot([xa, xa], [y0, yT], ls=(0, (2, 3)), color=BLUE, lw=1.6, alpha=0.7)
    ax.text(xa - 0.3, 0.5 * (y0 + yT) - 0.2, "no ancient\nlineage here", fontsize=11, color=BLUE,
            va="center", ha="right", linespacing=1.2)
    ax.text(xr - 0.3, 2.3, "reference r", fontsize=11, color=INK2, ha="center", va="top")
    ax.text(xm + 0.15, 2.3, "modern\ncontrol m", fontsize=11, color=INK2, ha="center", va="top",
            linespacing=1.1)
    ax.text(xa + 0.3, yT + 0.1, "ancient\ntarget", fontsize=11, color=BLUE, ha="left", va="center",
            linespacing=1.2)
    ax.annotate("", xy=(1.0, ytop), xytext=(1.0, y0),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5))
    ax.text(0.7, 6.0, "time before present", rotation=90, fontsize=11, color=MUTED,
            ha="center", va="center")
    ax.text(1.25, y0 + 0.02, r"$0$", fontsize=14, color=INK, va="bottom")
    ax.text(1.25, yT + 0.1, r"$T$", fontsize=16, color=BLUE, va="bottom")
    ax.text(6.2, 8.2, "above $T$ every pair\nof lineages coalesces\nat rate $c(t)$",
            fontsize=11, color=INK2, ha="left", va="center", linespacing=1.3)
    ax.text(5.0, 0.75, "$r$ and $m$ may coalesce anywhere in $(0, T)$;\nthe ancient lineage does not exist until $T$",
            fontsize=11, color=INK, ha="center", va="center", linespacing=1.3, fontweight="bold")
    ax.text(5.0, 9.65 - 0.2, "", fontsize=1)


# ---------------------------------------------------------------- panel b
def panel_b(ax):
    # time strip
    ax.text(0.5, 9.3, "ARG says this mutation is old:", fontsize=12, color=INK, fontweight="bold",
            va="center")
    y = 8.2
    ax.annotate("", xy=(9.6, y), xytext=(0.6, y),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5))
    for xx, lab, col in ((1.0, r"$0$", INK), (2.6, r"$T$", BLUE), (5.9, r"$A_0$", ORANGE_DK)):
        ax.plot([xx, xx], [y - 0.18, y + 0.18], color=col, lw=2)
        ax.text(xx, y - 0.35, lab, fontsize=15, color=col, ha="center", va="top")
    ax.plot([8.0], [y], marker="*", ms=17, color=ORANGE, mec=ORANGE_DK, zorder=5)
    ax.text(8.0, y + 0.35, "mutation\nage $a$", fontsize=11, color=ORANGE_DK, ha="center",
            va="bottom", linespacing=1.1)
    ax.text(4.25, y + 0.35, r"$A_0 \gg T$", fontsize=12, color=INK2, ha="center", va="bottom")
    # genotype matrix
    ncol = 10
    x0, dx = 3.5, 0.6
    rows = [("panel haplotype r\n(derived allele)", 5.7, "ref"),
            ("modern control m", 4.35, "mod"),
            ("ancient target", 3.0, "anc")]
    miss = {"mod": {4}, "anc": {2, 7}}
    ax.text(x0 + 0.5 * (ncol - 1) * dx, 6.75, "old sites, one column each (schematic)",
            fontsize=11, color=MUTED, ha="center", va="center")
    for lab, yy, kind in rows:
        col = BLUE if kind == "anc" else INK2
        ax.text(x0 - 0.45, yy, lab, fontsize=11, color=col, ha="right", va="center", linespacing=1.1)
        for j in range(ncol):
            xx = x0 + j * dx
            if kind != "ref" and j in miss[kind]:
                ax.plot([xx], [yy], marker="X", ms=11, color=ORANGE, mec=ORANGE_DK, mew=0.8, zorder=5)
            else:
                ax.plot([xx], [yy], marker="o", ms=9, color=col, zorder=5)
    ax.plot([0.7], [1.55], marker="o", ms=8, color=INK2)
    ax.text(0.95, 1.55, "shares the derived allele", fontsize=11, color=INK2, va="center")
    ax.plot([5.1], [1.55], marker="X", ms=10, color=ORANGE, mec=ORANGE_DK, mew=0.8)
    ax.text(5.35, 1.55, "target is ancestral:\nnot shared", fontsize=11, color=ORANGE_DK,
            va="center", linespacing=1.1)
    ax.text(5.0, 0.55, "ancient target fails to share a little more often",
            fontsize=11, color=INK, ha="center", va="center", fontweight="bold")


# ---------------------------------------------------------------- panel c
def panel_c(ax):
    # cumulative hazard cartoon, axes in lower-left
    ox, oy, w, h = 1.2, 1.3, 3.3, 3.3
    ax.annotate("", xy=(ox + w + 0.3, oy), xytext=(ox, oy), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5))
    ax.annotate("", xy=(ox, oy + h + 0.3), xytext=(ox, oy), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5))
    t = np.linspace(0, 1, 100)
    Hc = 1 - np.exp(-2.2 * t) * (1 + 0.0 * t)  # a concave cumulative hazard (cartoon)
    Hc = Hc / Hc[-1]
    ax.plot(ox + w * t, oy + h * Hc * 0.95, color=INK, lw=2.8)
    Ts = 0.5
    Hs = (1 - np.exp(-2.2 * Ts)) / (1 - np.exp(-2.2)) * 0.95
    xT, yH = ox + w * Ts, oy + h * Hs
    ax.plot([ox, xT], [yH, yH], color=ORANGE, lw=2.2, ls=(0, (4, 3)))
    ax.plot([xT, xT], [oy, yH], color=BLUE, lw=2.2, ls=(0, (4, 3)))
    ax.plot([xT], [yH], "o", ms=9, color=INK, zorder=5)
    ax.text(ox - 0.12, yH, r"$\hat H$", fontsize=16, color=ORANGE_DK, ha="right", va="center")
    ax.text(xT, oy - 0.15, r"$\hat T$", fontsize=16, color=BLUE, ha="center", va="top")
    ax.text(ox + w + 0.3, oy - 0.15, "age", fontsize=11, color=MUTED, ha="right", va="top")
    ax.text(ox + 0.15, oy + h + 0.7, r"$\int_0^{t} c(u)\,du$", fontsize=14, color=INK, ha="left", va="center")
    ax.text(0.5, 9.2, "$c(t)$: ARG lineage counts\nor a demes file", fontsize=11,
            color=MUTED, ha="left", va="center", linespacing=1.2)
    # equations on the right
    xe = 5.45
    ax.text(xe, 9.35, "each old mutation, $a > T$:", fontsize=11.5, color=MUTED, va="center")
    ax.text(xe, 8.3, r"$\dfrac{P_{\rm anc}}{P_{\rm mod}}=e^{\int_0^T c(t)\,dt}$", fontsize=19, color=INK, va="center")
    ax.text(xe, 6.85, "estimate from data:", fontsize=11.5, color=MUTED, va="center")
    ax.text(xe, 5.75, r"$f=\dfrac{\#\,\mathrm{not\ shared}}{\#\,\mathrm{pairs}}$",
            fontsize=16, color=INK, va="center")
    ax.text(xe, 3.9, r"$\hat H=\log\dfrac{f_{\rm anc}}{\bar{f}_{\rm mod}}$", fontsize=17, color=INK, va="center")
    ax.text(xe, 2.3, r"$\int_0^{\hat T} c(t)\,dt=\hat H$", fontsize=18, color=INK, va="center")
    ax.text(xe, 0.85, "pairs: (site, reference) with\nthe reference derived.\n$P$: chance a pair is not shared", fontsize=10,
            color=MUTED, va="center", linespacing=1.2)


SIM_CSV = ROOT / "results" / "nonsharing_ratio" / "sim_inferred_vs_true.csv"


def panel_d(fig, rect, demography, title, color, a0="20000"):
    rows = [r for r in csv.DictReader(open(SIM_CSV)) if r["demography"] == demography]
    x = np.array([float(r["true_age"]) for r in rows])
    y = np.array([float(r["T_arg_" + a0]) for r in rows])
    ax = fig.add_axes(rect)
    ax.plot([0, 6500], [0, 6500], color=MUTED, lw=1.3, ls="--", zorder=1)
    ax.scatter(x, y, s=34, color=color, alpha=0.85, edgecolor="white", lw=0.6, zorder=2)
    r = np.corrcoef(x, y)[0, 1]
    ax.text(0.04, 0.95, f"r = {r:.2f}\nmean bias {np.mean(y - x):+.0f} gen\nn = {len(x)} ancient samples",
            transform=ax.transAxes, va="top", fontsize=11.5, color=INK)
    ax.set_xlim(0, 6500)
    ax.set_ylim(-300, 6500)
    ax.set_xlabel("true age (generations)", fontsize=11.5, color=INK2)
    ax.set_ylabel("inferred age (generations)", fontsize=11.5, color=INK2)
    ax.set_title(title, fontsize=13, fontweight="bold", color=INK, loc="left")
    ax.tick_params(labelsize=10.5, colors=INK2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    HT, HB = 6.9, 5.1                      # inches: schematic row, data row
    FH = HT + HB
    k = HT / FH
    ty = lambda y: (1 - k) + y * k         # map the original 0-1 top-row coordinates into the taller figure
    fig = plt.figure(figsize=(18, FH), facecolor="white")
    fig.text(0.5, ty(0.955), "Dating ancient genomes by the non-sharing ratio", ha="center", fontsize=24,
             fontweight="bold", color=INK, family="serif")
    fig.text(0.5, ty(0.893), "An ancient lineage cannot coalesce with modern ones for its first $T$ generations; "
             "that deficit is the clock", ha="center", fontsize=13.5, color=INK2)
    W, H, y0 = 0.30, 0.70 * k, ty(0.075)
    xs = (0.02, 0.35, 0.68)
    for x, l, t in zip(xs, "abc", ("Ancient lineage starts late", "Old mutations, who shares", "Ratio gives age")):
        heading(fig, x, y0 + H + 0.012 * k, l, t)
    panel_a(panel_axes(fig, xs[0], y0, W, H))
    panel_b(panel_axes(fig, xs[1], y0, W, H))
    panel_c(panel_axes(fig, xs[2], y0, W, H))

    # row 2: simulation results
    fig.text(0.02, 0.385, "d", fontsize=20, fontweight="bold", color=BLUE, family="serif")
    fig.text(0.05, 0.387, "Works on simulated ancient samples, with ARGs and c(t) both inferred", fontsize=15,
             fontweight="bold", color=INK, va="baseline")
    yb, hb = 0.105, 0.225
    panel_d(fig, [0.06, yb, 0.22, hb], "constant", "Constant size", BLUE)
    panel_d(fig, [0.35, yb, 0.22, hb], "bottleneck", "Bottleneck (0.3 Ne, 500 to 300 gens ago)", ORANGE)
    fig.text(0.62, 0.325, "What was simulated", fontsize=12.5, fontweight="bold", color=INK, va="top")
    fig.text(0.62, 0.295,
             "msprime, 20 Mb, 5 replicates each; 26 panel + 10 held-out\n"
             "modern controls; 50 ancient samples per demography\n"
             "(true ages 0-6,000 gen, not in the ARG).\n"
             "ARG inferred by SINGER + POLEGON (5 draws, averaged);\n"
             "c(t) from the inferred ARG's lineage counts;\n"
             "A0 = 20,000 gen on the child-node time of the site.\n"
             "Using the true demography for c(t) changes mean bias\n"
             "by about 190 gen. Dashed line: inferred = true.",
             fontsize=11.5, color=INK2, va="top", linespacing=1.45)
    fig.text(0.5, 0.022, "a-c schematic; d simulated data (results/nonsharing_ratio/sim_inferred_vs_true.csv). "
             "A first real-data check did not recover ages; see docs/nonsharing_ratio/REAL_DATA_OFFSET_REPORT.md.",
             ha="center", fontsize=11, color=MUTED)
    fig.savefig(OUT / "graphical_abstract.png", dpi=160, facecolor="white")
    fig.savefig(OUT / "graphical_abstract.svg", facecolor="white")
    print("wrote", OUT / "graphical_abstract.png")


if __name__ == "__main__":
    main()
