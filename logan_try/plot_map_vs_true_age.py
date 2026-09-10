#!/usr/bin/env python3
"""Plot inferred MAP ages against supplied true/reference ages as SVG."""

import argparse
import csv
import html
import math
import re
from pathlib import Path


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def stats(rows):
    x = [row["reference_age"] for row in rows]
    y = [row["map_T"] for row in rows]
    xm, ym = sum(x) / len(x), sum(y) / len(y)
    cov = sum((a - xm) * (b - ym) for a, b in zip(x, y))
    ssx = sum((a - xm) ** 2 for a in x)
    ssy = sum((b - ym) ** 2 for b in y)
    pearson = cov / math.sqrt(ssx * ssy) if ssx and ssy else math.nan
    rmse = math.sqrt(sum((b - a) ** 2 for a, b in zip(x, y)) / len(x))
    return pearson, rmse


parser = argparse.ArgumentParser()
parser.add_argument("run_dir", type=Path)
parser.add_argument("ages_tsv", type=Path)
parser.add_argument("--epsilon", type=float, required=True)
parser.add_argument("--model-label", default="diffusion")
args = parser.parse_args()

estimates = read_tsv(args.run_dir / "ages_table.tsv")
ages = {row["sample_id"]: row for row in read_tsv(args.ages_tsv)}
missing = [row["sample"] for row in estimates if row["sample"] not in ages]
if missing:
    raise SystemExit(f"Missing true ages for: {', '.join(missing)}")

rows = []
for estimate in estimates:
    age = ages[estimate["sample"]]
    cal_limits = [float(value) for value in re.findall(
        r"\d+(?:\.\d+)?", age["cal_bp_95"])]
    raw_age = float(age["true_age"])
    reference_age = ((min(cal_limits) + max(cal_limits)) / 2.0
                     if len(cal_limits) >= 2 else raw_age)
    rows.append(
        {
            "sample": estimate["sample"],
            "sample_type": age["sample_type"],
            "radiocarbon_age_bp": raw_age,
            "reference_age": reference_age,
            "age_source": age["age_source"],
            "cal_bp_95": age["cal_bp_95"],
            "radiocarbon_error_1sigma": (
                float(age["radiocarbon_error_1sigma"])
                if age["radiocarbon_error_1sigma"] else math.nan
            ),
            "map_T": float(estimate["map_T"]),
            "ci95_lower_T": float(estimate["ci95_lower_T"]),
            "ci95_upper_T": float(estimate["ci95_upper_T"]),
        }
    )

joined = args.run_dir / "map_vs_true_age.tsv"
with joined.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

all_r, all_rmse = stats(rows)
ancient = [row for row in rows if row["sample_type"] == "ancient"]
anc_r, anc_rmse = stats(ancient)

width, height = 920, 800
left, right, top, bottom = 105, 45, 85, 120
plot_w, plot_h = width - left - right, height - top - bottom
axis_max = max(2000, int(math.ceil(max(max(r["reference_age"], r["map_T"], r["ci95_upper_T"]) for r in rows) / 500) * 500))


def px(value):
    return left + value / axis_max * plot_w


def py(value):
    return top + plot_h - value / axis_max * plot_h


parts = [
    '<?xml version="1.0" encoding="UTF-8"?>',
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
    "<style>",
    "text { font-family: Arial, Helvetica, sans-serif; fill: #20242a; }",
    ".grid { stroke: #dfe3e8; stroke-width: 1; }",
    ".axis { stroke: #20242a; stroke-width: 1.5; }",
    ".ci { stroke: #087e8b; stroke-width: 1.2; opacity: 0.28; }",
    ".age95 { stroke: #087e8b; stroke-width: 1.8; opacity: 0.45; }",
    ".ancient { fill: #087e8b; stroke: white; stroke-width: 1.2; opacity: 0.9; }",
    ".modern { fill: #d95f02; stroke: white; stroke-width: 1.2; opacity: 0.95; }",
    "</style>",
    '<rect width="100%" height="100%" fill="white"/>',
    f'<text x="{width/2}" y="32" text-anchor="middle" font-size="22" font-weight="bold">Genome-wide MAP age vs. reference age</text>',
    f'<text x="{width/2}" y="57" text-anchor="middle" font-size="13" fill="#59636e">18 ancient + 2 modern; 10 chromosomes; 10 ARG draws; {html.escape(args.model_label)}; epsilon = {args.epsilon:g}</text>',
]

tick_step = 500 if axis_max <= 3000 else 1000
for tick in range(0, axis_max + 1, tick_step):
    xx, yy = px(tick), py(tick)
    parts.extend([
        f'<line class="grid" x1="{xx:.2f}" y1="{top}" x2="{xx:.2f}" y2="{top+plot_h}"/>',
        f'<line class="grid" x1="{left}" y1="{yy:.2f}" x2="{left+plot_w}" y2="{yy:.2f}"/>',
        f'<text x="{xx:.2f}" y="{top+plot_h+27}" text-anchor="middle" font-size="12">{tick}</text>',
        f'<text x="{left-14}" y="{yy+4:.2f}" text-anchor="end" font-size="12">{tick}</text>',
    ])

parts.extend([
    f'<line class="axis" x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}"/>',
    f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}"/>',
    f'<line x1="{px(0):.2f}" y1="{py(0):.2f}" x2="{px(axis_max):.2f}" y2="{py(axis_max):.2f}" stroke="#6b7280" stroke-width="1.5" stroke-dasharray="7 6"/>',
    f'<text x="{px(axis_max*0.78):.2f}" y="{py(axis_max*0.82):.2f}" font-size="12" fill="#6b7280">1:1</text>',
])

for row in rows:
    xx, yy = px(row["reference_age"]), py(row["map_T"])
    lo, hi = py(row["ci95_lower_T"]), py(row["ci95_upper_T"])
    parts.append(f'<line class="ci" x1="{xx:.2f}" y1="{hi:.2f}" x2="{xx:.2f}" y2="{lo:.2f}"/>')
    cal_limits = [float(value) for value in re.findall(r"\d+(?:\.\d+)?", row["cal_bp_95"])]
    if len(cal_limits) >= 2:
        xlo, xhi = px(min(cal_limits)), px(max(cal_limits))
        parts.extend([
            f'<line class="age95" x1="{xlo:.2f}" y1="{yy:.2f}" x2="{xhi:.2f}" y2="{yy:.2f}"/>',
            f'<line class="age95" x1="{xlo:.2f}" y1="{yy-4:.2f}" x2="{xlo:.2f}" y2="{yy+4:.2f}"/>',
            f'<line class="age95" x1="{xhi:.2f}" y1="{yy-4:.2f}" x2="{xhi:.2f}" y2="{yy+4:.2f}"/>',
        ])
    label = html.escape(
        f'{row["sample"]}: reference={row["reference_age"]:.0f}, MAP={row["map_T"]:.0f}, '
        f'95% CI={row["ci95_lower_T"]:.0f}–{row["ci95_upper_T"]:.0f}'
    )
    parts.append(f'<circle class="{row["sample_type"]}" cx="{xx:.2f}" cy="{yy:.2f}" r="5.5"><title>{label}</title></circle>')

parts.extend([
    f'<text x="{left+plot_w/2}" y="{height-50}" text-anchor="middle" font-size="16">Reference age (cal-BP midpoint where available; modern = 50)</text>',
    f'<text x="28" y="{top+plot_h/2}" text-anchor="middle" font-size="16" transform="rotate(-90 28 {top+plot_h/2})">Inferred MAP age (years; 1 year/generation)</text>',
    f'<rect x="{left+18}" y="{top+18}" width="250" height="64" rx="4" fill="white" stroke="#c7cdd4" opacity="0.94"/>',
    f'<text x="{left+31}" y="{top+41}" font-size="13">All: r = {all_r:.2f}; RMSE = {all_rmse:.0f}</text>',
    f'<text x="{left+31}" y="{top+63}" font-size="13">Ancient only: r = {anc_r:.2f}; RMSE = {anc_rmse:.0f}</text>',
    f'<circle class="ancient" cx="{left+plot_w-190}" cy="{top+28}" r="5"/><text x="{left+plot_w-177}" y="{top+33}" font-size="12">Ancient (n={len(ancient)})</text>',
    f'<circle class="modern" cx="{left+plot_w-190}" cy="{top+50}" r="5"/><text x="{left+plot_w-177}" y="{top+55}" font-size="12">Modern (n={len(rows)-len(ancient)})</text>',
    f'<text x="{width/2}" y="{height-18}" text-anchor="middle" font-size="11" fill="#68727d">Vertical lines: inferred 95% intervals. Horizontal bars: calibrated age 95% intervals where reported.</text>',
    "</svg>",
])

svg = args.run_dir / "map_vs_true_age.svg"
svg.write_text("\n".join(parts) + "\n")
print(f"wrote {svg}")
print(f"wrote {joined}")
print(f"all: r={all_r:.4f}, RMSE={all_rmse:.1f}, n={len(rows)}")
print(f"ancient: r={anc_r:.4f}, RMSE={anc_rmse:.1f}, n={len(ancient)}")
