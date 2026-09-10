#!/usr/bin/env python3
"""Plot MAP sample-age estimates against radiocarbon ages."""

import argparse
import csv
import html
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RADIOCARBON = ROOT / "radiocarbon_selection.tsv"
parser = argparse.ArgumentParser()
parser.add_argument(
    "run_dir", type=Path, nargs="?",
    default=ROOT / "results" / "smoke_restart" / "1")
parser.add_argument("--label", default="Chromosome 1")
args = parser.parse_args()
OUT_DIR = args.run_dir.resolve()
ESTIMATES = OUT_DIR / "ages_table.tsv"
JOINED = OUT_DIR / "map_vs_radiocarbon.tsv"
SVG = OUT_DIR / "map_vs_radiocarbon.svg"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


rc_rows = {row["sample_id"]: row for row in read_tsv(RADIOCARBON)}
estimate_rows = read_tsv(ESTIMATES)
missing = [row["sample"] for row in estimate_rows if row["sample"] not in rc_rows]
if missing:
    raise SystemExit(f"Missing radiocarbon ages for: {', '.join(missing)}")

rows = []
for estimate in estimate_rows:
    rc = rc_rows[estimate["sample"]]
    rows.append(
        {
            "sample": estimate["sample"],
            "radiocarbon_age_bp": float(rc["radiocarbon_age_bp"]),
            "radiocarbon_error_1sigma": (
                float(rc["radiocarbon_error_1sigma"])
                if rc["radiocarbon_error_1sigma"]
                else math.nan
            ),
            "map_T": float(estimate["map_T"]),
            "ci95_lower_T": float(estimate["ci95_lower_T"]),
            "ci95_upper_T": float(estimate["ci95_upper_T"]),
        }
    )

with JOINED.open("w", newline="") as handle:
    fields = list(rows[0])
    writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

x = [row["radiocarbon_age_bp"] for row in rows]
y = [row["map_T"] for row in rows]
x_mean = sum(x) / len(x)
y_mean = sum(y) / len(y)
cov = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y))
ssx = sum((a - x_mean) ** 2 for a in x)
ssy = sum((b - y_mean) ** 2 for b in y)
pearson = cov / math.sqrt(ssx * ssy)
rmse = math.sqrt(sum((b - a) ** 2 for a, b in zip(x, y)) / len(x))

width, height = 920, 800
left, right, top, bottom = 105, 45, 75, 115
plot_w, plot_h = width - left - right, height - top - bottom
axis_max = 6000.0


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
    ".ci { stroke: #3b82a0; stroke-width: 1; opacity: 0.22; }",
    ".rcerr { stroke: #3b82a0; stroke-width: 1; opacity: 0.35; }",
    ".point { fill: #087e8b; stroke: white; stroke-width: 1.2; opacity: 0.88; }",
    ".boundary { fill: #d95f02; stroke: white; stroke-width: 1.2; }",
    "</style>",
    '<rect width="100%" height="100%" fill="white"/>',
    f'<text x="{width/2}" y="32" text-anchor="middle" font-size="22" font-weight="bold">{html.escape(args.label)} MAP age vs. radiocarbon age</text>',
    f'<text x="{width/2}" y="55" text-anchor="middle" font-size="13" fill="#59636e">50 ancient maize samples; one ARG draw; 1 year per generation</text>',
]

for tick in range(0, 6001, 1000):
    xx, yy = px(tick), py(tick)
    parts.extend(
        [
            f'<line class="grid" x1="{xx:.2f}" y1="{top}" x2="{xx:.2f}" y2="{top+plot_h}"/>',
            f'<line class="grid" x1="{left}" y1="{yy:.2f}" x2="{left+plot_w}" y2="{yy:.2f}"/>',
            f'<text x="{xx:.2f}" y="{top+plot_h+27}" text-anchor="middle" font-size="12">{tick}</text>',
            f'<text x="{left-14}" y="{yy+4:.2f}" text-anchor="end" font-size="12">{tick}</text>',
        ]
    )

parts.extend(
    [
        f'<line class="axis" x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}"/>',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}"/>',
        f'<line x1="{px(0):.2f}" y1="{py(0):.2f}" x2="{px(axis_max):.2f}" y2="{py(axis_max):.2f}" stroke="#6b7280" stroke-width="1.5" stroke-dasharray="7 6"/>',
        f'<text x="{px(4550):.2f}" y="{py(4750):.2f}" font-size="12" fill="#6b7280">1:1</text>',
    ]
)

for row in rows:
    xx, yy = px(row["radiocarbon_age_bp"]), py(row["map_T"])
    lo, hi = py(row["ci95_lower_T"]), py(row["ci95_upper_T"])
    parts.append(f'<line class="ci" x1="{xx:.2f}" y1="{hi:.2f}" x2="{xx:.2f}" y2="{lo:.2f}"/>')
    if math.isfinite(row["radiocarbon_error_1sigma"]):
        xlo = px(max(0, row["radiocarbon_age_bp"] - row["radiocarbon_error_1sigma"]))
        xhi = px(row["radiocarbon_age_bp"] + row["radiocarbon_error_1sigma"])
        parts.append(f'<line class="rcerr" x1="{xlo:.2f}" y1="{yy:.2f}" x2="{xhi:.2f}" y2="{yy:.2f}"/>')
    klass = "boundary" if row["map_T"] == 0 else "point"
    label = html.escape(
        f'{row["sample"]}: 14C={row["radiocarbon_age_bp"]:.0f} BP, '
        f'MAP={row["map_T"]:.0f}, 95% CI={row["ci95_lower_T"]:.0f}–{row["ci95_upper_T"]:.0f}'
    )
    parts.append(f'<circle class="{klass}" cx="{xx:.2f}" cy="{yy:.2f}" r="5.2"><title>{label}</title></circle>')

parts.extend(
    [
        f'<text x="{left+plot_w/2}" y="{height-47}" text-anchor="middle" font-size="16">Conventional radiocarbon age (14C years BP)</text>',
        f'<text x="28" y="{top+plot_h/2}" text-anchor="middle" font-size="16" transform="rotate(-90 28 {top+plot_h/2})">Inferred MAP age (generations; 1 year/generation)</text>',
        f'<rect x="{left+18}" y="{top+18}" width="224" height="63" rx="4" fill="white" stroke="#c7cdd4" opacity="0.94"/>',
        f'<text x="{left+31}" y="{top+42}" font-size="13">Pearson r = {pearson:.2f}</text>',
        f'<text x="{left+31}" y="{top+63}" font-size="13">RMSE = {rmse:.0f} years</text>',
        f'<circle class="boundary" cx="{left+plot_w-174}" cy="{top+28}" r="5"/><text x="{left+plot_w-160}" y="{top+33}" font-size="12">MAP at T=0 (n={sum(v == 0 for v in y)})</text>',
        f'<text x="{width/2}" y="{height-18}" text-anchor="middle" font-size="11" fill="#68727d">Vertical lines: inferred 95% intervals. Horizontal lines: radiocarbon ±1σ where reported. Radiocarbon ages are not calendar-calibrated.</text>',
        "</svg>",
    ]
)

SVG.write_text("\n".join(parts) + "\n")
print(f"wrote {SVG}")
print(f"wrote {JOINED}")
print(f"Pearson r={pearson:.4f}; RMSE={rmse:.1f}; n={len(rows)}")
