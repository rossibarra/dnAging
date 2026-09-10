#!/usr/bin/env python3
"""Select a reproducible 18-ancient/2-modern cohort from the 313-sample VCF."""

import csv
import random
import re
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parent
WORKBOOK = ROOT.parent / ".." / "original_files" / "meta_v10_9July2026.xlsx"
FAM = ROOT.parent / ".." / "original_files" / "430goodMays_no_cov_bias_G2.mind95.fam"
ORIGINAL = ROOT / "radiocarbon_selection.tsv"
OUT = ROOT / "cohort_eps006_18anc_2modern"
SEED = 20260908
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def workbook_rows(path):
    with ZipFile(path) as archive:
        strings = []
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        for item in root.findall(NS + "si"):
            strings.append("".join(node.text or "" for node in item.iter(NS + "t")))
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows = []
    for row in sheet.iter(NS + "row"):
        row_number = int(row.attrib["r"])
        values = {}
        for cell in row.findall(NS + "c"):
            column = re.match(r"[A-Z]+", cell.attrib["r"]).group()
            value_node = cell.find(NS + "v")
            value = "" if value_node is None else value_node.text or ""
            if cell.attrib.get("t") == "s" and value:
                value = strings[int(value)]
            elif cell.attrib.get("t") == "inlineStr":
                value = "".join(node.text or "" for node in cell.iter(NS + "t"))
            values[column] = value.strip()
        if row_number >= 4:
            rows.append((row_number, values))
    return rows


fam_rows = [line.split() for line in FAM.open()]
sample_by_alias = {}
for fields in fam_rows:
    sample = f"{fields[0]}_{fields[1]}"
    sample_by_alias.setdefault(fields[0], sample)
    sample_by_alias.setdefault(fields[1], sample)

rows = workbook_rows(WORKBOOK)
by_alias = {}
for row_number, values in rows:
    alias = values.get("B", "")
    if alias:
        by_alias.setdefault(alias, []).append((row_number, values))

with ORIGINAL.open(newline="") as handle:
    ancient_pool = {row["sample_id"]: dict(row)
                    for row in csv.DictReader(handle, delimiter="\t")}

# These are the two unambiguous dated individuals deliberately omitted from the
# original 50 only because their 400-BP ages were redundant for age coverage.
for alias in ("aBM_H_2", "aPM_L_177"):
    sample = sample_by_alias[alias]
    hits = by_alias[alias]
    ages = {float(values["M"]) for _, values in hits}
    if len(ages) != 1:
        raise SystemExit(f"{alias} does not have one unambiguous radiocarbon age")
    row_number, values = hits[0]
    ancient_pool[sample] = {
        "sample_id": sample,
        "radiocarbon_age_bp": f"{next(iter(ages)):g}",
        "radiocarbon_error_1sigma": values.get("N", ""),
        "cal_bp_95": values.get("S", ""),
        "calibration": values.get("U", ""),
        "workbook_row": str(row_number),
        "publication_name": values.get("A", ""),
        "fam_alias": alias,
    }

if len(ancient_pool) != 52:
    raise SystemExit(f"expected 52 unambiguous dated ancient samples; got {len(ancient_pool)}")

modern_pool = []
for alias, hits in by_alias.items():
    sample = sample_by_alias.get(alias)
    if sample is None:
        continue
    if any(values.get("M", "").lower() == "modern" for _, values in hits):
        modern_pool.append((sample, alias))
modern_pool = sorted(set(modern_pool))
if len(modern_pool) < 2:
    raise SystemExit("fewer than two metadata-confirmed modern VCF samples")

rng = random.Random(SEED)
ancient = rng.sample(sorted(ancient_pool), 18)
modern = rng.sample(modern_pool, 2)
cohort = [(sample, "ancient") for sample in ancient]
cohort += [(sample, "modern") for sample, _alias in modern]
rng.shuffle(cohort)

OUT.mkdir(parents=True, exist_ok=False)
(OUT / "samples.txt").write_text("\n".join(sample for sample, _ in cohort) + "\n")
with (OUT / "sample_ages.tsv").open("w", newline="") as handle:
    fields = ["sample_id", "sample_type", "true_age", "age_source",
              "radiocarbon_error_1sigma", "cal_bp_95", "fam_alias"]
    writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    modern_alias = dict(modern)
    for sample, sample_type in cohort:
        if sample_type == "modern":
            writer.writerow({
                "sample_id": sample, "sample_type": sample_type,
                "true_age": "50", "age_source": "assumed modern age",
                "radiocarbon_error_1sigma": "", "cal_bp_95": "",
                "fam_alias": modern_alias[sample],
            })
        else:
            row = ancient_pool[sample]
            writer.writerow({
                "sample_id": sample, "sample_type": sample_type,
                "true_age": row["radiocarbon_age_bp"],
                "age_source": "conventional radiocarbon age BP",
                "radiocarbon_error_1sigma": row["radiocarbon_error_1sigma"],
                "cal_bp_95": row["cal_bp_95"], "fam_alias": row["fam_alias"],
            })
(OUT / "selection.json").write_text(
    "{\n"
    f'  "random_seed": {SEED},\n'
    f'  "ancient_pool_size": {len(ancient_pool)},\n'
    f'  "modern_pool_size": {len(modern_pool)},\n'
    '  "n_ancient": 18,\n'
    '  "n_modern": 2\n'
    "}\n"
)
ages = [float(ancient_pool[sample]["radiocarbon_age_bp"]) for sample in ancient]
print(f"selected 18 ancient ages {min(ages):g}..{max(ages):g} BP and 2 modern at age 50")
print(f"seed={SEED}; wrote {OUT}")
