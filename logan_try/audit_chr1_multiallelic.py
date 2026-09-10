#!/usr/bin/env python3
"""Count chr1 sites recoverable from multiallelic panel VCF records."""

import json
import sys
from collections import defaultdict


panel_path, ancient_path, output_path = sys.argv[1:]
chrom_wanted = "1"
bases = {"A", "C", "G", "T"}

# Keep biallelic ancient allele pairs by position. This directly addresses the
# case in question: ancient REF/ALT is one pair within a multiallelic panel row.
ancient_pairs = defaultdict(set)
ancient_multiallelic_positions = set()
ancient_records = 0
with open(ancient_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t", 9)
        if fields[0] != chrom_wanted:
            continue
        ancient_records += 1
        pos = int(fields[1])
        ref = fields[3].upper()
        alts = tuple(fields[4].upper().split(","))
        if len(alts) != 1:
            ancient_multiallelic_positions.add(pos)
            continue
        alt = alts[0]
        if ref in bases and alt in bases and ref != alt:
            ancient_pairs[pos].add((ref, alt))


def called_counts(fields, ref, alts, ancient_ref, ancient_alt):
    """Return target-ALT and target-pair called counts; other ALTs are missing."""
    fmt = fields[8].split(":")
    if "GT" not in fmt:
        return 0, 0
    gt_i = fmt.index("GT")
    allele_by_index = (ref,) + alts
    target_alt = called = 0
    for sample in fields[9:]:
        values = sample.split(":")
        gt = values[gt_i] if gt_i < len(values) else "."
        for token in gt.replace("|", "/").split("/"):
            if token == ".":
                continue
            try:
                allele = allele_by_index[int(token)]
            except (ValueError, IndexError):
                continue
            if allele == ancient_ref:
                called += 1
            elif allele == ancient_alt:
                called += 1
                target_alt += 1
    return target_alt, called


baseline = set()
panel_multi_contains_pair = set()
panel_multi_same_ref_alt_present = set()
panel_multi_swapped_pair = set()
panel_multi_min20 = set()
panel_multi_min20_polymorphic = set()
same_position = set()
panel_records = panel_multiallelic_records = 0

with open(panel_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields[0] != chrom_wanted:
            continue
        panel_records += 1
        pos = int(fields[1])
        pairs = ancient_pairs.get(pos)
        if not pairs:
            continue
        same_position.add(pos)
        ref = fields[3].upper()
        alts = tuple(fields[4].upper().split(","))
        allele_set = {ref, *alts}
        if len(alts) == 1:
            panel_pair = frozenset((ref, alts[0]))
            if any(frozenset(pair) == panel_pair for pair in pairs):
                baseline.add(pos)
            continue

        panel_multiallelic_records += 1
        for ancient_ref, ancient_alt in pairs:
            if {ancient_ref, ancient_alt}.issubset(allele_set):
                panel_multi_contains_pair.add(pos)
                if ref == ancient_ref and ancient_alt in alts:
                    panel_multi_same_ref_alt_present.add(pos)
                if ref == ancient_alt and ancient_ref in alts:
                    panel_multi_swapped_pair.add(pos)
                alt_count, called = called_counts(
                    fields, ref, alts, ancient_ref, ancient_alt
                )
                if called >= 20:
                    panel_multi_min20.add(pos)
                    if 0 < alt_count < called:
                        panel_multi_min20_polymorphic.add(pos)

result = {
    "chrom": chrom_wanted,
    "ancient_records": ancient_records,
    "ancient_biallelic_positions": len(ancient_pairs),
    "ancient_multiallelic_positions_not_assessed": len(ancient_multiallelic_positions),
    "panel_records": panel_records,
    "same_position_with_biallelic_ancient": len(same_position),
    "current_biallelic_compatible_positions": len(baseline),
    "recoverable_panel_multiallelic_positions": len(panel_multi_contains_pair),
    "recoverable_same_ref_and_alt_present": len(panel_multi_same_ref_alt_present),
    "recoverable_ref_alt_swapped": len(panel_multi_swapped_pair),
    "recoverable_after_target_pair_min_n_20": len(panel_multi_min20),
    "recoverable_after_min_n_20_and_polymorphic": len(panel_multi_min20_polymorphic),
    "candidate_total_current_plus_recoverable": len(baseline | panel_multi_contains_pair),
    "candidate_gain_percent": (
        100.0 * len(panel_multi_contains_pair - baseline) / len(baseline)
        if baseline else None
    ),
}
with open(output_path, "w") as handle:
    json.dump(result, handle, indent=2)
    handle.write("\n")
print(json.dumps(result, indent=2))
