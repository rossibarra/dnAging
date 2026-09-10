import collections
import sys


panel_path, ancient_path = sys.argv[1:]
offsets = {
    "1": 0,
    "2": 308452472,
    "3": 552127664,
    "4": 790145432,
    "5": 1040475893,
    "6": 1266829343,
    "7": 1448186578,
    "8": 1633995495,
    "9": 1816406698,
    "10": 1979411443,
}
base = {"A": 0, "C": 1, "G": 2, "T": 3}
deltas = (-2, -1, 0, 1, 2)


def decode(fields):
    chrom = fields[0]
    if chrom not in offsets:
        return None
    ref = base.get(fields[3])
    alt = base.get(fields[4])
    if ref is None or alt is None:
        return None
    position = offsets[chrom] + int(fields[1])
    low, high = sorted((ref, alt))
    allele_code = low * 4 + high
    return chrom, position, allele_code


ancient = set()
with open(ancient_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        item = decode(line.split("\t", 5))
        if item is not None:
            _, position, allele_code = item
            ancient.add(position * 16 + allele_code)

matches = {delta: collections.Counter() for delta in deltas}
strand_only = {delta: collections.Counter() for delta in deltas}
panel_records = collections.Counter()
with open(panel_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        item = decode(line.split("\t", 5))
        if item is None:
            continue
        chrom, position, allele_code = item
        panel_records[chrom] += 1
        low, high = divmod(allele_code, 4)
        complement_low, complement_high = sorted((3 - low, 3 - high))
        complement_code = complement_low * 4 + complement_high
        for delta in deltas:
            if (position + delta) * 16 + allele_code in ancient:
                matches[delta][chrom] += 1
            elif (position + delta) * 16 + complement_code in ancient:
                strand_only[delta][chrom] += 1

for chrom in map(str, range(1, 11)):
    values = " ".join(
        f"delta_{delta:+d}={matches[delta][chrom]}"
        f"(+strand={strand_only[delta][chrom]})"
        for delta in deltas
    )
    print("CHROM", chrom, "panel_records", panel_records[chrom], values)
for delta in deltas:
    print(
        "TOTAL",
        f"delta_{delta:+d}",
        sum(matches[delta].values()),
        "strand_only",
        sum(strand_only[delta].values()),
    )
