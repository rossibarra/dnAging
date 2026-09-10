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


def encoded(fields):
    chrom = fields[0]
    if chrom not in offsets:
        return None
    ref = base.get(fields[3])
    alt = base.get(fields[4])
    global_position = offsets[chrom] + int(fields[1])
    allele = None if ref is None or alt is None else global_position * 16 + ref * 4 + alt
    swapped = None if ref is None or alt is None else global_position * 16 + alt * 4 + ref
    return chrom, global_position, allele, swapped


ancient_positions = set()
ancient_alleles = set()
ancient_records = collections.Counter()
with open(ancient_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.split("\t", 5)
        item = encoded(fields)
        if item is None:
            continue
        chrom, position, allele, _ = item
        ancient_records[chrom] += 1
        ancient_positions.add(position)
        if allele is not None:
            ancient_alleles.add(allele)

panel_records = collections.Counter()
same_position = collections.Counter()
exact = collections.Counter()
swapped = collections.Counter()
mismatch = collections.Counter()
with open(panel_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.split("\t", 5)
        item = encoded(fields)
        if item is None:
            continue
        chrom, position, allele, reverse = item
        panel_records[chrom] += 1
        if position not in ancient_positions:
            continue
        same_position[chrom] += 1
        if allele is not None and allele in ancient_alleles:
            exact[chrom] += 1
        elif reverse is not None and reverse in ancient_alleles:
            swapped[chrom] += 1
        else:
            mismatch[chrom] += 1

for chrom in map(str, range(1, 11)):
    print(
        "CHROM",
        chrom,
        "panel_records",
        panel_records[chrom],
        "ancient_records",
        ancient_records[chrom],
        "same_position",
        same_position[chrom],
        "exact",
        exact[chrom],
        "swapped",
        swapped[chrom],
        "mismatch",
        mismatch[chrom],
    )
print("ANCIENT_UNIQUE_POSITIONS", len(ancient_positions))
print("ANCIENT_UNIQUE_ALLELES", len(ancient_alleles))
