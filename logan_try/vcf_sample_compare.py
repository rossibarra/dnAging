import collections
import random
import sys


panel_path, ancient_path = sys.argv[1:]
chromosomes = {str(i) for i in range(1, 11)}
target_limit = 2000
rng = random.Random(20260908)
seen = collections.Counter()
reservoirs = {chrom: [] for chrom in chromosomes}

with open(panel_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.split("\t", 5)
        chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
        if chrom not in chromosomes:
            continue
        seen[chrom] += 1
        item = (pos, ref, alt)
        if len(reservoirs[chrom]) < target_limit:
            reservoirs[chrom].append(item)
        else:
            replacement = rng.randrange(seen[chrom])
            if replacement < target_limit:
                reservoirs[chrom][replacement] = item

targets = {chrom: {} for chrom in chromosomes}
for chrom, reservoir in reservoirs.items():
    for pos, ref, alt in reservoir:
        targets[chrom].setdefault(pos, set()).add((ref, alt))

observed = {chrom: collections.defaultdict(set) for chrom in chromosomes}
with open(ancient_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.split("\t", 5)
        chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
        if chrom in chromosomes and pos in targets[chrom]:
            observed[chrom][pos].add((ref, alt))

for chrom in map(str, range(1, 11)):
    positions_found = 0
    exact = 0
    swapped = 0
    mismatch = 0
    for pos, panel_alleles in targets[chrom].items():
        ancient_alleles = observed[chrom].get(pos, set())
        if not ancient_alleles:
            continue
        positions_found += 1
        for ref, alt in panel_alleles:
            if (ref, alt) in ancient_alleles:
                exact += 1
            elif (alt, ref) in ancient_alleles:
                swapped += 1
            else:
                mismatch += 1
                if mismatch <= 5:
                    print(
                        "MISMATCH",
                        chrom,
                        pos,
                        "panel",
                        (ref, alt),
                        "ancient",
                        sorted(ancient_alleles),
                    )
    print(
        "CHROM",
        chrom,
        "panel_positions_sampled",
        len(targets[chrom]),
        "positions_found",
        positions_found,
        "exact",
        exact,
        "swapped",
        swapped,
        "mismatch",
        mismatch,
    )
