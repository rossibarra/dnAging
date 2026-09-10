import collections
import sys


panel_path, ancient_path, samples_path = sys.argv[1:]
selected = [x.strip() for x in open(samples_path) if x.strip()]


def records(path):
    with open(path) as handle:
        sample_names = None
        for line in handle:
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
                yield ("HEADER", sample_names)
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            yield (
                fields[0],
                int(fields[1]),
                fields[3],
                fields[4],
                fields[8],
                fields[9:],
            )


panel = records(panel_path)
ancient = records(ancient_path)
phead = next(panel)
ahead = next(ancient)
assert phead[0] == ahead[0] == "HEADER"
ancient_index = {sample: i for i, sample in enumerate(ahead[1])}
missing = [sample for sample in selected if sample not in ancient_index]
assert not missing, f"selected samples missing from ancient VCF: {missing}"
selected_index = [ancient_index[sample] for sample in selected]

file_chromosomes = ("1", "10", "2", "3", "4", "5", "6", "7", "8", "9")
chrom_order = {chrom: i for i, chrom in enumerate(file_chromosomes)}
overlap_positions = collections.Counter()
allele_match = collections.Counter()
allele_mismatch = collections.Counter()
called = collections.Counter()
heterozygous = collections.Counter()


def key(record):
    return (chrom_order.get(record[0], 10**9), record[1])


def count_selected_genotypes(record):
    format_fields = record[4].split(":")
    try:
        gt_index = format_fields.index("GT")
    except ValueError:
        return
    for sample, index in zip(selected, selected_index):
        fields = record[5][index].split(":")
        if gt_index >= len(fields):
            continue
        alleles = fields[gt_index].replace("|", "/").split("/")
        values = [allele for allele in alleles if allele != "."]
        if values:
            called[sample] += 1
        if len(set(values)) > 1:
            heterozygous[sample] += 1


def grouped(iterator):
    record = next(iterator, None)
    while record is not None:
        group_key = key(record)
        group = [record]
        record = next(iterator, None)
        while record is not None and key(record) == group_key:
            group.append(record)
            record = next(iterator, None)
        yield group_key, group


panel_groups = iter(grouped(panel))
ancient_groups = iter(grouped(ancient))
panel_item = next(panel_groups, None)
ancient_item = next(ancient_groups, None)
while panel_item is not None and ancient_item is not None:
    panel_key, panel_group = panel_item
    ancient_key, ancient_group = ancient_item
    if panel_key < ancient_key:
        panel_item = next(panel_groups, None)
        continue
    if ancient_key < panel_key:
        for record in ancient_group:
            count_selected_genotypes(record)
        ancient_item = next(ancient_groups, None)
        continue
    chrom = panel_group[0][0]
    overlap_positions[chrom] += 1
    ancient_alleles = {(record[2], record[3]) for record in ancient_group}
    for record in panel_group:
        if record[2:4] in ancient_alleles:
            allele_match[chrom] += 1
        else:
            allele_mismatch[chrom] += 1
            if allele_mismatch[chrom] <= 5:
                print(
                    "MISMATCH",
                    chrom,
                    record[1],
                    record[2],
                    record[3],
                    "ancient",
                    sorted(ancient_alleles),
                )
    for record in ancient_group:
        count_selected_genotypes(record)
    panel_item = next(panel_groups, None)
    ancient_item = next(ancient_groups, None)

while ancient_item is not None:
    _, ancient_group = ancient_item
    for record in ancient_group:
        count_selected_genotypes(record)
    ancient_item = next(ancient_groups, None)

print("VCF_SAMPLES", len(phead[1]), len(ahead[1]), len(selected))
for chrom in map(str, range(1, 11)):
    print(
        "CHROM",
        chrom,
        "overlap_positions",
        overlap_positions[chrom],
        "allele_match",
        allele_match[chrom],
        "allele_mismatch",
        allele_mismatch[chrom],
    )
print("CALLED_RANGE", min(called.values()), max(called.values()))
heterozygous_samples = {sample: n for sample, n in heterozygous.items() if n}
print("SAMPLES_WITH_HETEROZYGOUS_CALLS", len(heterozygous_samples))
for sample, count in sorted(heterozygous_samples.items()):
    print("HET", sample, count)
