"""Draw the true three-epoch demography with demesdraw."""
import csv
from pathlib import Path
import demes
import demesdraw
import matplotlib.pyplot as plt

root = Path("msprime_three_epoch/demography")
with (root / "three_epoch.tsv").open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
builder = demes.Builder(time_units="generations", description="True three-epoch simulation demography")
builder.add_deme("population", epochs=[
    {"end_time": float(row["time_left"]),
     "start_size": float(row["effective_population_size"]), "size_function": "constant"}
    for row in reversed(rows)
])
graph = builder.resolve()
demes.dump(graph, root / "three_epoch.demes.yaml")
fig, axes = plt.subplots(1, 2, figsize=(10, 6))
demesdraw.tubes(graph, ax=axes[0], max_time=20000, labels=None,
                colours={"population": "#367a9c"}, title="Population-size schematic")
axes[0].set_ylabel("Generations before present")
for time, size in [(500, 20000), (5500, 5000), (15000, 50000)]:
    axes[0].text(0, time, f"Ne = {size:,}", ha="center", va="center",
                 fontsize=10, bbox={"facecolor": "white", "alpha": .8, "edgecolor": "none"})
axes[0].set_xlabel("Tube width proportional to diploid Ne")
demesdraw.size_history(graph, ax=axes[1], colours={"population": "#367a9c"},
                       title="True Ne trajectory")
axes[1].set(xlim=(0,20000), ylim=(0,55000), xlabel="Generations before present", ylabel="Diploid effective population size (Ne)")
axes[1].yaxis.label.set_rotation(90)
axes[1].yaxis.set_label_coords(-.15, .5)
for patch in list(axes[1].patches):
    patch.remove()
axes[1].hlines(50000, 10000, 20000, color="#367a9c", lw=2)
axes[1].axvline(1000, color="0.6", ls=":", lw=1)
axes[1].axvline(10000, color="0.6", ls=":", lw=1)
axes[1].set_xticks([0, 1000, 5000, 10000, 15000, 20000])
axes[1].tick_params(axis="x", rotation=45)
fig.suptitle("Three-epoch demography (older than 20,000 generations cropped)")
fig.tight_layout()
fig.savefig(root / "three_epoch_demesdraw.png", dpi=200)
fig.savefig(root / "three_epoch_demesdraw.pdf")
