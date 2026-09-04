# Ancient-sample simulation and ARG inference

We simulated a 10 Mb chromosome with `msprime` under a single constant-size
population with diploid effective size `Ne = 100,000`. The per-base,
per-generation recombination and mutation rates were both `1e-8`. Each
simulation contained 26 haploid samples at generation 0 and one haploid ancient
sample at generation 2,500. The mutated tree sequence was saved as
`known_arg.trees`, providing the known (simulated) ARG, and a VCF containing all
27 samples was also written. We used fixed ancestry/mutation seed pairs
`731029/731030` for `simarg` and `913021/913022` for `simarg2`.

SINGER inference used only the 26 modern haploid samples. We simplified the
simulated tree sequence to its generation-0 samples and wrote a separate
gzipped VCF because a standard VCF does not encode ancient sampling time and
the workflow expects an even number of haplotypes. Variants were treated as
polarised, and the mutation and recombination rates supplied to the workflow
were both `1e-8`.

We ran `singer-snakemake` with its default inference settings. The chromosome
was divided into approximately 1 Mb chunks; SINGER generated 100 MCMC ARG
samples per chunk with a thinning interval of 100 and 50% burn-in for summary
statistics. POLEGON then dated the inferred topologies using 100 samples,
thinning 10, and 50% burn-in. Corresponding samples from all chunks were merged
to produce 100 chromosome-wide inferred tree sequences (`.tsz` files), along
with diagnostics and coalescence-rate summaries. Runs were executed through
SLURM with 32 CPUs and 256 GB RAM. SINGER seeds were `1` for `simarg` and `2`
for `simarg2`.

Thus, `known_arg.trees` contains the 26 modern samples plus the ancient sample
and serves as simulation truth, whereas the SINGER ARG estimates are based on
the 26 modern samples alone.

