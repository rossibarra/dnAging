# Modelling notes

Why the model is the way it is: approximations we make deliberately, the conditions
under which they hold, and what would break them.

This is the third of three documents, and the division is intentional:

| file | question it answers |
|---|---|
| [README.md](./README.md) | how do I run it? |
| [MATH.md](./MATH.md) | what is the model, and how is it derived? |
| **NOTES.md** | why is this approximation acceptable, and when would it stop being? |

Anything that is a judgement call rather than a derivation belongs here.

---

## Ascertainment is exactly ignorable under nesting

**The concern.** The likelihood conditions on the ascertained site set as if
ascertainment carried no information about $T$ (MATH.md §2). In general it does:
keeping a site because it was polymorphic in a discovery panel of $n_B$ haplotypes
reweights the present-frequency posterior by

$$
P(A \mid X) \;=\; 1 - X^{n_B} - (1-X)^{n_B},
$$

which is not constant in $X$. Since $X$ is linked to $p_T$ through the frequency
trajectory, conditioning on the ascertained set can inform $T$.

**Why it does not bite here.** Our sites come from a discovery panel of **1500
sequenced lines that contains the ARG panel's lines**, and **no frequency filter was
applied** — any SNP present in that panel was used. Both facts are load-bearing.
Inference conditions on $1 \le d_0 \le 25$ — both alleles observed among the 26 ARG
haplotypes — and because those are a *subset* of the discovery panel, that implies
the site is polymorphic in the discovery panel. So $A$ has probability 1 given
conditioning we already perform, and

$$
E[p_T \mid d_0, t_i, A] \;=\; E[p_T \mid d_0, t_i]
$$

holds **identically**, not asymptotically, with no dependence on $n_B$.

Note this is a condition relating the *discovery panel* to the *ARG panel*. Where the
ancient samples come from is a separate assumption (they must not be in either panel;
see MATH.md §7).

**What it would cost if the panels were disjoint.** Not much, but not zero. Under a
neutral $1/X$ prior with $n = 26$, treating the 1500 as an independent panel shifts
$E[X \mid d_0]$ by $+1.7\%$ at $d_0 = 1$ and by under $0.05\%$ for $d_0 \ge 2$. The
controlling quantity is the posterior mass over $X$ below $1/n_B$, which for $d_0 = 1$
is $0.40$ at $n_B = 52$ but only $0.005$ at $n_B = 5008$ — so the effect is confined
to singletons and vanishes quickly in panel size.

**What would actually break it.** A *frequency threshold* at discovery instead of
plain polymorphism. A MAF or minor-allele-count cutoff is **not** implied by
$d_0 \ge 1$, so nesting gives no protection. With $n_B = 1500$:

| discovery rule | shift in $E[X \mid d_0{=}1]$ | at $d_0 = 2$ |
|---|---|---|
| polymorphic, nested panel | $0\%$ (exact) | $0\%$ |
| minor-allele count $\ge 5$ | $+8.5\%$ | $+0.4\%$ |
| MAF $\ge 1\%$ | $+25.5\%$ | $+2.7\%$ |

So reusing this pipeline on a site set built with a frequency cutoff — or on a
discovery panel that does not contain the ARG haplotypes — means the ascertainment
factor must be modelled rather than ignored. Singletons are where it shows up first.

**What ascertainment still costs, even when exactly ignorable.** It removes young and
rare-in-discovery alleles from the site set, and those are exactly the sites carrying
the lower-bound age signal. That is a loss of *power*, not a bias — expect broad
posteriors, and treat a tight interval on a single sample with suspicion.

---

## Panel missingness is ignorable only if it is allele-blind

**The concern.** $d_0$ is counted over however many of the 26 ARG-panel haplotypes
are *called* at a site, and the likelihood uses the moment plane built for that
exact $n$ (MATH.md §5, eq. 7). Choosing the plane by $n$ makes the binomial
sampling term right for the number of haplotypes observed, but it is **not** a
correction for missingness. What eq. (7) needs is that the called subset $R$ is an
allele-blind sample of the panel,

$$
P(d_0 \mid n, x_0, R) \;=\; \binom{n}{d_0} x_0^{\,d_0}(1-x_0)^{\,n-d_0},
$$

i.e. missing at random with respect to the allele a haplotype carries, given $n$.

**Why conditioning on $n$ is not enough.** The $n$-specific plane conditions on
*how many* haplotypes were called, never on *which*. If one allele is systematically
harder to call — reference bias in mapping or genotyping, an ALT-specific filter,
low depth correlated with the derived haplotype background — then within the called
subset that allele is under-represented, $E[d_0 \mid n] \neq n\,x_0$, and the bias
passes straight through into $d_0$. The plane is then the correct likelihood for the
wrong count: a downward-biased $d_0$ presents as a rarer allele, which the
age-conditioned trajectory reads as a different frequency history.

**Direction is not established, and is not simply signed.** Reference bias acts on
the REF/ALT axis, but the trajectory is conditioned on the *derived* count, and
polarity flips between draws ($d_0 = c_{\text{alt}}$ or $n - c_{\text{alt}}$,
MATH.md eq. 10). The same ALT under-calling therefore pushes $d_0$ down at
ALT-derived sites and up at ALT-ancestral ones, so it does not translate into a
uniform shift in $\hat T$. Do not quote a direction without measuring it.

**The condition.** Sound when panel dropout is allele-blind — uniform low depth,
random per-haplotype missingness, a filter applied on position rather than on
genotype. Not sound when calling success depends on the allele itself. Nothing in
the pipeline can detect the difference: `sites_panel_below_min_n` counts only the
sites that fell *below* `--min-n`, and says nothing about whether the calls that
survived above it are allele-blind.

**The check, and the mitigation.** The check is upstream in the panel VCF: compare
the panel ALT fraction across called counts $n$, and look for a systematic drift as
$n$ falls. If there is one, the mitigation is `--min-n 26` — use only fully called
panel sites, where there is no missingness left to be non-ignorable — at a cost in
sites, and therefore in power, exactly where ascertainment already costs power.

---

## REF/ALT harmonisation is strand-blind

Panel and ancient VCFs are joined on position, and the panel ALT count is harmonised
to the ancient VCF's orientation (MATH.md §6): identical alleles are used as-is,
exactly swapped alleles give $n - c_{\text{alt}}$, and anything else is skipped and
counted in `sites_allele_mismatch`.

That comparison is on base *identity* only — `_BASE` maps `A,C,G,T` with no
complement operation — so for the **palindromic** classes **A/T** and **G/C** a
strand flip is indistinguishable from a genuine REF/ALT swap:

| ancient | panel | code concludes | correct |
|---|---|---|---|
| A/T | A/T, same strand | same → $c_{\text{alt}}$ | yes |
| A/T | T/A, genuine swap | swapped → $n - c_{\text{alt}}$ | yes |
| A/T | T/A, **strand flip** | swapped → $n - c_{\text{alt}}$ | **no** — should be $c_{\text{alt}}$ |
| A/T | A/T, **flip + swap** | same → $c_{\text{alt}}$ | **no** — should be $n - c_{\text{alt}}$ |

**Non-palindromic sites fail safe.** A strand flip on an A/G site yields T/C, which
matches neither orientation, so it falls through to the skip branch. Losing a site
costs power; inverting $d_0$ corrupts the likelihood, so this asymmetry is the
behaviour we want.

**The useful diagnostic.** Because non-palindromic sites are skipped and palindromic
ones are not, a strand disagreement announces itself as a *large*
`sites_allele_mismatch` count — and in exactly that situation the A/T and G/C sites
that passed the check are the ones that are silently wrong. A high count is
therefore not merely "these VCFs were normalised against different references"; it
is a reason to distrust the palindromic sites that *did* match.

**The condition.** Orientation harmonisation is sound iff both VCFs report on the
same strand. That holds for this data set, so no strand handling is implemented. If
strand provenance is ever uncertain, the standard mitigation is to drop A/T and G/C
sites outright (via `--include-positions`) rather than to infer strand from allele
frequencies.

---

## The lower table boundary is exact only for $T \ge t_{\min}$

The mutation-age table starts at a strictly positive age $t_{\min}$ (`--age-min`,
default 10 generations) because the age axis is log-spaced and cannot include zero.
The store's branch intervals have no such floor: `below` is $0$ for a branch reaching
the present.

**Why that is not a problem for $T \ge t_{\min}$.** Branch intervals extending below
$t_{\min}$ are integrated only over the covered part but normalised by the true
branch length, which is exact — every uncovered mutation age is younger than the
sample, so it contributes $p_T = 0$ (MATH.md, eq. 10b). Sites are therefore **not**
dropped for reaching below the table. Dropping them would be actively harmful: it
would discard the youngest mutations, which are precisely the sites carrying the
lower-bound age signal, manufacturing the same power loss ascertainment causes.

**The residual, and its size.** For $T < t_{\min}$ the identity fails, because
uncovered ages in $(T, t_{\min})$ do contribute. The posterior below $t_{\min}$ is
therefore biased low. With the default `--t-min 0 --t-max 30000 --n-t 300` the grid
step is ~100 generations, so exactly one grid point ($T=0$, a present-day sample)
falls below $t_{\min}=10$. Lower `--age-min`, or raise `--t-min`, if sample ages
that young are ever of interest.

**What the earlier behaviour cost.** Before the renormalisation, a straddling branch
was inflated by $(\text{above}-\text{below})/(\text{above}-t_{\min})$:

| branch | inflation of $\bar p$ |
|---|---|
| $[0, 50]$ | $1.25\times$ |
| $[0, 1000]$ | $1.010\times$ |
| $[0, 10000]$ | $1.001\times$ |

Since $\bar p$ is nonzero only when $T < \text{above}$, the large errors required a
sample younger than a few tens of generations; for ancient samples at $T$ in the
thousands only branches with $\text{above} > T$ contribute, so the old error stayed
under 1%. It was worth fixing exactly rather than bounding, because the fix is a
single factor.

**Diagnostic.** `sites_age_clipped_low` in `run.json` counts sites where some draw's
branch reached below $t_{\min}$. These sites are *used*, not skipped — the counter
records how often the renormalisation is doing work, so a large value is expected
for young-mutation-rich site sets and is not a warning.

**Do not mirror the correction at the upper boundary.** Uncovered ages above the
table are *older* than the sample and do contribute, so the same rescaling there
would be wrong. That end is handled by capping at `--mutation-age-max`. The code
comment in `phi_lookup` says so; the natural instinct is to symmetrise it.
