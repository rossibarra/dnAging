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
