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

**Why it does not bite here.** Our sites are ascertained as polymorphic in a
discovery panel of ~1500 haplotypes that **contains** the 26 ARG-panel haplotypes.
Inference conditions on $1 \le d_0 \le 25$ — both alleles observed among the 26 —
and because those 26 are a *subset* of the 1500, that implies the site is polymorphic
in the discovery panel. So $A$ has probability 1 given conditioning we already
perform, and

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
