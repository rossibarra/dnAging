# TODO

Deferred work is now tracked as
[GitHub issues](https://github.com/rossibarra/dnAging/issues). Items still
graduate from the working documents — [working_diffusion.md](working_diffusion.md),
[working_betabinom.md](working_betabinom.md), and
[working_insertion.md](working_insertion.md) — once a decision has been made
about *what* to do but the work itself is out of scope for the current pass; open
an issue for each rather than adding it here.

## Open items moved from this file

| Issue | Topic |
|---|---|
| [#1](https://github.com/rossibarra/dnAging/issues/1) | Insertion likelihood is pseudo-haploid only; README does not say so and the code does not enforce it |
| [#2](https://github.com/rossibarra/dnAging/issues/2) | Approach 1: heterozygotes silently become derived under `--ploidy 1` (het policy, counter, warning) |
| [#3](https://github.com/rossibarra/dnAging/issues/3) | Approach 1: check moment inequalities and genotype probabilities under `--ploidy 2` instead of clipping |
| [#4](https://github.com/rossibarra/dnAging/issues/4) | End-to-end msprime validation of $E[p_T \mid d_0, t_i]$ and the second moment |
| [#5](https://github.com/rossibarra/dnAging/issues/5) | Insertion: localise age bias with inferred SINGER+POLEGON ARGs and inferred Ne |

The former "Repository provenance gaps" section was resolved before the move:
every file it listed as untracked is now committed, and `.gitignore` carries the
`msprime_exact_time_ne50k/` entry. The full text of the old list is in git
history (commit f03b4954c).
