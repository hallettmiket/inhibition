---
id: D0070
title: Neither metric's value converges, but consensus preserves rank order where frequency does not
date: 2026-08-06
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/pose_consensus.py
  - shared/composite_rank.py
  - docs/ranking_rationale.md
  - decisions/D0068-enrichment-depends-on-search-effort-and-energy-beats-it-at-convergence.md
evidence:
  - '15 crystallographic Cys113 positives, each docked at 200 and 2000 runs, BOTH metrics computed from the SAME .dlg at each effort'
  - 'frequency (viable-NAC fraction): median 0.185 -> 0.067, median |change| 0.118'
  - 'consensus (pairwise agreement, top_n=10 held fixed): median 0.533 -> 0.356, median |change| 0.289'
  - 'rank agreement across efforts: frequency Spearman rho = -0.047; consensus Spearman rho = +0.568'
  - 'both shift systematically downward with effort: Wilcoxon p = 0.0002 (frequency), p = 0.0059 (consensus)'
  - 'top_n held FIXED at 10 across both efforts, per pose_consensus.require_same_n'
---

# Consensus preserves rank order; frequency does not preserve anything

## What was asked

`shared/pose_consensus.py` was built on an argument, not a measurement: that
agreement among the **top-N** poses should be more stable than the viable-NAC
**fraction**, because the fraction is computed over every run and therefore
inherits how hard the search looked, while a top-N window does not. Its own
docstring said this was untested. `scripts/consensus_convergence.py` tests it.

Both metrics are computed from the **same docking run** at each effort, so the
comparison isolates the metric rather than confounding it with run-to-run
scatter. `top_n` is held fixed at 10 across efforts, because a consensus at
N = 10 and one at N = 50 are different quantities.

## The result, both halves

| | 200 runs | 2,000 runs | median \|change\| | rank agreement across efforts |
|---|---|---|---|---|
| frequency | 0.185 | 0.067 | **0.118** | **ρ = −0.047** |
| consensus | 0.533 | 0.356 | **0.289** | **ρ = +0.568** |

Both decline systematically with search effort (Wilcoxon p = 0.0002 and 0.0059).
**Consensus does not escape D0068** — its absolute value is search-dependent too,
and moves further in absolute terms.

But **frequency's rank correlation across efforts is −0.047**. The ordering
produced at 200 runs carries no information about the ordering at 2,000. Not a
weak relationship: none. Consensus at +0.568 is a moderate, usable one.

## The reading, and the part of it that is post-hoc

`consensus_convergence.py` fixed its reading in advance as: *consensus is the fix
for D0068 only if it moves LESS than the frequency on these same dockings.*

**That criterion was ambiguous and I am resolving the ambiguity after seeing the
data**, which is worth stating plainly rather than presenting the favourable
reading as though it had been specified. "Moves less" did not say *in what*:

- **In absolute value, consensus moves MORE** (0.289 vs 0.118). On the letter of
  the pre-registration, it fails.
- **In rank order, consensus moves far less** (ρ = +0.568 vs −0.047).

Rank is the relevant measure *for a ranking* — the composite consumes an ordering
and never an absolute level, and D0068's damage was to the ordering. So the rank
comparison is the one that bears on the decision. But it is a choice of measure
made with the numbers visible, and a reader is entitled to weigh it accordingly.
The absolute-change column is reported beside it for that reason.

## What follows

- **Neither metric may be quoted as an absolute value without its run count.**
  D0068's rule stands and now covers consensus as well.
- **Consensus is admitted as a ranking component**, on the rank evidence, and
  `shared/composite_rank.py` already consumes it as one component among several
  with its own uncertainty rather than as a score.
- **Frequency's ρ = −0.047 is the more important number here.** It is a stronger
  statement of D0068 than D0068 itself made: the 200-run ranking is not merely
  imprecise at the head, its *order* is uninformative about a better-sampled
  ordering of the same molecules. The screen's value is as a FILTER — the top 300
  are statistically indistinguishable from known crystallographic binders — and
  that survives, because a filter needs a threshold rather than an order.
- **The mechanism is the same for both.** At low effort an under-converged search
  returns poses that resemble each other because it has not looked anywhere else;
  at high effort it finds genuinely distinct minima. Agreement falls because the
  search got better, not because the molecule changed.

## What this does not settle

- **n = 15**, the same ceiling as everything else on this project.
- Only two efforts were compared. Whether consensus's rank stabilises further
  beyond 2,000 runs, or continues to drift, is unmeasured.
- `top_n = 10` was not itself varied. The rank stability may depend on it, and
  `require_same_n` exists precisely because that parameter is part of the
  quantity.
