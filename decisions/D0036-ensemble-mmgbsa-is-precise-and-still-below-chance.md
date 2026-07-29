---
id: D0036
title: Better sampling does not rescue MM-GBSA — the ensemble is precise and still below chance
date: 2026-07-29
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/mmgbsa_ensemble.py
  - scripts/run_mmgbsa_ensemble.py
  - decisions/D0031-class-matched-decoys-remove-the-apparent-covalent-enrichment.md
  - decisions/D0032-mmgbsa-gate-and-the-power-floor-on-negative-verdicts.md
evidence:
  - '167 candidates rescored per-frame over 2 ns GB implicit-solvent MD, 0 failures'
  - 'gate set complete: 82 of 83 ligands, 2 actives, 2 chemotypes'
  - 'docking ROC-AUC 0.537; MM-GBSA single-structure 0.350; MM-GBSA ensemble 0.394'
  - 'propagating each candidate MEASURED SEM into the metric: AUC 95% [0.356, 0.463], P(AUC>0.5) = 0.002'
  - 'Sulfopin dG -7.58 +/- 0.28; 50 of 80 decoys score better'
  - 'Juglone dG -7.76 +/- 0.51; 47 of 80 decoys score better'
  - 'decoy dG median -8.91, sd 6.35; beating 80% of decoys needs dG < -11.38'
  - 'one-frame vs 90-frame Spearman on the gate set is 0.283, so a SINGLE structure there is noise-dominated'
  - 'the 90-frame mean is not: SEM 0.28-0.51 kcal/mol against a decoy spread of 6.35'
  - 'verdict remains UNDERPOWERED: 2 actives < 3 floor, 2 chemotypes < 6 floor'
---

# The ensemble is precise, and still below chance

## The prediction was wrong

Before the run I expected the ensemble to show that actives and class-matched
decoys were **never resolvable** at this noise level — a wide interval
straddling 0.5, converting D0031/D0032's negatives into "the measurement could
never have decided". That is not what came back.

Propagating each candidate's own measured SEM into the metric gives a 95%
interval of **[0.356, 0.463]** and **P(AUC > 0.5) = 0.002**. The ensemble is
not too imprecise to decide. It decides, and it ranks known actives *below*
property-matched decoys.

## The numbers

| metric | ROC-AUC | active ranks (of 82) |
|---|---|---|
| docking `affinity_kcal` | 0.537 | 29, 49 |
| MM-GBSA single-structure | 0.350 | 44, 63 |
| **MM-GBSA ensemble** | **0.394** | **48, 52** |

| active | dG (kcal/mol) | decoys scoring better |
|---|---|---|
| Sulfopin | -7.58 +/- 0.28 | 50 / 80 |
| Juglone | -7.76 +/- 0.51 | 47 / 80 |

Both actives sit mid-pack. The decoy distribution has median -8.91 and sd 6.35,
so beating 80% of decoys requires dG below **-11.38**; the actives are at -7.6
and -7.8. They are not marginally misplaced — they are where you would expect a
random member of the set to land.

## Why this is a stronger result than D0032, not a repeat of it

D0032's MM-GBSA number rested on a single minimised structure, and the
one-frame-versus-ninety control shows why that was fragile: on **this** set a
single frame recovers the 90-frame ranking at Spearman only **0.283**. A single
structure on class-matched decoys is noise-dominated, exactly where the true
differences are smallest.

Averaging 90 frames removes that. The SEM falls to 0.28-0.51 kcal/mol against a
decoy spread of 6.35 — better than an order of magnitude of headroom. So the
earlier negative was **not** an artefact of poor sampling, which was the most
obvious remaining explanation for it. Fixing the sampling moves the point
estimate from 0.350 to 0.394 and leaves it firmly below chance.

That closes off the cheapest available defence of the method on this target.

## What the measurement supports, and what it does not

**Supports:** the below-chance point estimate is not measurement noise. Within
this gate set, ensemble MM-GBSA orders Sulfopin and Juglone below the median
class-matched decoy, with an uncertainty far too small to explain the gap.

**Does not support:** "MM-GBSA fails on Pin1." That is an inference from **two**
molecules and **two** chemotypes, and the gate's own floors (3 actives, 6
chemotypes) exist precisely to refuse it. The verdict stays UNDERPOWERED.

These two statements are not in tension, and keeping them apart is the whole
point. *Measurement precision* and *statistical power for a general claim* are
different quantities. This run has the first in abundance and the second not at
all. A confident number computed from two molecules is still a number computed
from two molecules.

## What would change the answer

Not more sampling. That was this experiment, and it is now spent — going from
one structure to 90 frames moved the point estimate by +0.044 while shrinking
the error bar by roughly a factor of five.

What remains:

1. **More actives.** Still the binding constraint, still a literature problem
   rather than a compute one. Two is below the floor and no amount of GPU time
   changes that.
2. **Explicit solvent.** This tier is GB implicit — no water structure, no
   viscosity, a too-fast conformational clock. It is the tier the T5 spec did
   NOT ask for. Whether explicit solvent behaves differently is untested here,
   and it is the one remaining fidelity axis that has not been ruled out.
3. **The pocket itself.** Two independent estimators of very different cost now
   both fail to separate one known active from same-chemotype decoys on this
   receptor. That continues to point at the target — shallow, solvent-exposed,
   the regime where structure-based scoring is weakest — rather than at either
   scoring function.

## A note on what the error bar bought

The uncertainty was not decoration. Without it the honest reading of 0.394
would have been "below chance, but who knows" — indistinguishable from
D0032's 0.350. With per-candidate SEM propagated into the metric, the reading
becomes "below chance, and measurement error is not why". The error bar is what
converted an ambiguous negative into a specific one.
