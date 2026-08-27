---
id: D0098
title: Rank on warhead engagement geometry — the incumbent column predicts the MD outcome no better than chance, and mode aggregates fail because modes are mixtures
date: 2026-08-27
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/engagement_rank.py
  - shared/nac_criterion.py
  - scripts/rank_v2.py
  - exp/22_warhead_engagement_ranking/run_all.py
  - exp/22_warhead_engagement_ranking/splitter_effect.py
  - decisions/D0097-the-contact-splitter-is-built-and-selectable-and-five-pose-frameworks-are-mapped.md
evidence:
  - '@tt8804: "now we need to fix ranking accordingly, focusing on warhead optimal engagement only for now" and "we can rank by poses(modes) geometry only and have an option to rank by ligand by average geometry score across modes"'
  - 'TARGET: frac_attack_ready, the fraction of swept MD frames in which the warhead sits in the near-attack window, over 147 modes with an ok sweep result'
  - 'THE INCUMBENT PREDICTS NOTHING: conditional_eb rho = -0.015, p = 0.855, top-20 hit rate 45% against a 44% base rate (lift 1.02x)'
  - 'enrichment and viable_fraction are slightly NEGATIVE: rho = -0.043, p = 0.605, lift 0.90x'
  - 'THE SIMULATED POSE PREDICTS IT STRONGLY: anchor quality of the pose actually swept, rho = +0.652, p = 3.5e-19; the binary start_attack_ready gives rho = +0.462, top-20 hit 95%, lift 2.15x, and correlates with mode size at only +0.013'
  - 'start_dist_a alone gives rho = -0.664 -- closer is better, and it is the single strongest column measured'
  - 'EVERY MODE-LEVEL AGGREGATE IS WEAK: anchor_quality_max +0.130, p90 +0.124, q75_mean +0.109, mean -0.031, median -0.063'
  - 'anchor_quality_max correlates with mode size at +0.638 -- it is substantially a size proxy, and argmax anchoring is separately measured as the WORST representative selector (6.7% crystal recovery against 33.3% for the medoid of the well-anchored quartile)'
  - 'WHY THE AGGREGATES FAIL: a shipped mode is a mixture. Median within-mode spread of per-pose engagement is 0.776 on a scale that runs 0 to 1, and 93% of modes span more than half of it'
  - 'THE EXPLANATION MAKES A PREDICTION AND IT HOLDS: same 10 molecules, same 300 runs, same seed, only the splitter changed -- spread falls 0.658 -> 0.288 (2.3x) and modes spanning more than half the scale fall from 71% to 27%'
  - 'SILENT DEFECT FOUND AND FIXED: nac_criterion.anchor_quality returned 0.0 for an unmapped mechanism name instead of raising, so a typo would rank every pose of the affected molecules last with no error anywhere'
  - 'BLOCKED: 0 of 561 persisted representative files for nac_v5 carry pose_idx -- they predate the #76 fix -- so the representative geometry cannot be recovered for all 4,432 modes from disk'
runbook: python exp/22_warhead_engagement_ranking/run_all.py; python exp/22_warhead_engagement_ranking/splitter_effect.py
---

# D0098 — rank on engagement geometry

## The measurement the project had never made

`rank_validated` is False on every shortlist because the enrichment gate fired
(D0041) — but that gate asks whether docking separates actives from inactives.
This asks something answerable from data already on disk: **does a mode's static
score predict how much of an MD trajectory it spends in attack geometry?**

147 swept modes carry `frac_attack_ready`. Base rate above the 1% floor: 44%.

| metric | ρ | p | top-20 hit | lift | ρ vs size |
|---|---:|---:|---:|---:|---:|
| **the simulated pose's anchor quality** | **+0.652** | 3.5e-19 | — | — | — |
| `start_attack_ready` (binary) | +0.462 | 3.7e-09 | **95%** | **2.15×** | +0.013 |
| `anchor_quality_max` | +0.130 | 0.12 | 60% | 1.36× | **+0.638** |
| `anchor_quality_p90` | +0.124 | 0.13 | 60% | 1.36× | +0.366 |
| mode size (control) | +0.102 | 0.22 | 35% | 0.79× | +1.000 |
| **`conditional_eb` — THE INCUMBENT** | **−0.015** | **0.86** | 45% | **1.02×** | −0.094 |
| `enrichment` / `viable_fraction` | −0.043 | 0.61 | 40% | 0.90× | −0.155 |

**The column the pipeline ranks on today is indistinguishable from chance**, and
the frequency statistics are slightly negative. Geometry of one pose is strongly
predictive. `anchor_quality_max` looks like the exception until its size
correlation is read: +0.638, and argmax anchoring is separately the worst
representative selector measured (6.7% against 33.3%).

## Why every mode-level aggregate fails

Not because the statistics are wrong. Because **a shipped mode is a mixture**:
the median mode's poses span **0.776** of the anchor-quality scale, which itself
only runs 0 to 1, and **93%** span more than half of it. A group containing both
an ideal attack geometry and a hopeless one has no summary worth ranking, and
averaging it is how a real signal (+0.652 for one pose) becomes noise (−0.031
for the mean of its group).

This is D0088 restated in the ranking's own currency, and it is why *fixing
ranking* and *fixing splitting* were never separable.

## The explanation predicts something, and it holds

Same ten molecules, same 300 runs, same seed, only `--split-method` changed:

| | groups | poses/group | spread (median) | spanning >½ the scale |
|---|---:|---:|---:|---:|
| `warhead_dbscan` | 42 | 64.6 | 0.658 | 71% |
| `contact_linkage` | 962 | 3.1 | **0.288** | **27%** |

A **2.3× reduction** in exactly the quantity the explanation blames. It does
*not* show the aggregate now predicts the outcome — that needs modes swept under
the new splitter, and none have been.

## What was built

`shared/engagement_rank.py`:

* `mode_engagement(poses, statistic)` — one score per group, plus its **spread**,
  because a group whose members disagree must be visible as such;
* `rank_modes(...)` — geometry only, ranked within warhead class (an SN2 backside
  attack and a perpendicular approach do not span the same angles);
* `rank_ligands(modes, how="mean"|"best"|"median")` — the per-ligand option, with
  `n_modes` always reported beside the mean, because that denominator moves with
  docking depth (D0092).

**The pose-count gate is not applied by default.** It exists because a *frequency*
over three poses means nothing; an engagement score is a property of one pose and
is as estimable in a group of one as in a group of fifty. A caller that wants the
old gate asks for it.

Unknown statistics and aggregations raise rather than defaulting — a typo must
not silently produce a differently-ordered table.

## A silent defect found on the way

`nac_criterion.anchor_quality` returned **0.0** for an unmapped mechanism name
instead of raising. A typo would have scored every pose of the affected molecules
zero, ranking them last, with no exception anywhere — the symptom would have been
a warhead class that simply never appeared near the top. It now raises, following
`canonical_class()`'s rule that returning an unmapped value unchanged is how it
becomes its own category. Verified first that every mechanism in the per-pose
table is mapped.

## What is NOT settled

* **Range restriction, and it is not fixable here.** The 147 modes were *selected*
  for sweeping by `conditional_eb`. Every correlation above is measured inside the
  band the incumbent already liked. Metric-against-metric comparison is fair — all
  face the same restriction — but no absolute number is a population estimate.
* **`frac_attack_ready` is reachability of attack geometry**, not binding and not
  reactivity. Nothing here makes a shortlist `rank_validated`.
* **The representative's geometry cannot be recovered for all 4,432 modes.** All
  561 persisted representative files for `nac_v5` predate the #76 fix and carry no
  `pose_idx`, so they cannot be joined to their own measurements — the exact
  defect #76 was raised to close, with the fix in code and the data never
  regenerated. Until a re-screen, the strongest metric is computable only for
  modes that were swept.
* **No MD has been run under `contact_linkage`**, so the aggregate's predictive
  power under tight groups is a prediction, not a result.
