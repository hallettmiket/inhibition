---
id: D0049
title: Ranking is size-decorrelated by local median within heavy-atom strata, not by ligand efficiency
date: 2026-08-01
status: accepted
approach: shared
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/rank_shortlist.py
  - scripts/measure_size_correlation.py
  - tests/test_rank_within_class.py
evidence:
  - 'D0043 baseline reproduced exactly on the current frames: T_1 raw rho -0.617, LE -0.938'
  - 'ligand efficiency is WORSE than the raw score in 5 of 6 pools: T_1 -0.938 vs -0.617; T_2/atra -0.591 vs -0.262; T_2/du_xu -0.654 vs +0.119; T_2/guo -0.681 vs -0.088; T_4 -0.480 vs -0.025'
  - 'the one exception is T_3: LE -0.208 vs raw -0.695'
  - 'D0047 moved T_3 from -0.479 (pre-fix column) to -0.695; and T_4 from +0.181 to -0.025'
  - 'T_2 size bias is seed-dependent and changes sign: atra -0.262, guo -0.088, du_xu +0.119'
  - 'straight-line residual in kcal/mol leaves T_4 at -0.247'
  - 'residual refitted in RANK space leaves ~0.26 under heteroscedastic scatter'
  - 'adopted binned local-median residual: T_1 +0.028, T_2/atra +0.032, T_2/du_xu -0.002, T_2/guo +0.031, T_3 0.000, T_4 +0.007'
  - 'shortlist churn at top-10: T_1 8/10 kept, T_2/atra 4/10, T_2/du_xu 2/10, T_2/guo 7/10, T_3 7/10, T_4 25/70'
---

# Rank on how good a molecule is *for its size*

## Context

D0043 established that our rankings are partly a molecular-size sort, measured
Spearman rho = -0.617 (T_1) and -0.479 (T_3) against heavy-atom count, and
recorded the fix as **open by design** rather than quietly choosing one. It also
measured ligand efficiency at -0.938 and rejected it as the mirror-image bias.

Issue #9 item 4, from the meeting with Ian, asks us to "normalise to ligand
efficiency for docking to avoid favouring larger molecules". That is the right
instinct about the problem and the wrong instrument, and three inputs had
genuinely changed since D0043 was written, so it was re-measured rather than
answered by citation:

1. **D0047** — `affinity_kcal` was the CNN-selected pose's affinity. T_3 and
   T_4's correlations were computed on the pre-fix column.
2. **Ligands are protonated at pH 7.4**, not docked as drawn.
3. **T_2 is five seed neighbourhoods**, not ATRA alone.

## Decision

**Rank on the size-decorrelated residual**: within equal-population strata of
heavy-atom count, subtract the stratum's **median** score. Ranking then asks
*how good is this molecule compared with others of its size*.

**Ligand efficiency is not adopted as a sort key.** It remains computed and
displayed, because chemists read it and D0043 rejects it as a *ranking* metric,
not as a reported quantity.

## Why not the simpler options

| candidate | result |
|---|---|
| ligand efficiency | worse than the raw score in 5 of 6 pools |
| straight-line residual, kcal/mol | leaves T_4 at rho -0.247 — the dependence is monotone but not linear |
| residual refitted in rank space | leaves ~0.26 — scatter around the trend widens with size, and a line through rank space cannot absorb that |
| **binned local median** | **every arm to \|rho\| <= 0.032** |

The first two were implemented and rejected *by the test*, not by argument:
`test_decorrelation_actually_removes_the_size_dependence` failed on both.

Local centring assumes nothing about the shape of the size-score relationship,
and the median makes it robust to the outliers a docking score reliably
produces. The residual stays in **kcal/mol** — "this much better than the
typical molecule of its size" — rather than becoming a plausible float in
unstated units, which is the failure shape `how_this_project_breaks.md`
catalogues.

## Consequences

**Shortlists churn substantially**, most where the raw size bias was furthest
from zero — T_2/du_xu keeps only 2 of its top 10, because its raw correlation
was *positive* (+0.119) so decorrelation moves it most. This is a larger change
than D0047's.

`rank_raw_metric`, `rank_metric_used` and `rank_size_decorrelated` are carried
on every ranked frame, so a molecule that moved can be attributed to this
decision rather than leaving a reader to guess. `decorrelate_size=False`
restores the pre-decision ordering exactly, and a test pins that.

Below `MIN_DECORRELATION_N` (30) docked rows there is no trustworthy fit; those
frames fall back to the raw metric with `rank_size_decorrelated=False` stamped
rather than ranking on a fit from a handful of points.

**What this does NOT do.** It does not make the ranking valid. The enrichment
gate has still fired (D0041) and pose recovery is still 5% (D0046); every row
still carries `rank_validated = False`. Removing a known bias from an
unvalidated score leaves an unvalidated score — this makes the ordering *less
wrong in one identified way*, which is not the same as right.

## What would change this

A scorer that passes a powered gate (#4 Phase 2.1) would make the whole
question different: decorrelation is worth doing on a metric that discriminates
and is close to pointless on one that does not. Revisit then.
