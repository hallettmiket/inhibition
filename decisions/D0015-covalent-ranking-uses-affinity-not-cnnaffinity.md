---
id: D0015
title: Covalent ranking uses gnina affinity, not CNNaffinity
date: 2026-07-27
status: partially_withdrawn
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null   # justification revised by D0028
affects:
  - config/choreography.yaml
  - shared/covalent_protocol.py
  - docs/approaches/t3.md
  - docs/approaches/t4.md
evidence:
  - 'covalent affinity_kcal: ROC-AUC 0.815, CI [0.667, 0.931] EXCLUDES 0.5, EF1% 16.7'
  - 'covalent CNNaffinity: ROC-AUC 0.707, CI [0.408, 0.921] INCLUDES 0.5, EF1% 0.0'
  - 'EF1% of 0.0 means not one known active landed in the top 1% of the CNNaffinity ranking'
  - 'gnina itself warns CNN scoring is not calibrated for covalent docking'
  - 'both verdicts capped at UNDERPOWERED by the 4-chemotype floor'
runbook: null
---

## REVISION NOTICE (2026-07-28, D0028)

**The conclusion stands; the evidence below does not.** Re-measuring on
adduct-form ligands (D0022) gives `affinity_kcal` ROC-AUC **0.718** with a CI of
**[0.483, 0.944] that INCLUDES 0.5** — where this record cites 0.815 with an
interval that excluded it. The excluded-interval comparison was this decision's
stated reason for preferring affinity over CNNaffinity, and it is gone.

`affinity_kcal` remains the rank metric because it still beats CNNaffinity on
every statistic (AUC 0.718 vs 0.392, EF1% 19.0 vs 0.0, BEDROC 0.333 vs 0.146)
and because gnina warns CNN scoring is uncalibrated for covalent docking. Do not
cite the interval. See D0028, which also documents a class-imbalance confound in
the decoy set affecting both measurements.

## Context

D0011 recorded that the Rev 3 spec makes gnina `CNNaffinity` T_3's rank metric
and T_4's secondary, while gnina itself prints a warning on every covalent run
that CNN scoring is not calibrated for covalent docking. Rather than decide on
the warning alone, D0011 accepted the empirical route: run both metrics through
the enrichment gate and let the measurement decide.

M3 ran 6 actives against 294 warhead-bearing decoys.

## Decision

**T_3 and T_4 rank on gnina's Vina-style `affinity` (kcal/mol, lower better).**
`CNNaffinity` is carried as an advisory annotation, never as a rank metric.

The two metrics separate cleanly on the evidence that matters. `affinity_kcal`
enriches: its confidence interval excludes 0.5 and it puts actives in the top
1% (EF1% 16.7). `CNNaffinity` does neither — its interval includes 0.5, which is
consistent with no enrichment at all, and an EF1% of 0.0 means **not one known
active reached the top 1%** of its ranking.

That is independent confirmation of the tool's own warning, arrived at without
relying on it.

## Consequences

T_3's output contract changes: the rank metric is now kcal/mol and
**lower-is-better**, where the spec had it dimensionless and higher-is-better.
The direction is stated explicitly in several places and all of them move
together. Section 7's within-covalent re-score follows the same metric.

The protocol keeps `cnn_scoring: rescore` so the advisory number is still
produced; only its ROLE changes. The fingerprint is therefore unaffected by this
decision.

Both verdicts remain UNDERPOWERED at 4 independent chemotypes, so this is a
comparison between two metrics on the same data rather than a claim that
covalent docking is validated. Under D0012 the ranking carries forward with its
uncertainty displayed.
