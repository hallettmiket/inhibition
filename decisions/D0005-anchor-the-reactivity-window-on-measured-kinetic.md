---
id: D0005
title: Anchor the reactivity window on measured kinetics, and never rank by LUMO
date: 2026-07-27
status: accepted
approach: t4
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - data/reference/pin1_reactivity_kinetics_1.csv
  - config/gates.yaml
evidence:
  - 'Pearson r = 0.396 between intrinsic k and Pin1 labeling across 8 compounds'
  - 'k spans 13.6x: 0.005 to 0.068 M-1 s-1'
  - '4e: near-lowest k (0.007) but highest labeling (97%)'
  - '4a: among highest k (0.030) but lowest labeling (17%)'
  - 'values digitized from Reddi 2023 Figure 5C, ~1 significant figure'
runbook: docs/runbooks/resolving_unverified_structures.md
---

## Context
The spec bounds T_4's reactivity window using computed LUMO. Figure 5C of Reddi 2023 supplies MEASURED second-order rate constants for Sulfopin and 4a-4g, which prompted checking whether reactivity predicts engagement.

## Decision
The window is anchored on measured kinetics where available, with computed LUMO used to place NEW warheads on that calibrated scale rather than as the source of truth. T_4 must not rank candidates by LUMO: within the precedented range, electrophilicity does not predict engagement - recognition does.

## Consequences
The reactivity window is a SAFETY filter for condition (ii), not a potency signal. Ranking by LUMO would have been an easy and invisible mistake. The kinetics values are figure-digitized and flagged as such; exact values need the SI tables.
