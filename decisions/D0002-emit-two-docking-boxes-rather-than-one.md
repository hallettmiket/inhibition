---
id: D0002
title: Emit two docking boxes rather than one
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallet'
origin: adversary
supersedes: []
superseded_by: null
affects:
  - config/receptor.yaml
  - shared/receptor_prep.py
evidence:
  - 'QT7 is covalent at Cys113, so a tight box is centred on the warhead sub-pocket'
  - 'covalent box 20 A (t3,t4); expanded box 26 A (t1,t2)'
  - 'Cys113 SG is 4.26 A from the box centre, inside both'
runbook: docs/runbooks/receptor_selection.md
---

## Context
Adversary finding M5: a box drawn around a covalent ligand is centred on the warhead sub-pocket. That is right for the covalent approaches, which attack that atom, and wrong for the non-covalent ones, which would be biased toward a sub-pocket they have no reason to prefer.

## Decision
receptor_prep.py emits box.json (20 A, covalent, T_3/T_4) and box_expanded.json (26 A, full PPIase pocket, T_1/T_2). Each records which approaches use it.

## Consequences
Non-covalent and covalent dock scores are computed in different volumes and are therefore NOT directly comparable even before the tool difference is considered. This reinforces the no-authoritative-cross-approach-join design.
