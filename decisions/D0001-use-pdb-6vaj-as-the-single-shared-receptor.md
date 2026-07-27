---
id: D0001
title: Use PDB 6VAJ as the single shared receptor
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallet'
origin: spec
supersedes: []
superseded_by: null
affects:
  - config/receptor.yaml
  - shared/receptor_prep.py
  - config/sources.yaml
evidence:
  - 'LINK record: SG CYS A 113 - C10 QT7 A 201 at 1.78 A (a real covalent bond length)'
  - 'resolution 1.42 A'
  - 'TITLE: CRYSTAL STRUCTURE ANALYSIS OF HUMAN PIN1'
  - 'QT7 ligand present, 16 atoms'
  - 'sha256 820fd5969131bef8... pinned in sources.lock.json'
runbook: docs/runbooks/receptor_selection.md
---

## Context
All four approaches must dock into the identical prepared receptor or their scores are not comparable. A covalently-bound reference ligand also defines the box for free and hands over the exact attachment atom.

## Decision
6VAJ is the shared receptor for the whole choreography. Verified first-hand rather than on the spec's assertion (this closed adversary finding M6). It is hash-pinned; changing it mid-choreography is forbidden.

## Consequences
Every cross-approach comparison depends on this file. Any future retarget swaps this one entry and re-runs receptor_prep, but all previously computed comparisons are invalidated and must be recomputed.
