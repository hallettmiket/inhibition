---
id: D0003
title: Retain and report unknown heteroatoms rather than stripping by default
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/receptor_prep.py
evidence:
  - '6VAJ heteroatoms: 132 HOH, 31 PG4, 16 QT7, 10 SO4'
  - 'PG4 (cryoprotectant) was retained by the first run and reported, then added to the strip list'
  - 'PG4 nearest atom was 22.65 A from the box centre, so no result was affected'
  - 'final prep: 1215 protein atoms kept, 0 unrecognized retained'
runbook: docs/runbooks/receptor_selection.md
---

## Context
Stripping every HETATM would silently remove structural cofactors and metals; keeping everything that is not water would silently leave cryoprotectants occupying pocket volume. Both failures are invisible - docking succeeds either way.

## Decision
receptor_prep.py strips only an explicit list of solvent/buffer/cryoprotectant codes. Anything unrecognized is KEPT and counted in prep_log.json as other_het_atoms_retained. A non-zero count is a signal to go look.

## Consequences
Preparing a new target requires a heteroatom inventory pass rather than a blind run; the runbook makes that a step. The mechanism proved itself immediately by catching PG4.
