---
id: D0006
title: Tier warhead structures by evidence and default to VERIFIED only
date: 2026-07-27
status: accepted
approach: t4
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - data/reference/warhead_classes_2.csv
  - shared/warhead_library.py
evidence:
  - '4 tiers: VERIFIED, VERIFIED_CLASS_ONLY, NEEDS_DESIGN, UNVERIFIED'
  - 'enumerable now: chloroacetamide, sulfamate_acetamide, sulfonate_acetamide'
  - 'BDHI is VERIFIED_CLASS_ONLY (PubChem CID 21983498) - usable as a window anchor, not enumerable'
  - 'naphthoquinone is NEEDS_DESIGN - juglone and KPT-6566 are intact quinones yielding no attachable fragment'
  - 'prior run: 6 of 16 warhead classes collapsed to inert amides once attached'
runbook: docs/runbooks/resolving_unverified_structures.md
---

## Context
The warhead set must be data so the choreography can go wide later by adding rows. But breadth without provenance is how a library ends up half-dead, and a status field nothing enforces is just a comment.

## Decision
warhead_classes_*.csv carries a structure_status per class; warhead_library.enumerable() defaults to VERIFIED only and logs a warning naming every class a caller widens to. window_anchor_classes() separately admits VERIFIED_CLASS_ONLY, because a window needs the chemotype, not the compound.

## Consequences
T_4 currently has 3 enumerable classes. BDHI and naphthoquinone are blocked on attachment-regiochemistry design, which is a chemist's call and not a coding task. Going wide later is a CSV edit plus an explicit widening argument.
