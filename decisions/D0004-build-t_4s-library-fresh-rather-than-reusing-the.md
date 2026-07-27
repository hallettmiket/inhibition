---
id: D0004
title: Build T_4's library fresh rather than reusing the prior run
date: 2026-07-27
status: accepted
approach: t4
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/gates.yaml
  - data/reference/warhead_classes_2.csv
evidence:
  - 'prior library: /data/lab_vm/refined/pin1_acr_screen, owner hemam, built 2026-07-21'
  - 'reference set assembled 2026-07-26, five days later'
  - 'so the prior library could not have been grounded in the reference set'
  - 'prior library was 7104 rows = 16 warhead classes x 444 R-groups'
runbook: null
---

## Context
A 7,104-member combinatorial library from a prior run existed and was of good quality (unique SMILES, core-verified, built on the graph engine). Reusing it would have saved real work.

## Decision
Do not reuse it. Build T_4's warhead and R-group libraries fresh, grounded in the frozen reference set. The inherited library_size: 7104 was removed from gates.yaml; library size is now a derived, pinned output.

## Consequences
Loses the prior work but gains a library that can actually satisfy the spec's grounding requirement, and removes a dependency on another group's directory tree (owner hemam, group ssmd-u-biodatsci-otherslab) that could move or vanish.
