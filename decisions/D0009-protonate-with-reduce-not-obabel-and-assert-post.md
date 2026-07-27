---
id: D0009
title: Protonate with reduce, not obabel, and assert post-conditions
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/receptor_prep.py
  - config/receptor.yaml
evidence:
  - 'obabel -p 7.4 on 6VAJ dropped 28 of 150 residues, including Cys113'
  - 'it also renamed residues (LYS A 6 -> TRP) and invented chains A I J M N'
  - 'reduce -BUILD preserves all 150 residues, chain A only, all names identical'
  - 'obabel -xr without -p is fine for pdbqt: Cys113 SG at exact raw coordinates'
  - 'guard verified against the corrupt output: catches all four failure modes'
runbook: docs/runbooks/receptor_selection.md
---

## Context

The first `receptor_prep.py` protonated with `obabel -p 7.4`. That is not a
receptor-preparation tool: on 6VAJ it renumbered residues from 1, **renamed**
them, invented four extra chain IDs, and silently dropped 28 of 150 residues —
including **Cys113**, the catalytic residue that T_3 and T_4 both target.

The output still looked like a protein. It would have docked without error and
produced plausible, meaningless scores. It was hash-pinned into a manifest and
described in published documentation before anyone noticed.

It was caught only incidentally: the DiffSBDD smoke test tried to look up
pocket residue 59 and got a `KeyError`.

The deeper failure was in verification, not tool choice. The module checked that
Cys113 existed in the **input** and never re-checked the **output**, and the
atom counts it logged were taken before protonation — so the log looked healthy.

## Decision

Protonate with `reduce -BUILD` (AmberTools; Word et al. 1999), which preserves
the heavy-atom record. Convert to PDBQT with `obabel -xr` and **no** `-p` flag,
which is safe.

Add post-conditions checked on the **outputs**, not the inputs: residue count
preserved, chain set unchanged, no residue renamed, and the catalytic residue
present with the right name in both the PDB and the PDBQT. Failure raises and
refuses to produce a receptor.

## Consequences

M1 had to be redone; the earlier prepared receptor and its recorded hashes were
invalid. `reduce` lives in the `amber_md` env rather than `cheminf`, so the
binary is searched for rather than assumed on PATH.

The general lesson is recorded in the receptor runbook: **verify the artifact
you produced, not the one you consumed.** Any preparation step that can silently
drop part of a structure needs a post-condition, not a log line.
