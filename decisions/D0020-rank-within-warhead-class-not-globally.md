---
id: D0020
title: T_4 ranks within warhead class, with a per-class quota, not globally
date: 2026-07-27
status: accepted
approach: t4
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - approaches/t4_combinatorial/04_rank_within_class.py
  - config/approaches/t4_combinatorial.yaml
  - tests/test_rank_within_class.py
  - integration/app/DECISIONS_TAB_SPEC.md
evidence:
  - 'each warhead class docks against a DIFFERENT covalent_lig_atom_pattern, so pose ensembles are not drawn from one distribution'
  - 'warhead heavy-atom count varies across classes (20-44 in the shortlist fixture) and a Vina-style score tracks size'
  - 'on a synthetic frame where one class has the best raw affinities, a global top-9 excludes an entire class; the quota shortlist represents all three'
  - 'D0015 fixed affinity_kcal as the rank metric after measuring it on 6 actives + 294 decoys; no equivalent measurement exists for ligand efficiency'
runbook: null
---

## Context

T_4 fixes the sulfopin core and varies the warhead and the R-group. The question
it exists to answer is therefore *which warhead chemistry works on this core* —
a comparison across chemotypes.

The obvious implementation is to sort all 1,683 docked survivors on
`affinity_kcal` and take the best. That would answer a different question badly.

gnina's affinity is comparable only among molecules docked the same way, and
these are not. Each class is docked against its own
`covalent_lig_atom_pattern`, so the search is constrained differently per class
and the resulting pose ensembles are not samples from a common distribution. On
top of that, a Vina-style score grows with heavy-atom count, and the classes
differ substantially in size — the naphthoquinones carry roughly twice the heavy
atoms of the acrylamide. A global sort would largely report which warhead is
biggest and greasiest, and would hand most of the shortlist to one chemotype.

## Decision

**Rank within warhead class; each class contributes a fixed quota (3) to the
shortlist.**

Three supporting choices:

1. **Ligand efficiency is computed but advisory.** `LE = -affinity / HAC` is the
   standard size correction and reviewers will ask for it, so it is reported.
   It is *not* the rank metric. D0015 fixed the rank metric on `affinity_kcal`
   because the enrichment gate measured it on this target (ROC-AUC 0.815, EF1%
   16.7); nothing comparable has been measured for LE here. Substituting an
   unmeasured metric because it is conventional would discard the only
   calibration the choreography has.

2. **Classes with few successful docks are flagged, not dropped.** Below 20
   docks, "best in class" is most of the class and the rank is not selective.
   Such rows carry `rank_is_selective = False` and say so in
   `shortlist_reason`.

3. **`OUTSIDE_WINDOW` classes still contribute** (D0019). Excluding them at
   ranking time would silently re-impose the kinetics filter that D0019
   deliberately relaxed — a veto reintroduced through the back door of a
   selection step.

Nothing is rejected here. This stage stamps and orders.

## Consequences

The shortlist is a designed cross-chemotype comparison, not a league table. It
is deliberately *not* the set of nine best-docking molecules, and the GUI must
not present it as such: each row carries `class_rank`, `class_n_docked`,
`class_percentile` and `shortlist_reason` so a reader can see it was chosen as
best-of-its-class rather than best-overall.

Cross-class affinity comparison remains unlicensed downstream. If the
integration phase wants to rank warhead chemistries against each other, the
evidence for that is MM-GBSA on the true covalent adduct (step 9), which models
the bonded complex explicitly, not the docking score.

The per-class quota is in `config/approaches/t4_combinatorial.yaml`, so widening
the shortlist is a config change with a manifest record, not an edit to the
ranking code.
