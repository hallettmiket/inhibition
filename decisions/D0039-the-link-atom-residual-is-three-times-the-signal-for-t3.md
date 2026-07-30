---
id: D0039
title: The link-atom residual is measured, and for T_3 it is three times the interaction energy
date: 2026-07-30
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/merge_ensemble_dg.py
  - integration/app/app.py
  - decisions/D0037-the-junction-dihedral-was-the-sp3-analogue-and-dG-was-not-an-interaction-energy.md
evidence:
  - 'median |residual|/|interaction|: T_1 0.00, T_2 0.00, T_4 0.16, T_3 2.88'
  - 'T_3 interaction -5.05 kcal/mol, internal residual -14.78 kcal/mol (n=25 candidates)'
  - 'T_4 interaction -12.95 kcal/mol, internal residual +0.77 kcal/mol (n=27)'
  - 'T_1 -20.16 / -0.00 (n=24); T_2 -21.51 / -0.00 (n=24)'
  - 'all 182 ensemble results carried the split; the frames did not, because the merge predated it'
  - 'no shortlist, gate or rank consumes any dG_ensemble column — grep over approaches/, shared/, scripts/ returns nothing'
---

# The contamination has a size, and it is not small

## What D0037 established and what it left open

D0037 found that the ensemble `dG` was not an interaction energy. Single-
trajectory MM-GBSA is *defined* by its bonded terms cancelling between the three
legs; the covalent decomposition cuts the Cys113 SG–C bond and caps both sides
with hydrogen, and the caps do not cancel. The scorer was split into an
interaction part and an internal residual on that basis.

What D0037 did not do is say **how big** the residual is. The split existed in
`shared/mmgbsa_ensemble.py` and in all 182 per-candidate results, but
`scripts/merge_ensemble_dg.py` predated it, so no frame carried the new columns
and no consumer could see them. For three days every reported ensemble `dG` was
the full potential difference, and nobody — including me — knew whether the
contamination was 1% or 300%.

## The measurement

| approach | covalent | interaction | internal residual | median \|residual\|/\|interaction\| |
|---|---|---|---|---|
| T_1 | no | −20.16 | −0.00 | **0.00** |
| T_2 | no | −21.51 | −0.00 | **0.00** |
| T_4 | yes | −12.95 | +0.77 | 0.16 |
| T_3 | yes | −5.05 | −14.78 | **2.88** |

Medians over the candidates carrying an ensemble result (24, 24, 27, 25).

Two things follow immediately.

**The non-covalent approaches validate the method.** T_1 and T_2 have no link
atom, and their residual is 0.00 — not small, *zero* to the precision reported.
That is what single-trajectory MM-GBSA requires, and it is the control that says
the split is measuring what it claims rather than manufacturing a difference.

**T_3's ensemble dG is mostly artefact.** A residual of −14.78 against an
interaction energy of −5.05 means roughly three quarters of the number is the
hydrogen caps. Anyone reading `dG_ensemble_kcal` for a T_3 candidate was reading
a quantity dominated by a modelling device.

## Why T_3 and not T_4

Both are covalent and both use the same decomposition, so a 2.88 against a 0.16
is not explained by the method alone. The obvious candidate is the junction
parameters: T_3 and T_4 draw on different warhead classes, and the `cc`/`cd`
dihedral introduced in D0037 — the one GAFF2 provides no generic for, where the
2.430/2-fold/0° choice is mine rather than the force field's — governs the
naphthoquinone chemotypes.

**This is a hypothesis, not a finding.** It is stated here so the sensitivity
check that D0037 already called for has a specific prediction to test: if the
`cc`/`cd` term is responsible, varying it should move T_3's residual and leave
T_1/T_2 at zero and T_4 roughly where it is. Until that runs, the honest
statement is that T_3's residual is large and its cause is unestablished.

## What this does and does not change

**Does not change any ranking.** No shortlist, gate or rank consumes any
`dG_ensemble` column — verified by grep across `approaches/`, `shared/` and
`scripts/`, which returns nothing. The ensemble values have always been
descriptive. This is the second time that has made a correction cheap (D0038 was
the first), and it is not luck: keeping an unvalidated quantity out of the gate
is what makes it survivable when the quantity turns out to be wrong.

**Does change what the interface shows.** The frames now carry the interaction
and residual columns, and the dossier reads the ratio per approach and either
tells the reader to use the interaction column or states that the terms cancel.
The threshold is 0.10. The warning quotes the measured ratio rather than
describing it, because "there is a residual" is compatible with both 0.16 and
2.88 and the reader cannot act on the sentence without the number.

## The general lesson

The split was implemented, tested and correct, and it reached nothing. A
correction that stops at the module that computes it is not a correction — it is
a corrected value sitting beside the uncorrected one that everything actually
reads. The gap here was one dictionary in a merge script, and it survived three
days precisely because both the module and its tests were right.

What would have caught it sooner: asking, after any fix, *which file does a
reader open, and does the fix appear in it?* The scorer was never the file
anyone opened.
