---
id: D0024
title: Regiochemistry re-decided on adduct-form poses — bdhi_c4 and naphthoquinone_benzo
date: 2026-07-27
status: accepted
approach: t4
decided_by: '@mhallet'
origin: implementation
supersedes: [D0021]
superseded_by: null
affects:
  - approaches/t4_combinatorial/05_regiochemistry_comparison.py
  - config/approaches/t4_combinatorial.yaml
  - data/reference/warhead_classes_4.csv
evidence:
  - 'bdhi_c4 median -3.79 vs bdhi_c5 -2.87; paired median difference -0.58 kcal/mol, Wilcoxon p = 0.0035'
  - 'bdhi pose success no longer separates the arms: 18 vs 29 discordant pairs, McNemar p = 0.14'
  - 'bdhi no-pose collapsed once the bromine was removed: c4 50% -> 14%, c5 37% -> 20%'
  - 'naphthoquinone_benzo over c2: 69 discordant pairs vs 6, McNemar p = 1.16e-14, unchanged from the pre-redock run'
  - 'naphthoquinone_c2 still fails to pose for 96% of R-groups'
  - 'convergence control passed: the three SN2 acetamides give identical best (-8.16) and median (-5.02)'
runbook: null
---

## Context

D0021 decided both regiochemistries on poses that turned out to have been docked
in the pre-reaction form (D0022). Its BDHI call was withdrawn; its naphthoquinone
call was argued to be unaffected. The re-dock on adduct-form ligands settles
both.

## Decision

**Carry `bdhi_c4` and `naphthoquinone_benzo`.**

**BDHI reverses.** `bdhi_c4` now leads on median (−3.79 vs −2.87) where it
previously trailed (+0.01 vs −2.22). The reversal came through the loser: C4's
retained bromine had blocked half its poses, and that read as a geometric
failure of the C4 attachment when it was really a failure to remove a leaving
group. Verdict **UNDERPOWERED** — the affinity difference is real but modest
(p = 0.0035, rank-biserial −0.247) and the pose-success endpoint no longer
separates the arms (p = 0.14), because with the artifact gone most poses succeed
in both. This is a weak preference, and worth revisiting if MM-GBSA disagrees.

**Naphthoquinone is unchanged, STRONG.** Benzo over C2, 69 discordant pairs
against 6, p = 1.2e-14; C2 still finds no pose for 96% of R-groups. D0021's
withdrawal notice predicted exactly this on the grounds that Michael acceptors
carry no leaving group. The prediction holding is mild independent evidence the
diagnosis in D0022 was correct.

## Consequences

The BDHI attachment is now a *weak* call rather than a decided one. It should be
carried as provisional and revisited at MM-GBSA, which models the bonded complex
explicitly rather than inferring geometry from a docking constraint.

The general lesson is recorded rather than left implicit: **a downstream
comparison inherits every defect of the poses it reads.** Both D0021 calls were
computed correctly, tested appropriately, and reported with calibrated
confidence — and one of them was still backwards, because the input was wrong in
a way no statistic on that input could reveal. The controls that caught it were
chemical (valence, interatomic distance), not statistical.
