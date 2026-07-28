---
id: D0032
title: MM-GBSA does not rescue the ranking — and a negative verdict needs as much power as a positive one
date: 2026-07-28
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/run_mmgbsa_gate.py
  - shared/enrichment_gate.py
  - config/gates.yaml
  - decisions/D0031-class-matched-decoys-remove-the-apparent-covalent-enrichment.md
evidence:
  - 'MM-GBSA on the class-matched gate set: ROC-AUC 0.140, CI [0.060, 0.240], EF1% 0.0, BEDROC 0.000'
  - 'docking on the SAME 51 ligands: ROC-AUC 0.440, EF1% 0.0 — MM-GBSA is 0.300 WORSE'
  - 'Sulfopin ranks 44 of 51 by MM-GBSA dG and 29 of 51 by docking affinity'
  - 'Sulfopin dG -15.76; decoy dG median -21.90, best -46.45'
  - 'dG vs heavy-atom count Spearman -0.291; affinity vs heavy-atom count -0.409, so the size artefact does NOT explain it'
  - 'coverage: 51 of 83 ligands scored — all 32 naphthoquinone_c2 failed on the sp2 junction gap'
  - 'only 1 active survives, so the result is UNDERPOWERED in both directions'
  - 'the gate initially graded this FAIL, because its power floor did not govern the FAIL branch'
---

# MM-GBSA does not rescue the ranking

## What was asked

D0031 left the build with no working discriminator: docking is at chance
against class-matched decoys. The natural reading of a fidelity ladder is
that MM-GBSA *refines* a docking shortlist — but if the rung below is at
chance, refinement is the wrong frame. MM-GBSA is an **independent**
estimator, so the question was whether it can do what docking could not,
tested on the same actives, the same class-matched decoys and the same
graded gate.

## What came back

| metric | ROC-AUC | EF1% | Sulfopin's rank |
|---|---|---|---|
| docking `affinity_kcal` | 0.440 | 0.0 | 29 of 51 |
| MM-GBSA `dG_kcal` | **0.140** | 0.0 | **44 of 51** |

MM-GBSA is **0.300 ROC-AUC worse than docking** on the identical
ligands. Sulfopin — a verified nanomolar covalent Pin1 inhibitor — is
ranked below 43 of 50 property-matched chloroacetamide decoys, with a dG
of -15.76 against a decoy median of -21.90.

**The obvious artefact does not explain it.** MM-GBSA is notorious for
tracking molecular size, but here dG correlates with heavy-atom count at
Spearman -0.291 while docking's affinity correlates at -0.409. The more
size-biased metric is the one that did better.

## The verdict is UNDERPOWERED, and that matters

Only **one** active survives: the sp2 junction gap killed all 32
naphthoquinone systems, taking Juglone with them. One active cannot
support a negative claim, so the honest verdict is UNDERPOWERED in both
directions and the below-chance point estimate is *reported, not
claimed*.

## The gate defect this exposed

The gate graded it **FAIL** — "demonstrably anti-correlated with known
actives". It should not have been able to.

`_verdict()` carried the docstring "Power is checked BEFORE the point
estimates", and did not do that: the FAIL branch sat **above** the
chemotype floor. So a damning verdict could be returned from evidence the
same function would refuse to call STRONG.

With one active the interval is not what it appears to be. ROC-AUC
reduces to the fraction of decoys that one molecule beats, and the
bootstrap resamples that same active every iteration — so the CI
describes only decoy sampling. It looks *tight* precisely because the
quantity that matters, variation between actives, has no way to enter it.

**Decision: the power floor governs FAIL too.** `min_actives_for_verdict`
(3) joins `min_independent_chemotypes_for_verdict` (6), and below either
the verdict is UNDERPOWERED whichever way the estimate points. A
below-chance estimate is still printed, with a note that it *would* grade
FAIL given adequate power. FAIL remains reachable — a test pins that,
because a floor that made it unreachable would be its own defect.

This is the same asymmetry the project has guarded against in one
direction only. D0012 established that scarce actives must not manufacture
a confident PASS. The mirror image is just as wrong: scarce actives must
not manufacture a confident *rejection*.

## What this does and does not license

**Does:** MM-GBSA is not the discriminator this build needs. Nothing
justifies promoting it over docking, and the fidelity ladder's top rung
does not fix the bottom rung's problem.

**Does not:** conclude that MM-GBSA fails on Pin1. One active, one
chemotype, single-structure minimisation, no ensemble and no entropy.
The measurement is too weak to convict.

**The most useful reading** is that two methods of very different cost
both fail to separate one known active from same-chemotype decoys on this
receptor. That points at the pocket — shallow, solvent-exposed, the
regime where structure-based scoring is weakest — rather than at either
scoring function.

## What would actually change the answer

1. **Extend the junction to sp2 attachment.** It costs 32 of 83 gate
   ligands and 5 of T_4's 7 chemotypes. Restoring Juglone alone doubles
   the actives and the chemotypes.
2. **More actives.** Still the binding constraint, and still a literature
   problem rather than a compute one.
3. **Ensemble MM-GBSA over short explicit-solvent MD.** The single-structure
   protocol has no entropy and no error bar. OpenMM 8.5.2 with CUDA is
   installed and reads our existing prmtop files directly.
