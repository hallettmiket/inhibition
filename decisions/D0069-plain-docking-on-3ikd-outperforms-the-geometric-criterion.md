---
id: D0069
title: Plain docking on 3IKD separates covalent binders from measured inactives, and outperforms the geometric criterion
date: 2026-08-06
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - docs/ranking_rationale.md
  - decisions/D0065-warhead-presentation-geometry-ranks-covalent-candidates.md
  - decisions/D0041-docking-enrichment-is-not-significant.md
evidence:
  - 'unbiased AutoDock-GPU on 3IKD, rigid receptor, no reactive potential, no reactive typing, same 26 A box centre (14.036 7.186 -2.108), 2000 runs'
  - 'same 75 molecules as the geometric validation: 15 crystallographic Cys113 positives, 60 warhead-matched AID 504891 measured inactives'
  - 'PLAIN docking best_dg: AUC 0.783, p=0.0004 (positives median -6.66 kcal/mol, negatives -5.71)'
  - 'per class plain: chloroacetamide AUC 0.900 p=0.0004; snar 0.850 p=0.062; michael 0.725 p=0.088'
  - 'reactive docking best_dg on the same molecules: AUC 0.824 p=0.0001'
  - 'reactive docking geometry (enrichment, 2000 runs): AUC 0.672 p=0.0207'
  - 'plain vs reactive best_dg: Spearman rho +0.880 -- the reactive potential adds little to the energy'
  - 'plain -dG vs enrichment: rho +0.211 -- energy and geometry are largely independent signals'
---

# Plain docking on the corrected receptor works, and beats the geometry

## The result

Asked to separate 15 crystallographically-verified Cys113 binders from 60
warhead-matched **measured** inactives — the same 75 molecules the geometric
framework was validated on — **ordinary AutoDock docking with no reactive
potential and no reactive typing** gives:

| method | AUC | p |
|---|---|---|
| **plain docking (`best_dg`)** | **0.783** | 0.0004 |
| reactive docking (`best_dg`) | 0.824 | 0.0001 |
| reactive docking, geometry (`enrichment`) | 0.672 | 0.0207 |

Per class, plain docking: **chloroacetamide 0.900** (p = 0.0004), SNAr 0.850
(p = 0.062), Michael 0.725 (p = 0.088).

**Plain docking outperforms the geometric criterion this branch was built to
provide** (0.783 vs 0.672), and the reactive potential adds little on top:
plain and reactive `best_dg` correlate at ρ = **0.880**.

## Why this does not contradict D0041

D0041 measured **Vina** enrichment on **6VAJ** against **property-matched
computational decoys**, and found ROC-AUC 0.599, CI [0.311, 0.874], EF1% 0.0.

Every one of those three differs here:

| | D0041 | this |
|---|---|---|
| receptor | 6VAJ (sulfopin-induced fit) | **3IKD** (chemist-prepared) |
| program | Vina | AutoDock4 / AutoDock-GPU |
| negatives | property-matched decoys, *assumed* inactive | warhead-matched, **measured** inactive |

**The receptor is the likely explanation**, and it is not a new hypothesis: D0059
already measured 3IKD giving 2.6× the pose recovery of 6VAJ (best-of-9 15.9% →
41.5%). This says the same change also restores enrichment. The docking score
was not broken; it was being asked about a pocket induced-fit around a different
ligand.

## What this means for the ranking framework

`docs/ranking_rationale.md` opens with "rank on whether a molecule can orient to
form the bond, **not** on how good the bond would be", and justifies discarding
affinity by five measured failures. **On the corrected receptor that
justification does not hold.** The affinity estimate separates actives from
inactives better than the geometry does.

This is not a claim that the geometry is worthless. It carries real, significant
signal (D0065, and chloroacetamide AUC 0.756 at convergence), it is **largely
independent** of the energy (ρ = 0.211), and the pre-registered stage-4 rule
combining them beat either alone for chloroacetamide. Two independent signals
that each work is a better position than one.

It is a claim that **the framework's founding premise — that affinity carries no
signal here — was a property of 6VAJ and did not survive the receptor change**,
and that nobody re-tested it after D0059 invalidated the 6VAJ measurements. That
omission is this record's real subject.

## What must happen before this is acted on

1. **Re-run the enrichment gate on 3IKD.** D0041's verdict is formally invalid
   under D0059 and has never been replaced. This result is 75 molecules; the gate
   has the machinery for a proper measurement.
2. **The negatives are still HTS inactives.** Weak per compound, and shared with
   the geometric validation, so the *comparison* is sound even though the
   absolute level is uncertain.
3. **15 positives.** Same ceiling as everywhere else in this project.
4. **Do not discard the geometry.** It is independent signal, and the combined
   rule already outperformed both for the one class with enough positives to
   test.

## The pattern worth recording

An assumption was measured once, on a receptor later shown to be wrong, and then
inherited by everything built afterwards without being re-measured. The receptor
change (D0059) invalidated the measurement that justified discarding affinity,
and the invalidation was recorded — `attach_gate` even fails closed on it — but
the *conclusion* it supported was never revisited.

This is the project's signature defect operating on a belief rather than a value:
carried forward by inheritance rather than by identity, populated and plausible,
and wrong.
