---
id: D0068
title: Enrichment depends on search effort, and at convergence the docking energy discriminates better than the geometry
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
  - shared/nac_criterion.py
  - scripts/nac_rank.py
evidence:
  - 'refine of the top 300 at 2000 runs: median enrichment fell from 2.91x (200 runs) to 0.96x on the SAME molecules; Spearman rho 0.364; top-50 overlap 23/50'
  - '15 crystallographic positives, NOT selected on score, fell the same way: 2.27x -> 0.97x. So this is not selection bias.'
  - 'n_poses equals nrun exactly at both settings (200 and 2000), so the pose population is all runs and nothing is being clustered away'
  - 'with 10x the search: median S...C distance 3.34 -> 3.79 A, median angle 44.6 -> 57.5 deg, viable fraction 0.235 -> 0.077'
  - 'same 75 molecules (15 pos / 60 neg): enrichment AUC 0.787 (p=0.0003) at 200 runs vs 0.672 (p=0.0207) at 2000 runs'
  - 'per class at 2000 runs: chloroacetamide AUC 0.756 p=0.016; michael 0.750 p=0.065; snar 0.575 p=0.388'
  - 'best_dg alone on the same 75 molecules at 2000 runs: AUC 0.824, p=0.0001 -- better than the geometry'
---

# Enrichment is a function of search effort, and energy beats it at convergence

## What was found

The shortlist refinement was meant to sharpen the top 300 by docking them 10×
harder. Instead it revealed that **the metric does not converge**.

The same 300 molecules scored **2.91× median at 200 runs and 0.96× at 2,000**.
Spearman ρ between the two was 0.364 and only 23 of the top 50 survived.

The obvious reading is winner's curse — select the top of a noisy estimator and
it regresses. **That is not what happened.** The 15 crystallographic positives,
which are chosen by crystallography and never by score, fell exactly the same
way: **2.27× → 0.97×**.

## The mechanism

`n_poses` equals `nrun` exactly at both settings, so nothing is being clustered
away and the population really is every run. What changes is *where the search
ends up*:

| on the same molecules | 200 runs | 2,000 runs |
|---|---|---|
| median S···C distance | 3.34 Å | **3.79 Å** |
| median approach angle | 44.6° | **57.5°** |
| viable fraction | 0.235 | **0.077** |

**More search effort finds lower-energy poses that are less reaction-competent.**
The reactive well (D0064: `r_eq` 3.2 Å, `eps` 1.0) is not deep enough to dominate
the true energy minimum. At low `nrun` the search has not converged, and
under-optimised poses sit near the well by accident; at high `nrun` AutoDock finds
the actual optimum, which is not a near-attack conformation.

So **the viable fraction is a property of the search, not only of the molecule.**
Absolute enrichments are not converged quantities and must never be quoted as if
they were.

## Does the discrimination survive? Partly — and the energy does better

Measured on the same 75 molecules (15 crystallographic positives, 60
warhead-matched measured inactives) at both settings:

| | AUC | p |
|---|---|---|
| enrichment, 200 runs | 0.787 | 0.0003 |
| **enrichment, 2,000 runs** | **0.672** | 0.0207 |
| **`best_dg` alone, 2,000 runs** | **0.824** | **0.0001** |

Per class at convergence: chloroacetamide **0.756** (p = 0.016), Michael 0.750
(p = 0.065), SNAr 0.575 (n.s.).

**The geometric signal is real and it weakens under convergence, and the plain
docking energy from the same runs discriminates better than it.**

## What this does and does not overturn

**Does not overturn:** the validation itself. D0065's AUCs compared positives and
negatives **at the same `nrun`**, so the comparison was like-for-like and the
robustness run (10 disjoint negative draws, chloroacetamide 0.908) stands as a
statement about 200-run enrichment. Geometry does carry signal, at every effort
tested.

**Does overturn:** the framework's central claim as written. `ranking_rationale`
says to rank on whether a molecule can orient to form the bond **rather than** on
affinity, on the strength of five measurements showing the docking score carries
no signal. On this receptor, with this program, asked to separate molecules, the
energy carries **more** signal than the geometry. The premise that licensed
discarding it does not hold here.

**Does invalidate:** the production shortlist. The top 300 selected at 200 runs
regress to a median of 0.96×, and 23 of 50 top-ranked survive re-measurement. A
shortlist drawn from the 200-run screen is not usable.

## Consequences

1. **No shortlist is issued from the 200-run screen.** It ranked 5,769 candidates
   on a quantity that does not converge.
2. **`nrun` becomes part of the metric's definition**, not a tuning knob. Any
   enrichment must be quoted with the run count that produced it.
3. **The outstanding control is now the important one**: these energies come from
   runs whose *sampling* was biased toward the warhead–sulfur contact, so an
   **unbiased AutoDock run** is required before `best_dg` can be called an
   affinity signal rather than a proximity artefact. If it survives that, the
   covalent arm should rank on energy, or on energy gated by geometry (stage 4's
   C2 rule already beat enrichment alone for chloroacetamide, AUC 0.953 vs
   0.908 at 200 runs).
4. **D0065 is not withdrawn but is now conditional** — its numbers describe
   200-run enrichment, and the interpretation section must carry this record's
   caveat.
