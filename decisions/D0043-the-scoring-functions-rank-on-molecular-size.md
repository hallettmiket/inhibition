---
id: D0043
title: The shortlists are size-selected — every scoring function here ranks partly on heavy-atom count
date: 2026-07-30
status: accepted
approach: shared
decided_by: '@mhallet'
origin: adversary
supersedes: []
superseded_by: null
affects:
  - shared/pocket_size.py
  - approaches/t3_reinvent/01_generate.py
  - approaches/t2_atra_crem/01_generate.py
  - decisions/D0041-the-first-verdict-docking-does-not-demonstrably-enrich.md
evidence:
  - 'CORRECTED 2026-07-31 — on each approach RANK metric, all lower-is-better:'
  - 'spearman(heavy_atoms, vina_affinity) T_1 = -0.617 (n=3233): bigger scores better'
  - 'spearman(heavy_atoms, affinity_kcal) T_3 = -0.479 (n=4080): bigger scores better'
  - 'spearman(heavy_atoms, vina_affinity) T_2 = -0.230 (n=1882): bigger scores better'
  - 'spearman(heavy_atoms, affinity_kcal) T_4 = +0.181 (n=1683): bigger scores WORSE'
  - 'SUPERSEDED: the original +0.745 T_3 / +0.305 T_4 used cnn_affinity, which is not the rank metric'
  - 'T_3 generated median 25 heavy atoms; T_3 SHORTLIST median 39, max 51'
  - 'T_3 scaffold is 12 heavy atoms, so shortlisted R-groups have median 27 — over twice the scaffold'
  - 'ligand_efficiency over-corrects: spearman(heavy_atoms, LE) T_1 = -0.938'
  - 'the 55-heavy-atom pocket ceiling removes 48/4803 T_1, 0/1882 T_2, 2/5396 T_3, 0/1782 T_4'
  - 'pocket cavity volume on 6VAJ within 6 A of QT7: 1018 A^3'
---

# The molecules are not too big; the ranking prefers big molecules

> ## Correction, 2026-07-31
>
> **The T_3 and T_4 numbers below were computed on `cnn_affinity`, which is not
> their ranking metric.** Both rank on `affinity_kcal`; `cnn_affinity` is
> explicitly flagged uncalibrated for covalent docking (D0011) and is carried as
> advisory only. Recomputed on the metric that actually orders the shortlists:
>
> | approach | rank metric | rho | direction |
> |---|---|---|---|
> | T_1 | `vina_affinity` | −0.617 | bigger scores **better** |
> | T_3 | `affinity_kcal` | **−0.479** | bigger scores **better** |
> | T_2 | `vina_affinity` | −0.230 | bigger scores **better** |
> | T_4 | `affinity_kcal` | **+0.181** | bigger scores **worse** |
>
> All four metrics here are lower-is-better, so the signs are directly
> comparable. Two things change:
>
> **T_3's size dependence is weaker than reported** — −0.479, not the +0.745
> quoted below. The finding survives; its magnitude does not.
>
> **T_4 runs the other way.** The claim below that "all four mean the same
> thing: bigger molecules score better" is FALSE for T_4 on its real ranking
> metric: larger molecules score slightly worse. The size bias holds for T_1,
> T_2 and T_3 only.
>
> The error was mine and it is the project's recurring one — a column chosen by
> plausible name rather than by checking which one the ranker reads. It was
> caught because `rank_shortlist.rank()` refuses any metric not declared in
> `LOWER_IS_BETTER` with its direction, rather than assuming one. The guard was
> written for exactly this and it worked.
>
> Everything else below — the shortlist/generated size gap, the ligand-efficiency
> over-correction, the pocket ceiling — is unaffected: those were computed on
> heavy-atom counts and on T_1/T_2, not on the covalent metrics.

## What was reported and what I first assumed

Issue #1: *"the decorated r group is quite massive, resulting in most of the
generated molecules being above the reasonable size for a pin1 inhibitor ...
Some of these almost look like protacs."* The proposed cause was that the
enumeration runs multiple decoration cycles and keeps building on one R-group.

Neither half survives checking, and the truth is more consequential than either.

**LibInvent runs once.** `n_smiles: 20000`, one attachment point, no iteration
and no reseeding. There are no cycles to limit.

**T_3's generated molecules are not large.** Median **25** heavy atoms, p90 32.
For reference the 6VAJ co-crystal ligand QT7 is 16 and ATRA is 22. Nothing about
the generated pool is protac-like.

## Where the large molecules actually come from

The SHORTLIST, not the generator.

| | generated median | shortlist median | shortlist max |
|---|---|---|---|
| T_3 | 25 | **39** | 51 |

T_3's scaffold is 12 heavy atoms, so a shortlisted molecule carries a median
**27-atom R-group — more than twice the scaffold it decorates.** That is exactly
what was observed. It is a selection effect, and the selector is the score.

## The measurement

Spearman correlation between heavy-atom count and the ranking metric, over every
scored candidate:

| approach | metric | rho | n |
|---|---|---|---|
| **T_3** | cnn_affinity | **+0.745** | 4080 |
| T_1 | vina_affinity | **-0.617** | 3233 |
| T_4 | cnn_affinity | +0.305 | 1683 |
| T_2 | vina_affinity | -0.230 | 1882 |

The signs differ only because Vina affinity is better when more negative and
CNN affinity is better when larger. **All four mean the same thing: bigger
molecules score better.** For T_3, size accounts for roughly 55% of the rank
variance. The ranking is, to a first approximation, a size ranking.

This is not a surprising defect — more atoms make more contacts, and neither
Vina's scoring function nor gnina's CNN is normalised for that — but it had not
been measured here, and it is large.

## Why this matters more than the filter it came from

D0041 measured the non-covalent gate at **WEAK**: AUC 0.599, CI [0.311, 0.874],
**EF1% 0.0**. This is a mechanism for that result. If the score is substantially
a size ranking, then whether it enriches for binders depends on whether the
known actives happen to be larger than their property-matched decoys — and the
decoys are property-matched on molecular weight, which is precisely the axis
that would neutralise it. A size-driven score should perform near chance against
size-matched decoys. It does.

## Ligand efficiency is not the fix

The obvious correction over-corrects. `spearman(heavy_atoms,
ligand_efficiency)` is **-0.938** for T_1: LE divides by heavy-atom count so
completely that it becomes a smallness ranking rather than a size-neutral one.
Swapping one for the other trades a bias for its mirror image.

What is needed is a size-decorrelated score — for example ranking on the
residual of score against heavy atoms, or scoring within size strata. That is a
design decision, not a patch, and it is recorded here as open rather than
quietly chosen.

## The pocket ceiling that started this

Derived properly and kept, though it does almost nothing.

A grid cavity calculation on 6VAJ, restricted to within 6 A of QT7 and requiring
enclosure on 4 of 6 directions, gives **1018 A^3** — ~27 heavy atoms at typical
55% packing, ~34 at tight 70%. An earlier pass over the whole 20 A docking box
returned 2185 A^3 and 59-75 atoms; it was measuring surface grooves and open
solvent. Pin1's PPIase site is shallow, so a box centred on the ligand leaves
the pocket quickly.

Set at **55** heavy atoms — ~1.6x tight packing, clearing every known
non-peptidic binder including BJP-06-005-3 at 52, and excluding only the
peptidic macrocycles. It removes 48/4803 T_1, 0/1882 T_2, 2/5396 T_3, 0/1782
T_4.

It is worth keeping as a guard against future runaway lineages, and it is worth
saying plainly that it does not address the problem it was requested for. A
ceiling cannot fix a ranking that prefers large molecules; it can only cap how
large the preferred ones get.

## The lesson

The request was for a filter, the diagnosis pointed at the generator, and the
cause was in the scorer. Measuring the thing that was actually complained about
— the size distribution of what reached the shortlist, against what was
generated — took one query and pointed somewhere nobody was looking. A filter
built on the original assumption would have removed 2 molecules out of 5396 and
been reported as done.
