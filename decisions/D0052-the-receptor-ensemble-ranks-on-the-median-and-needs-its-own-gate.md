---
id: D0052
title: The T_2 receptor ensemble reports per-receptor scores, ranks on the MEDIAN, and must pass its own enrichment gate first
date: 2026-08-01
status: proposed
approach: t2
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/noncovalent_dock_run.py
  - shared/receptor_prep.py
  - shared/rank_shortlist.py
  - config/approaches/t2_atra_crem.yaml
evidence:
  - 'ensemble decided in #6 item 6: 6VAJ + 3IKG + 3IKD + 9INR, all five seeds into all four'
  - 'best-across-receptors is a maximum over 4 correlated draws; its upward bias scales with the width of a ligand mean rotatable bonds: liu_2024_c3 10.65, du_xu 4.81, guo_pfizer 5.15'
  - 'liu_2024_c3 also carries mean 43.7 heavy atoms against 28.5-28.6 for du_xu/guo_pfizer'
  - 'D0049 removed a size bias of rho -0.617 (T_1) / -0.695 (T_3) from the ranking'
  - "gate verdicts are measured per receptor and box: D0016's ROC-AUC 0.535 is a 6VAJ number"
  - 'D0046: cross-docking into 6VAJ recovers 5.0% of poses vs 16.2% self-docking'
  - 'measured throughput on this box: du_xu 9,736 in 9h56m (~980/h), guo_pfizer 8,670 in 9h04m (~956/h), atra 1,882 in 1h16m (~1,486/h), liu_2024_c3 16,806 still running at 22.7h (<=740/h)'
---

# Pre-registering how four receptors become one number

## Context

#6 item 6 decided the ensemble and correctly demanded that the **combination
rule be recorded before anyone looks at results** — the D0045 discipline. This
record does that, and adds two things the original task list did not
anticipate, both found while implementing the multi-seed GUI.

## Decision 1 — report all four, rank on the MEDIAN

Every candidate carries `vina_affinity_<pdb>` for all four receptors, plus:

* `vina_affinity_ensemble_median` — **the rank metric**
* `vina_affinity_ensemble_best` — carried as a labelled column, never sorted on
* `vina_affinity_ensemble_argbest` — which receptor produced the best score

**Why not best-across-receptors, the usual choice.** Best-across is a *maximum
over four correlated draws*, so its upward bias grows with the width of a
ligand's score distribution — and that width scales with conformational
flexibility. Our pools differ enormously on exactly that axis:

| seed | mean rotatable bonds | mean heavy atoms |
|---|---|---|
| du_xu | 4.81 | 28.6 |
| guo_pfizer | 5.15 | 28.5 |
| **liu_2024_c3** | **10.65** | **43.7** |

Ranking on best-across would therefore hand a systematic advantage to the
flexible, larger pool for a reason that has nothing to do with binding —
reintroducing the class of artefact D0049 has just removed, one level up, and
doing it *invisibly* because "we docked into an ensemble" sounds like a
refinement. The median is robust to one receptor being a poor fit for a given
ligand, which is the actual thing an ensemble is for.

Best-across is still computed, because "which receptor does this ligand
prefer" is a real question — it is simply not the sort key.

## Decision 2 — the ensemble metric does not inherit 6VAJ's gate

`rank_shortlist.attach_gate` looks a verdict up by (stratum, metric). An
ensemble metric is a **new metric**, and the enrichment gate has never been run
on it. Attaching D0016/D0041's verdict — measured on 6VAJ with
`box_expanded.json` — to a median-over-four-receptors score would report a
number about one receptor as though it described the ensemble.

That is precisely the defect D0051 fixed one level down, where a verdict the
check did not anticipate defaulted to validating the ranking. **Run the
enrichment gate on the ensemble metric before it ranks anything.** If the gate
has no entry for it, `attach_gate` already returns `UNGATED`, which D0051 makes
non-validating — so the failure mode is now safe, but shipping an ungated
ranking without saying so is still not acceptable.

This is also the honest test of the ensemble's value: if ensemble cross-docking
does not beat 6VAJ's 5.0% pose recovery (D0046), the ensemble has not bought
what it was adopted to buy, and that is a result worth recording.

## Decision 3 — the ensemble score is size-decorrelated like any other

D0049 applies unchanged: rank on the size-decorrelated residual of the
ensemble median, not on the raw median.

## What the original task list underestimated

`shared/noncovalent_dock_run.py` is single-receptor **by construction**, not by
configuration:

| | today | needed |
|---|---|---|
| receptor | `RECEPTOR_PDBQT` module constant | parameter |
| box | `BOX_EXPANDED` module constant | per receptor, from its own reference ligand |
| `run_vina_gpu` | `(ligand_dir, out_dir, gpu)` | takes a receptor + box |
| pose directory | `poses_{LIGAND_PREP_TAG}` | must also carry a receptor tag |
| `collect_modes` | one row per `candidate_id` | keyed on (candidate, receptor) |
| frame | one `vina_affinity` | four + three ensemble columns |

**The pose-directory collision is the dangerous one.** Four receptors writing
`poses_ph7.4/` would overwrite each other's output, and `collect_modes` would
parse whatever was written last while every manifest recorded a successful run.
It is the same defect the ligand prep cache already carries a tag to prevent —
a path keyed on less than its inputs — and it must be fixed the same way, with
a new directory rather than a refreshed one.

## Cost, measured rather than estimated

#6's task list says ~120 GPU-hours, ~24 h on five cards. Measured throughput on
this box makes that optimistic, because it assumes a uniform rate and the pools
are not uniform:

| seed | molecules | measured / projected per receptor |
|---|---|---|
| atra | 1,882 | 1.3 h |
| du_xu | 9,736 | 9.9 h |
| guo_pfizer | 8,670 | 9.1 h |
| potter_astex | 7,376 | ~7.5 h |
| **liu_2024_c3** | **16,806** | **~38 h** (10.65 rotatable bonds) |

≈ **66 GPU-hours per receptor**, so ≈ **265 GPU-hours** for four. On five free
cards that is 2-3 days wall-clock, and **liu_2024_c3 is the critical path**: it
alone is ~152 GPU-hours and parallelises only across the four receptors.

Worth deciding whether the ensemble runs on all five seeds or on the shortlists
only. A shortlist-only ensemble (~125 molecules per seed) costs minutes and
answers "does the ensemble change who is on top", which may be the question
actually worth the compute given that the ranking is unvalidated either way.
