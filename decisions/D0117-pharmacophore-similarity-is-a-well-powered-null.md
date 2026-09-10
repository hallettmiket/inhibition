---
id: D0117
title: Pharmacophore similarity does not rank on this target, and this time the null is well powered
date: 2026-09-09
status: accepted
approach: t2
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/pharmacophore.py
  - scripts/pharmacophore_gate.py
  - scripts/rank_pharmacophore.py
  - /data/lab_vm/append_only/inhibition/00_shared_substrate/enrichment_gate.token
evidence:
  - 'AID 504891: 34 assayed actives vs 361,354 ASSAYED inactives, no decoy construction'
  - 'ROC-AUC 0.5679, 95% CI [0.473, 0.658] — includes 0.5; verdict WEAK'
  - '27 independent chemotypes against a floor of 6 — the power floor is CLEARED, so this is not UNDERPOWERED'
  - 'leave-one-chemotype-out AUC spans 0.538-0.598: no single cluster carries it'
  - 'sensitivity, model set +PiB: AUC 0.5542, CI [0.466, 0.646] — same verdict, so the model-set rule is not load-bearing'
  - 'model set and gate actives are DISJOINT: ingest recorded 0 reference binders called active in AID 504891'
  - 'EF1% 0.00 carries almost no information: chance gives 0 actives in the top 1% 71% of the time at n_actives=34'
  - 'both degree-2 pools ranked and stamped rank_validated=False'
  - 'Guo pool: 69.1% of all 30,000 molecules and 25/25 of the shortlist are nearest to ONE model molecule (Potter-Astex)'
  - 'ATRA pool: Du-Xu is nearest for 11.5% pool-wide but 15/25 of the shortlist'
  - 'size decorrelation: rho(metric, HAC) +0.308 -> -0.016 (ATRA), +0.192 -> -0.022 (Guo)'
  - 'median ph2d_max_similarity 0.101 (ATRA) vs 0.266 (Guo) — not comparable across seeds'
  - 'frames: 02_t2_atra_crem_degree2/D2_9.parquet, 02e_t2_guo_crem_degree2/D2_2.parquet'
  - 'leave-one-out recovery vs the same 361,354 inactives: median percentile 18.9% (chance 50%)'
  - '28,082 MEASURED INACTIVES outscore Sulfopin when it is held out of the model set'
  - 'only 2/15 held-out binders reach the top 1%, and they are ONE chemotype (Reddi-4d/4g, ECFP4 0.741, identical scores)'
  - 'ATRA held out lands at the 99.94th percentile; Juglone at 99.47th'
runbook: null
---

## Context

#4 Phase 2.2 asks for a scorer orthogonal to docking, and says to measure it on
a recovery task **before** it ranks anything. Two facts made now the moment:

* The T_2 degree-2 pools — 60,000 molecules — have never been ranked by
  anything. The ATRA pool's docking segfaulted on all six chunks (D0060) and
  the Guo pool was never docked. The permissions fix D0060 called for has not
  landed in the month since.
* Phase 1's measured-inactive datasets were ingested on 2026-08-05 and never
  used, because Phase 2.1 was going to push **Vina** through them and Vina
  needs the GPU.

A 2D pharmacophore scorer needs no receptor, no pose, no GPU. It is the one
level of theory that could run against those 361,354 assayed negatives today.

## What was measured

Gobbi 2D pharmacophore fingerprint; score is the maximum Tanimoto to a model
set of published Pin1 binders. **This is similarity to known actives, not a
match against a pharmacophore hypothesis** — the hypothesis form needs the
co-crystals and 3D conformers and is deferred.

The model set follows a rule **written down before it was applied** (D0045's
discipline, so the answer could not be chosen): parses, fewer than 4 backbone
amide bonds, at most 45 heavy atoms. That retains 15 binders including
Liu-2024-C3 at 41 heavy atoms, and drops four peptides and CPPs carrying 28 to
103 pharmacophore features that no ~30-heavy-atom CReM product can approach.

**The train/test split is real and was not arranged.** The model set is the
literature binders; the gate actives are AID 504891's qHTS hits.
`ingest_measured_inactives` had already recorded *"reference binders also
called ACTIVE here: 0"*. Nothing in the model set can be found by looking for
itself.

| | ROC-AUC | 95% CI | EF1% | BEDROC | verdict |
|---|---|---|---|---|---|
| declared rule (15 binders) | **0.5679** | [0.473, 0.658] | 0.00 | 0.085 | **WEAK** |
| sensitivity, +PiB (16) | 0.5542 | [0.466, 0.646] | 0.00 | 0.080 | WEAK |

## The decision

`ph2d_max_similarity` is registered as a rank metric (higher-is-better, D0115)
and **does not validate a ranking**. The verdict is WEAK; only STRONG validates
(D0051's allowlist), so both degree-2 pools carry `rank_validated = False`.

The stratum is written into the shared gate token as
`pharmacophore_aid504891`. `write_token` merges per stratum and per metric, so
the `non_covalent` and `covalent` verdicts every other approach ranks on are
untouched — verified after writing.

## Why this null is worth more than the previous ones

**It is the first one that is not underpowered.** Every earlier negative on
this target came with an interval wide enough to drive a bus through, and
D0045 / §4 of the orientation document record the covalent stratum as
UNDERPOWERED by construction — validated Pin1 chemistry is scarce, and that
scarcity is a property of the target rather than a fixable defect.

This gate had 34 actives spanning **27 independent chemotypes against a floor
of 6**, and 361,354 negatives that were *assayed rather than assumed*. Compare:

| | ROC-AUC | 95% CI | CI width | power |
|---|---|---|---|---|
| Docking enrichment (D0041) | 0.599 | [0.311, 0.874] | 0.563 | underpowered |
| **Pharmacophore similarity** | **0.568** | **[0.473, 0.658]** | **0.185** | **27 chemotypes, cleared** |

The point estimates are similar. The intervals are not. D0041 could not
distinguish "does not work" from "we cannot tell"; this can, and the answer is
that the scorer sits within noise of not enriching with the interval narrow
enough for that to mean something. Leave-one-chemotype-out spans 0.538–0.598,
so no single cluster is carrying it either.

**This also disposes of the decoy-construction defence.** D0041's null could
always be blamed on property-matched decoys being unrealistically easy or
unrealistically hard. These negatives were measured in a real assay. The
explanation is not available here.

## Two things NOT to read off this

**EF1% = 0.00 is not the damning number it looks like.** With 34 actives, the
top 1% is 3,614 molecules and chance alone puts 0.34 actives there. A random
ordering returns EF1% = 0.00 **71% of the time**, and the statistic can only
take the values 0, 2.94, 5.88 — clearing the gate's own `ef_1pct_min: 5.0`
requires 2 actives in the top 1%. At this active count EF1% has almost no
resolving power and the ROC-AUC interval is what carries the verdict. Recorded
because a bare "EF1% 0.0" beside a WEAK verdict reads as much stronger evidence
than it is.

**The result is about ligand-based SIMILARITY, not about pharmacophores.** A 3D
hypothesis derived from the Pin1 co-crystals asks a different question — where
a feature sits relative to the pocket, rather than whether the molecule carries
a similar feature arrangement to some known binder. This measurement does not
speak to it, and the module docstring says so at the point of use.

## The shortlists, and what is odd about them

Both pools are ranked within themselves and **never across**. The two seeds
sit at different absolute similarity levels for reasons that are properties of
the seed rather than of the chemistry — median `ph2d_max_similarity` is
**0.101** for ATRA (a polyene with almost no pharmacophore features) against
**0.266** for Guo, and ATRA's *best* molecule (0.393) would rank below Guo's
median. A merged ordering would be a seed sort. Same trap `config/seeds.yaml`
already records for docking scores across seeds.

Each pool's own seed is held out of its model set. Measured before this was
built: Guo's pool sits at median 0.291 to its own seed against 0.269 to any
other binder, seed nearest for **58.8%** of molecules; ATRA's at 0.053 against
0.107, nearest for 16.8%. The circularity is real for one pool and absent in
the other, which is why the hold-out was applied to both rather than to the one
where it looked necessary.

### A max over a set collapses to one molecule, and the two pools do it differently

`ph2d_nearest_active` is written to every row precisely so this is visible in a
`value_counts()` rather than requiring someone to suspect it. It should be
looked at, because the metric's name promises more than it delivers:

| | nearest model molecule, WHOLE pool | top 25 |
|---|---|---|
| ATRA | KPT-6566 17.1%, Du-Xu 11.5%, Liu-2024-C3 10.9%, Reddi-4g 10.8%, Reddi-4d 10.6% | **Du-Xu 15/25**, Potter 5, Pu 3, Juglone 2 |
| Guo | **Potter-Astex 69.1%**, EGCG 16.7%, Liu-2024-C3 6.7%, Du-Xu 4.6% | **Potter-Astex 25/25** |

**The Guo pool's ranking is a Potter-Astex similarity search wearing a set's
name.** With Guo's own seed correctly held out, Potter-Astex is simply the
nearest remaining chemotype to Guo's, and it takes 69% of the pool and *all* of
the shortlist. "Maximum Tanimoto over 14 binders" is doing the work of one
comparison. That is not wrong — it is what a max is — but a reader who takes
the metric at its name will believe the shortlist is supported by fourteen
molecules' worth of evidence when it is supported by one.

The ATRA pool spreads across the model set pool-wide (no member above 17.1%)
and then concentrates at the top anyway: **15 of its top 25 are nearest to
`Du-Xu-naphthalenecarboxamide`, against 11.5% pool-wide.**

And Du-Xu is the wrong molecule to have carrying a shortlist. It has an
explicit provenance caveat in `config/seeds.yaml` — the source paper emphasises
a thiazole series and the deposited structure has no thiazole ring, so its
Ki 6 nM is attached to a structure nobody has checked against the SI. The
caveat says it is fine as a *neighbourhood seed* but **"NOT admissible as a
gate active or a potency anchor"**. Using it as a pharmacophore model molecule
is nearer the second than the first. This does not change the verdict — the
gate is WEAK with or without any single member, and leave-one-chemotype-out
confirms no cluster carries it — but a shortlist that is 60% one unverified
molecule's neighbourhood is worth less than its position implies.

### The size decorrelation did engage, and it was needed

A Tanimoto to a fixed model set has a real size dependence: a larger molecule
has more feature pairs and more chances to match. Measured, then removed by
D0049's local-median treatment:

| pool | rho(raw metric, HAC) | rho(residual, HAC) |
|---|---|---|
| ATRA | **+0.308** | −0.016 |
| Guo | **+0.192** | −0.022 |

Both land inside D0049's |rho| ≤ 0.034 band. Ranking on the raw Tanimoto would
have been partly a heavy-atom sort, which is the same defect D0043 found in the
docking score — reached by a completely different route, which is worth
noticing on its own.

## The confirmation that removes the last objection

The AID gate could still be answered with *"those 34 qHTS hits are a different
population from the literature binders — of course a model built on one does
not rank the other."* So the scorer was also run on the reference set's **own
molecules**, which that objection cannot reach.

**Leave-one-out recovery.** Hide one of the 15 model molecules, score it
against the **361,354 measured inactives** using the remaining 14, and read off
where it lands. Held-out binder and background go through the same call, so a
difference between them cannot come from a difference in how they were scored.

| held out | score | measured inactives scoring at least as high | percentile |
|---|---:|---:|---:|
| Reddi-2023-4d | 0.400 | 160 | **0.04%** |
| Reddi-2023-4g | 0.400 | 211 | **0.06%** |
| Sulfopin | 0.267 | **28,082** | 7.77% |
| Potter-Astex | 0.234 | 50,424 | 13.95% |
| Du-Xu | 0.237 | 53,215 | 14.73% |
| Guo-Pfizer | 0.237 | 55,742 | 15.43% |
| Liu-2024-C3 | 0.234 | 60,812 | 16.83% |
| Pu-benzylguanine | 0.219 | 68,306 | 18.90% |
| Ieda-2019-(S)-2 | 0.228 | 69,932 | 19.35% |
| EGCG | 0.223 | 74,143 | 20.52% |
| Tian-chloropyrimidine | 0.184 | 149,075 | 41.25% |
| Liu-2022-ZL-Pin13 | 0.150 | 215,957 | 59.76% |
| KPT-6566 | 0.136 | 263,134 | 72.82% |
| Juglone | 0.040 | 359,453 | **99.47%** |
| ATRA | 0.029 | 361,152 | **99.94%** |

**Median 18.9% against a chance median of 50%.** So the scorer is not noise —
it carries real signal. It is simply nowhere near strong enough to rank on, and
the table says exactly how it fails:

* **28,082 molecules that were assayed and found INACTIVE score higher than
  Sulfopin**, this project's nanomolar covalent lead and the best-characterised
  binder in the set. That is the whole result in one line.
* **The only two recoveries in the top 1% are the same molecule twice.**
  Reddi-4d and Reddi-4g score **identically** (0.40048543689320387) because
  each is recovering the other: they differ only in tert-butyl vs cyclohexyl,
  ECFP4 Tanimoto **0.741**, far above the 0.4 threshold `cluster_chemotypes`
  uses — one chemotype, not two. Reported as "2 of 15" and worth **one**. This
  is precisely the failure `chemotype_ids` and leave-one-chemotype-out exist to
  prevent in the gate, appearing here in the open where it can be seen.
* **ATRA lands at the 99.94th percentile** — 361,152 of 361,354 measured
  inactives outscore it. A T_2 seed. The scorer ranks a real (weak, promiscuous)
  Pin1 binder below essentially every molecule in a screening deck. Juglone,
  another genuine binder, is at 99.47%.

The last point also explains the two pools' behaviour in §"The shortlists":
ATRA has almost no pharmacophore features, so nothing in its neighbourhood can
score, while Guo's neighbourhood collapses onto whichever peer is nearest.

**A ligand-based scorer that cannot recover a known binder when fourteen of its
peers are in the model set will not recover a novel CReM product.** That is the
finding, and it is measured on molecules whose answer was already written down
— the route `how_this_project_breaks.md` scores as the best-performing of all.

Recorded as a recovery measurement, NOT a second gate: no verdict is graded and
no token is written, so there is only one verdict about this scorer for a reader
to act on.

## What this does not close

The scorer is measured and it does not rank. That leaves #4's Phase 2.2 half
done: the 3D hypothesis route is untouched, and it is the version of "use a
pharmacophore" that the literature actually means. It is now cheaper to
justify, because the gate it would have to clear runs in minutes (D0116) and
the bar it has to beat is written down.
