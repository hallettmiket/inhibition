---
id: D0090
title: The pose cloud does not saturate, and the reason is that the whole cloud fits inside the scoring function's error bar
date: 2026-08-17
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - exp/5_mode_saturation/run_all.py
  - shared/pose_cluster.py
  - config/target.yaml
  - docs/build_plan_next.md
  - decisions/D0088-modes-come-from-pose-similarity-alone-clustered-with-hdbscan.md
evidence:
  - '@tt8804: "I dont understand how so many different poses can be such low energy??"'
  - 'ENERGY SPAN WITHIN ONE MOLECULE''S CLOUD: median 3.96 kcal/mol best-to-worst (IQR 3.56-4.33); best -7.34, median pose -5.59'
  - '3% of a molecule''s poses lie within 0.5 kcal/mol of its best, 15% within 1.0, 63% within 2.0, 95% within 3.0'
  - 'empirical docking scores carry ~2-3 kcal/mol error; blinded alchemical FEP is 2.44 kcal/mol RMSE and is rejected in this project as too imprecise (state_of_the_project §4)'
  - 'spearman(energy, viable) = -0.093 over 236,313 poses -- mean energy of viable poses -5.76 against -5.51 for the rest, a 0.25 kcal/mol difference'
  - 'within t4_716800c125a7, ten HDBSCAN modes have median energies from -4.58 to -6.32, a 1.7 kcal/mol span across ten DISTINCT places'
  - 'SATURATION, 6,000 poses docked and 5,390 PoseBusters-valid, ladder to 5,390: HDBSCAN modes fit y = a*n^b at b = 0.977 (R2 0.999); modes+singletons b = 1.025; covering number b = 0.879 at 1.0 A, 0.753 at 1.5 A, 0.628 at 2.0 A'
  - 'linear fit R2 = 1.000 against logarithmic 0.816 for mode count -- there is no plateau and no asymptote at any tested resolution'
  - 'at production depth 500 poses need 275 representatives at 2.0 A; de-duplication is 45%, and at 1.0 A it is 3%'
  - 'HDBSCAN noise fraction RISES with depth, 33% -> 40%, where filling in known regions would make it fall'
  - 'docking cost measured at 3 x 2,000 runs in 32 seconds, so depth is not budget-limited'
  - 'THE REACTIVE RECEPTOR IS NOT THE CAUSE (exp/12, 7 molecules, 2,000 runs each, both arms): the PLAIN receptor is FLATTER. Mean reactive-minus-plain energy span +1.89 kcal/mol and within-2-kcal fraction -0.307, consistent in 7 of 7'
  - 'plain-receptor clouds put 25-75% of poses within 2 kcal/mol of the best (mean 53%) against 3.5-42% for the reactive arm'
  - 'saturation is unchanged by the receptor: b at 3.5 A is 0.28-0.45 in both arms, mean reactive-minus-plain +0.072'
runbook: null
---

# D0090 — the cloud does not saturate, and the energies say why

## The question

@tt8804, on the saturation result: *"I don't understand how so many different
poses can be such low energy??"*

They are not. **They are indistinguishable, which is a different thing, and the
distinction is the finding.**

## What was measured

**The whole cloud fits inside the scoring function's error bar.**

| | |
|---|---:|
| energy span within one molecule's cloud, best to worst | **3.96 kcal/mol** (median) |
| poses within 0.5 kcal/mol of that molecule's best | 3% |
| poses within 1.0 | 15% |
| poses within 2.0 | **63%** |
| poses within 3.0 | **95%** |

An empirical docking score carries roughly **2–3 kcal/mol** of error. This
project already treats that magnitude as disqualifying: `state_of_the_project`
§4 rejects alchemical FEP partly because a blinded benchmark put it at
**2.44 kcal/mol RMSE**. AutoDock's function is not better than FEP.

So 63% of a molecule's poses sit within the tool's own uncertainty of the best
one, and 95% within a generous reading of it. **The energies cannot order these
poses.** "So many different poses at such low energy" is really "the landscape
this function reports is flat to within its own resolution."

Two corollaries, both measured:

* **Energy does not predict attack geometry.** ρ(energy, viable) = **−0.093**
  over 236,313 poses. Viable poses average −5.76 against −5.51 for the rest — a
  **0.25 kcal/mol** separation, an order of magnitude inside the error.
* **Distinct places are not energetically distinct.** Ten HDBSCAN modes of
  `t4_716800c125a7` have median energies spanning **1.7 kcal/mol**, from −4.58
  to −6.32. Ten different placements, one energy to within the noise.

## Why this explains the non-saturation

A search saturates when the landscape has basins: runs fall into them, the same
basins keep being found, and new sampling returns known answers. **A flat
landscape has no basins**, so every run terminates somewhere slightly new and
the count of distinct places grows with the number of looks.

That is exactly what was measured (`exp/5`, 5,390 PoseBusters-valid poses):

| metric | power-law exponent *b* in *a·n^b* | R² |
|---|---:|---:|
| HDBSCAN modes | **0.977** | 0.999 |
| modes + singletons | **1.025** | 1.000 |
| covering number @ 1.0 Å | 0.879 | 0.998 |
| covering number @ 1.5 Å | 0.753 | 0.995 |
| covering number @ 2.0 Å | **0.628** | 0.992 |

Mode count is **linear** — the linear fit is R² 1.000 against 0.816 for
logarithmic. The covering number is sublinear but a power law with a positive
exponent **has no asymptote**: doubling the poses at 2 Å still returns ~1.5×
the distinct places.

And the noise fraction **rises** with depth, 33% → 40%. If sampling were filling
in known regions it would fall.

## What this costs the plan

* **There is no 95% coverage depth to find.** §2.4e proposed deriving
  `docking.n_runs` from the depth at which coverage is 95% complete. Log growth
  would have permitted that; n^0.63 does not. The honest replacement is a
  diminishing-returns statement with the resolution named as a choice: *"at 2 Å
  and 500 runs we hold 275 distinct placements, and doubling the depth adds
  ~50% more."*
* **The search is not bounded and not local**, at least to 11× production depth.
* **De-duplication is weaker than assumed.** At production depth, 500 poses need
  **275** representatives at 2 Å — 45% collapse — and **97%** of poses are their
  own representative at 1 Å. Clustering is not mostly removing repeats.
* **`min_cluster_size` cannot fix it.** D0088 chose HDBSCAN because it never
  produces a bag; it never produces one because it has **no length scale** —
  only an absolute count of 3. Confirmed empirically here: it subdivides
  linearly forever. Any bounded count needs an explicit resolution, which is a
  decision we make rather than one the data yields.

## The reactive receptor was checked, and it is NOT the cause — `exp/12`

The obvious objection was that we soften the van der Waals parameters
(`R_EQ_12 = 3.2`, `EPS_12 = 1.0`) so the warhead can approach Cys113, which
flattens the landscape near the warhead by construction. Seven molecules were
docked at 2,000 runs into **both** the reactive receptor and the plain 3IKD.

**The plain receptor is flatter, in 7 of 7.**

| | reactive − plain (mean) | reading |
|---|---:|---|
| energy span across the cloud | **+1.89 kcal/mol** | reactive spans MORE |
| fraction within 2 kcal/mol of best | **−0.307** | reactive discriminates MORE |
| saturation exponent b at 3.5 Å | +0.072 | essentially unchanged |

On the plain receptor, **25–75% of poses (mean 53%) sit within 2 kcal/mol of the
best**, against 3.5–42% for the reactive arm. Stock AutoDock on this pocket is
*less* discriminating than our modified setup, not more.

So the flatness is a property of the scoring function on this pocket, and the
reactive parameterisation mildly improves it. **The finding stands and is
strengthened.** Saturation is untouched: b sits at 0.28–0.45 in both arms.

Stated as the experiment states it: this compares the reactive SETUP against the
plain one, and the reactive arm also carries reactive atom typing and a flexible
sidechain. It supports *"our reactive setup is not what flattens the
landscape"*, not a claim about the softened vdW term in isolation.

## What must be checked before this hardens

**One molecule for the saturation curve.** Everything above is `t4_716800c125a7`, which has four
rotatable bonds. A more rigid ligand should saturate faster. The deep cloud is
persisted and docking is 9 seconds per 2,000 runs, so two or three more
molecules is cheap and should precede any general claim.

## Consequence for the release

This does not block the PoseBusters quota (D0089), which stands on validity and
equal denominators rather than on saturation. It does remove the argument that
was going to justify the quota's *depth*, and it weakens the case that
clustering is mainly de-duplication — at 1 Å resolution there is almost nothing
to de-duplicate.
