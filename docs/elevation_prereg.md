# Pre-registration — which ranking metric selects for physical stability?

*Written and committed **before** any elevation simulation was run. @tt8804,
2026-08-06. The point of writing it first is that with four groups and several
possible readouts, a result chosen after the fact can be made to say almost
anything — and D0045 exists in this project because that has happened.*

---

## The question

Two metrics rank the same molecules and disagree almost completely:

- **enrichment** — the fraction of docking runs reaching a near-attack
  conformation. Does **not** converge (D0068): the same molecules fall from
  2.91× to 0.96× at 10× the search effort, and its rank correlation across
  efforts is **ρ = −0.117**.
- **consensus** — whether a molecule's **top-10 poses by energy** agree.
  Rank-stable across efforts (**ρ = +0.568**, D0070).

A sanity check @tt8804 asked for exposed how far apart they are: **397 molecules
have a single binding mode with ≥0.90 pose agreement, and only 4 of them clear
enrichment > 5.70.** The enrichment cut captures **1%** of the well-aligned
molecules and misses five chloroacetamides — the one warhead class with a
validated criterion (AUC 0.908).

So: which of the two is selecting for something physical?

## The design, and the three confounds that shaped it

Each confound was measured, not assumed.

**1. Consensus is partly rigidity** (ρ = −0.259 vs rotatable bonds, p = 4×10⁻⁸⁹).
Rigid molecules have fewer ways to sit *and* are trivially more stable under
dynamics. Groups are therefore **matched on rotatable-bond count** — median 4.0
across all three BDHI groups — not merely balanced in size.

**2. Enrichment and warhead class are confounded.** High-enrichment cells are
BDHI-dominated; low-enrichment cells are acrylamide-dominated. The two cannot be
separated across the whole set, **only within BDHI**. The core comparison is
therefore BDHI-only, which is a limit on what any result can claim.

**3. High enrichment is nearly a subset of high consensus.** The
high-enrichment/low-consensus cell holds 5 molecules. Enrichment is not an
independent axis but a stricter filter inside consensus, so the informative
contrast is A vs B rather than a four-cell factorial the data cannot fill.

### Groups

| group | n | enrichment | consensus | rotb | what it isolates |
|---|---|---|---|---|---|
| **A** hi-enr / hi-cons, BDHI | 8 | 5.63 | 1.000 | 4.0 | — |
| **B** lo-enr / hi-cons, BDHI | 8 | 1.90 | 1.000 | 4.0 | **vs A: enrichment alone** |
| **D** lo-enr / lo-cons, BDHI | 8 | 1.68 | 0.200 | 4.0 | **vs B: consensus alone** |
| **V** hi-cons, chloroacetamide | 5 | 1.42 | 1.000 | 5.0 | the validated class |
| **REF** crystallographic positives | ≤8 | — | — | — | **the anchor** |

**REF is what makes this interpretable.** Crystallographic positives are the only
molecules we know actually react with Cys113. Without them, "A beats D" compares
two arbitrary groups; with them the question becomes *which group behaves like a
molecule that really reacts*.

## Readouts, both fixed now

**Tier 1 — equilibration survival (free).** `gromacs_explicit` applies no
position restraints during NVT/NPT, so 300 ps of unrestrained dynamics already
runs before any bias. The warhead's displacement from its docked position across
that window is a physical stability signal we were discarding. Reported as
Δd(warhead→Cys113 SG) between the docked pose and the start of production.

**Tier 2 — BPMD.** The bias needed to push the warhead out of the near-attack
window (`shared/bpmd.py`). Note tier 2 starts from the **post-equilibration**
pose, so tier 1 also qualifies tier 2: a molecule that drifted 3 Å during
equilibration is not having its docked pose tested at all.

## Readings, fixed in advance

| observation | conclusion |
|---|---|
| **B ≈ A**, both > D | **Consensus is the filter; enrichment adds nothing.** The shortlist should be drawn on consensus, and the 397 are the real candidate pool. |
| **A > B**, both > D | Enrichment adds something real beyond consensus. Both belong in the ranking. |
| **A ≈ B ≈ D** | Neither metric predicts stability. BPMD is measuring something orthogonal, and the ranking has no physical support from this experiment. |
| **D ≥ A, B** | Something is wrong with the design or the metric direction; report as a failure, do not reinterpret. |
| **V ≈ REF** | Consensus-selection inside the validated class finds molecules that behave like known binders — the strongest available result, and the one that would justify a synthesis shortlist. |

**Statistics.** n = 8 per BDHI group supports only large effects. Group
comparisons by Mann-Whitney with the effect size reported alongside; **no
significance claim will be made from n = 5** (group V). Every value is reported
with its replica spread, because a stability score without its spread invites an
ordering the data cannot carry.

## What this cannot settle

- **It is BDHI-only** for the A/B/D contrasts. BDHI is the class that collapses
  at convergence (7.23× → 0.60×) and has **zero crystallographic positives**, so
  a positive result here does not transfer to other classes without re-testing.
- **BPMD tests stability, not reactivity.** A stable near-attack pose is
  necessary for the reaction, not sufficient. Nothing here measures whether a
  molecule reacts.
- **n = 5 for the validated class**, because only five chloroacetamides clear the
  consensus bar. That arm is descriptive.
