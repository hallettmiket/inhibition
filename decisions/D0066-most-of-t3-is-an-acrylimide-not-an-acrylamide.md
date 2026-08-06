---
id: D0066
title: 77% of T_3 is an acrylimide, not an acrylamide, and the ranking mildly prefers them
date: 2026-08-05
status: proposed
approach: T_3
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/nac_rank.py
  - docs/ranking_rationale.md
evidence:
  - 'T_3 scaffold as generated: [*]N(C(=O)C=C)C1CCS(=O)(=O)C1 -- all 4,086 kept candidates share it, and the attachment point is on the NITROGEN'
  - '3,138 of 4,086 kept T_3 candidates (76.8%) match [NX3](C(=O)C=C)C(=O)[#6]: the acrylamide nitrogen carries a second acyl group'
  - 'sulfopin (the parent) is NOT an imide -- its nitrogen carries one acyl and two alkyls'
  - 'neither crystallographic Michael positive is an imide: 9INN COC(=O)/C=C/C(=O)N[C@@H]1CCS(=O)(=O)C1, 9JF6 COC(=O)/C=C\\C(=O)N(Cc1cc2ccccc2s1)[C@H]3CCS(=O)(=O)C3'
  - 'ranking prefers them: imide median enrichment 1.59x vs 1.41x non-imide, AUC 0.580, p=3.3e-05 (n=937 vs 267)'
  - 'top 25 by enrichment is 88% imide, top 50 84%, against 78% of the pool'
---

# Most of T_3 is an acrylimide, and the ranking mildly prefers them

## What was found

Scrutinising the top-ranked T_3 candidates, every one had the same feature: the
acrylamide nitrogen carries a **second acyl group**.

    T_3 scaffold      [*]N(C(=O)C=C)C1CCS(=O)(=O)C1
    top hit           C=CC(=O)N(C(=O)C=C(NC(=O)C=Cc1ccc2c(c1)OCO2)c1ccccc1)C1CCS(=O)(=O)C1
                           ^^^^^^^^^^^^ two acyls on one nitrogen

It is not a property of the top of the list. **3,138 of 4,086 kept T_3
candidates — 76.8% — are N,N-diacyl.** The scaffold's attachment point sits on
the nitrogen, so LibInvent was asked to decorate exactly there, and three times
in four it attached an acyl group.

**Neither the parent nor the ground truth shares this.** Sulfopin's nitrogen
carries one acyl and two alkyls — a plain amide. Both crystallographic Michael
acceptors at Cys113 (9INN, 9JF6) are plain amides too.

## Why it matters

An N-acyl amide is an **imide**, and the second carbonyl competes for the
nitrogen lone pair. The nitrogen can no longer donate into the acrylamide
carbonyl, which leaves the β-carbon substantially more electrophilic than in a
normal acrylamide, and the imide itself hydrolytically labile in aqueous buffer.

The practical consequences are the classic covalent-inhibitor failure modes:
a warhead that reacts with everything it meets, and one that may not survive
dilution into an assay long enough to reach the target at all.

**This is a chemistry judgement and is flagged, not decided here.** The reasoning
above is a non-chemist's reading; it is exactly the kind of nameable structural
liability issue #12 §B4 asks about ("what else should we be filtering that we are
not?"), and it needs the Lu lab's answer before anything is filtered.

## The ranking mildly prefers them, and that is probably real

| | median enrichment | > 2× |
|---|---|---|
| imide (n = 937) | **1.59×** | 25% |
| non-imide (n = 267) | 1.41× | 20% |

AUC 0.580, p = 3.3 × 10⁻⁵. A real preference, and a modest one.

The likely cause is geometric rather than artefactual: a second acyl locks the
nitrogen planar and restricts rotation about the N–C(O) bond, pre-organising the
acrylamide instead of letting it sample freely. That is a genuine
warhead-presentation advantage and the criterion is right to see it.

**It is still a problem in practice**, because the criterion is deliberately
blind to reactivity. `docs/ranking_rationale.md` rests on step 2's rate being a
property of the warhead *class*, so that all molecule-to-molecule variation lives
in step 1 — and an imide is a **different class** from an amide, with a different
intrinsic rate. The assumption that licenses ranking on geometry alone is exactly
the assumption these compounds break.

## Consequence

**A shortlist taken straight off this ranking would be ~88% imide.** Any T_3
shortlist must therefore be reported with its imide fraction stated, and split
imide / non-imide, until the chemistry judgement comes back. Nothing is filtered
yet — five previously proposed rules were discarded for rejecting molecules
someone had actually made (#12 §B3), and that discipline holds here.

## What this does not settle

- Whether an acrylimide is acceptable to make and assay at all. Chemist's call.
- Whether the geometric preference survives the robustness run. Measured on 1,204
  of 4,086 scored so far; to be re-checked on the full set.
- Whether T_3 should be **re-generated** with the attachment point moved off the
  nitrogen, which is the upstream fix if imides are ruled out. That would
  invalidate the T_3 arm rather than filter it, so it is not proposed lightly.
