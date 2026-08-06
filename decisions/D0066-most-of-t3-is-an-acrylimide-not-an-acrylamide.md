---
id: D0066
title: 97% of T_3 carries a second electron-withdrawing group on the acrylamide nitrogen
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
  - 'CORRECTED COUNT, all activating groups on the acrylamide nitrogen: N-acyl 76.6%, N-urea 9.6%, N-sulfonyl 8.9%, N-carbamate 1.8% -- 96.8% of 2,403 scored T_3 candidates carry at least one'
  - 'genuinely plain N-alkyl / N-aryl acrylamides: 76 of 2,403 (3.2%)'
  - 'CORRECTED PREFERENCE: plain acrylamides score HIGHER than activated ones (median 1.78x vs 1.59x, AUC 0.398). The earlier 0.580 compared N-acyl imides against a pool dominated by the weaker sulfonyl/urea classes, not against plain amides.'
  - '51 of the 76 plain acrylamides reach the crystallographic Michael range (>=1.59x); 8 exceed the best crystal positive (2.69x)'
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


---

# Correction, 2026-08-05 (same day)

Two things above were wrong, and the corrected version is both worse and better.

## The scope was understated: 97%, not 77%

The original SMARTS, `[NX3](C(=O)C=C)C(=O)[#6]`, requires the second carbonyl to
be bonded to **carbon**. That misses three groups that are just as activating,
for the same reason — they compete for the same nitrogen lone pair:

| group on the acrylamide N | share of T_3 |
|---|---|
| N-acyl (a true imide) | 76.6% |
| N-urea / carbamoyl | 9.6% |
| N-sulfonyl | 8.9% |
| N-carbamate | 1.8% |
| **any of them** | **96.8%** |

**Only 76 of 2,403 scored T_3 candidates (3.2%) are genuinely plain N-alkyl or
N-aryl acrylamides** — the chemotype sulfopin is, and the chemotype both
crystallographic Michael acceptors are. The earlier "non-imide" count of 562 was
too generous by 486, because the sulfonyl, urea and carbamate compounds were
being counted as clean.

## The ranking does NOT prefer the activated compounds

Redone against genuinely plain acrylamides rather than against a mixed pool:

| | median enrichment |
|---|---|
| any activating group on N (n = 2,327) | 1.59× |
| **genuinely plain (n = 76)** | **1.78×** |

AUC 0.398 — the plain compounds score **higher**. The earlier figure (imides
beating non-imides, AUC 0.580) was a real comparison but a misleading one: the
"non-imide" side was dominated by sulfonyl and urea compounds, which score lower
than both. Comparing the activated chemotype against the one we actually want
reverses the direction.

**So the concern that the ranking promotes chemically undesirable compounds is
withdrawn.** It does not. The desirable chemotype scores best.

## What the real problem is

Not the ranking — the **generation**. T_3 produced 97% activated acrylamides
because its scaffold puts the attachment point on the nitrogen and LibInvent
almost always hung an electron-withdrawing group there. The approach explored
overwhelmingly the wrong chemotype, and the ranking is finding the right one in
spite of that.

**There is still a usable shortlist.** Of the 76 plain acrylamides, **51 reach
the crystallographic Michael range (≥ 1.59×) and 8 exceed the best crystal
positive (2.69×)**:

    2.94x  C=CC(=O)N(Cc1ccc(F)cc1)C1CCS(=O)(=O)C1
    2.82x  C=CC(=O)N(CC1CCCCO1)C1CCS(=O)(=O)C1
    2.82x  C=CC(=O)N(c1ccc(C)cc1)C1CCS(=O)(=O)C1
    2.82x  C=CC(=O)N(Cc1ccc(C)o1)C1CCS(=O)(=O)C1
    2.69x  C=CC(=O)N(Cc1ccc2c(c1)OCO2)C1CCS(=O)(=O)C1

## T_4 does not have this problem at all

The same test on T_4, which builds combinatorially from a curated R-group library
rather than generating decorations:

| | N-activated share |
|---|---|
| T_3 acrylamide (LibInvent, attachment on N) | **96.8%** |
| T_4 acrylamide (combinatorial, curated R-groups) | **0.0%** (0 of 187) |
| T_4 chloroacetamide | **0.0%** (0 of 187) |

Every T_4 acrylamide is a plain N-alkyl or N-aryl amide —
`C=CC(=O)N(c1ccccc1)[C@@H]1CCS(=O)(=O)C1` and the like. The constrained design
protected the chemotype; the generative freedom on the nitrogen destroyed it.

**This is the clearest head-to-head the project has between the two approaches**,
and it goes against the more sophisticated one. T_4 is 41% the size of T_3 and
essentially all of it is the chemotype we want, against 3% for T_3.

The consequence for shortlisting is therefore stronger than before, not weaker:
**a T_3 shortlist should be drawn from the plain acrylamides**, which is 3% of
the approach's output. Whether the activated ones are usable at all remains the
chemistry judgement for #12 §B4.
