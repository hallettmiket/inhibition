---
id: D0048
title: Two synthesizability rules adopted — acyl phosphate and stereogenic phosphorus, both T_1-only
date: 2026-07-31
status: accepted
approach: shared
decided_by: '@mhallet'
origin: adversary
supersedes: []
superseded_by: null
affects:
  - shared/synthesizability.py
  - shared/alerts.py
  - scripts/reshortlist_synthesizable.py
evidence:
  - 'both rules checked against every parseable reference binder: 0 rejected'
  - 'rejected candidate rule `any phosphorus [P]`: kills 4 binders (Wildemann, Guo-Pfizer, Liu-Pei, Jiang-Pei)'
  - 'rejected candidate rule `alkyl phosphate [CX4][OX2][PX4]`: kills the same 4'
  - 'rejected candidate rule `catechol`: kills EGCG'
  - 'scope: 236/4803 T_1 candidates fire; 0/1882 T_2, 0/5396 T_3, 0/1782 T_4'
  - 't1_db179d172dda (was rank #10) fires both and left the shortlist'
  - 'SAscore rated that molecule 3.88 — the LOWEST (easiest) of T_1 top 10'
  - 'its SAscore fragment term is +1.201: every ECFP4 environment is individually common'
  - 'AiZynthFinder: UNSOLVED after 3000 iterations / 7968 routes; 183 of 197 tree molecules still carry the acyl phosphate'
  - 'nearest purchasable analogue: max Tanimoto 0.51 across 4.5B Enamine REAL compounds'
---

# Two rules that survive the reference set, and three that did not

## What was adopted

**`acyl_phosphate`** — `[CX3](=[OX1])[OX2][PX4]`. A mixed
carboxylic-phosphoric anhydride: a high-energy acyl-transfer group biology uses
*because* it is transient. In a 5-membered ring it is doubly activated, since
cyclic phosphates hydrolyse ~10^6–10^8 faster than acyclic ones from ring
strain and both hydrolysis modes relieve it.

**`stereogenic_phosphorus`** — a predicate, not SMARTS, because "four different
substituents on P" is not a substructure. P-stereogenic synthesis needs a chiral
auxiliary or a resolution even where most money has been spent on it.

`Rule` gained an optional `predicate` for the second, and `_n_stereogenic_
phosphorus` uses `FindPotentialStereo` rather than `FindMolChiralCenters`
because CIP labelling raises *"Digraph generation failed: more than 100000
nodes"* on the large peptidic reference binders — **a rule that crashes on a
known binder is worse than no rule.**

## What was rejected, recorded so nobody re-proposes it

| candidate | known binders killed |
|---|---|
| `any phosphorus` `[P]` | **4** — Wildemann, Guo-Pfizer, Liu-Pei, Jiang-Pei |
| `alkyl phosphate` `[CX4][OX2][PX4]` | the same 4 |
| `catechol` | 1 — EGCG |
| any 7-membered carbocycle | 0, but see below |

**Pin1 is a phosphate-binding enzyme.** A phosphorus ban bans the target's own
pharmacophore. `acyl_phosphate` is safe *only* because it demands a carbonyl
carbon on the ester oxygen, which no reference phosphate has. That single atom
is the whole difference between a usable rule and one that deletes four known
binders.

The 7-ring case is deliberately **not** a rule: 523 of 4,803 T_1 candidates
(10.9%) contain one, which is a real signal that DiffSBDD's 3D→graph bond
inference may be closing 7-rings where 6-rings belong — but benzodiazepines,
tropanes and colchicine all have 7-rings, and the molecule that prompted this
has a 6-5-7 cyclohepta[b]indole core that is a real, purchasable ring system.
**Instrument the rate as a generator-health metric; do not filter on it.**

## Why the rules were needed: SAscore cannot see a combination

`t1_db179d172dda` carried an acyl phosphate, a catechol and a stereogenic
phosphorus. **SAscore rated it 3.88 — the lowest, i.e. easiest, of T_1's entire
top 10** (range 3.88–5.89), and all seven existing rules passed it.

The reason is structural, not a tuning problem: its SAscore fragment term is
**+1.201**, because every ECFP4 environment in it (aryl-OMe, C-O-P, P=O,
phenol) is individually *common*. **It is the combination that has no
precedent, and a fragment-additive score cannot represent a combination.**
Median SAscore of T_1's 24 acyl-phosphate molecules is 4.38 against 4.14
overall — no discrimination at all.

Independent confirmation: AiZynthFinder left it **unsolved after 3,000
iterations and 7,968 routes**, with 183 of the 197 molecules in its retained
trees still carrying the acyl phosphate — the search never disconnected it.
Nearest purchasable analogue is Tanimoto 0.51 across 4.5 billion Enamine REAL
compounds.

## Scope: this is a DiffSBDD-specific failure

| approach | candidates firing |
|---|---|
| **T_1** | **236 / 4,803** |
| T_2 | 0 / 1,882 |
| T_3 | 0 / 5,396 |
| T_4 | 0 / 1,782 |

T_2, T_3 and T_4 are seeded from real molecules and generate by fragment swap,
decoration or enumeration, so neither motif can arise. T_1 is the zero-seed arm
and carries no synthetic information by construction.

## Related change: T_1's alerts were never gated

`alert_gate_pass` was `True` for all 4,803 T_1 rows — not because they passed,
but because `disqualifying` defaults to `()` and no caller ever passes it, so
the gate is a no-op by construction. The column read as "passed the alert
filter" and meant "no filter was configured."

**PI decision: report, do not gate.** T_1's alerts stay advisory — a catechol
rule would reject EGCG, and two earlier alert-derived rules died the same way.
But the column now says which it is: `alert_gate_pass` is `NA` when no gate ran,
and `alert_gate_applied` records the fact.

## Consequence

Shortlists rebuilt. T_1 lost 3 rule-failures including `t1_db179d172dda`, the
molecule that prompted the rules and which had reached rank #10 with MM-GBSA
already computed on it.
