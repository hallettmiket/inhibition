---
id: D0029
title: Three of T_4's nine warhead classes are one class after the reaction
date: 2026-07-28
status: accepted
approach: t4
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - approaches/t4_combinatorial/04_rank_within_class.py
  - approaches/t4_combinatorial/06_mmgbsa.py
  - approaches/t4_combinatorial/05_regiochemistry_comparison.py
  - config/approaches/t4_combinatorial.yaml
  - integration/app/DECISIONS_TAB_SPEC.md
evidence:
  - '1,683 docked rows collapse to 1,309 unique adducts (dock_id)'
  - '100% of chloroacetamide, sulfamate_acetamide and sulfonate_acetamide share their adduct with the other two: 187 adducts, each reached by all three routes'
  - 'the other six classes have 0% cross-class adduct sharing'
  - 'affinity_kcal is identical to the last decimal across all 187 triples: max |chloro - sulfamate| = 0.0, max |chloro - sulfonate| = 0.0'
  - 'the shortlist quota gave these three classes 9 slots for 3 unique molecules; the 27-candidate shortlist is 21 molecules'
  - '374 of 1,683 covalent docks (22%) recomputed a pose already computed'
  - '05_regiochemistry_comparison compared only bdhi_c4/c5 and naphthoquinone_c2/benzo, so no reported comparison is affected'
---

# Three of T_4's nine warhead classes are one class after the reaction

## What was found

T_4 enumerates nine warhead classes. Three of them —
`chloroacetamide`, `sulfamate_acetamide` and `sulfonate_acetamide` —
differ only in their leaving group:

| class | warhead | leaving group | adduct |
|---|---|---|---|
| `chloroacetamide` | `C(=O)CH2-Cl` | `Cl-` | `C(=O)CH2-S-Cys113` |
| `sulfamate_acetamide` | `C(=O)CH2-OSO2NH2` | sulfamate | `C(=O)CH2-S-Cys113` |
| `sulfonate_acetamide` | `C(=O)CH2-OSO2R` | sulfonate | `C(=O)CH2-S-Cys113` |

The leaving group is, by definition, *gone* once the bond to Cys113 has
formed. So all three deliver the **same bound species**. Since D0022
made docking operate on the adduct form, all three now dock the same
molecule, and `dock_id` — which hashes the adduct SMILES — correctly
assigns them one identity.

This is not a defect, and it was not news. `03_covalent_dock.py` has
documented it since the D0022 re-dock — "all three give an IDENTICAL
adduct — verified for all 198 R-groups" — and deliberately docks them
once. The docking layer had it right.

**The defect is that the fact stopped there.** Every stage after
docking still counts to nine.

## How it surfaced

Two shortlisted candidates in different warhead classes produced
byte-identical MM-GBSA inputs and an identical minimised complex energy
(-12304.39 kcal/mol). That looked like a pose-lookup bug — the working
directory is keyed on `candidate_id`, the pose is fetched by `dock_id`,
and those are deliberately not the same thing. It was not a bug: the
two candidates genuinely are one molecule, exactly as the docking stage
said they would be.

The only real difference between the two `complex.prmtop` files was
tleap's `DATE =` stamp.

That a correct, documented, load-bearing fact was recorded in one
module's docstring and honoured nowhere else is the transferable
lesson here. A chemistry fact that changes how results must be counted
belongs in the library or the config, where every stage reads it — not
in the prose of the one stage that happened to discover it.

## Consequences

**The class quota is triple-counting one chemotype.** Every class gets
three shortlist slots. The acetamide family therefore received nine
slots for three molecules, and the shortlist that MM-GBSA is scoring is
27 rows over 21 distinct molecules. Six of the 27 MM-GBSA runs
recompute an identical system.

**Any "diversity across warhead classes" claim overstates by two.**
T_4 covers **seven** post-reaction chemotypes, not nine. The
integration GUI must say seven.

**D0020 is not violated, it is confirmed.** D0020 forbids comparing
affinity across warhead classes. Here three classes are the same class,
so their scores are comparable — and are exactly equal. The one case
where cross-class comparison is legitimate is the case where there is
no cross-class comparison to make.

## What stays distinct

The three classes remain genuinely different **before** the reaction,
and that difference is not cosmetic:

- leaving-group ability sets the labelling rate, so the reactivity
  triage (step 2, D0019) legitimately separates them;
- they differ in off-target liability and in hydrolytic stability;
- they are three different syntheses.

The right model is **one adduct class reached by three warhead routes**.
That is more useful than nine classes, not less: for each shortlisted
adduct we can now say there are three independent synthetic routes to
the same bound molecule, and the choice among them is a kinetics and
selectivity decision rather than a binding one.

## Decision

1. Post-reaction ranking, quotas and diversity counts group by
   **adduct class** (7), not warhead class (9).
2. Warhead class is retained on every row and carried into the GUI as
   the *route*, because the reactivity triage needs it and because
   three routes to one adduct is a result worth showing.
3. The shortlist is re-derived on adduct class once the current MM-GBSA
   run completes; the six freed slots go to molecules not yet scored.
   The in-flight results stay valid — they are correct energies for the
   molecules they name, merely redundant.
4. Docking should key its cache on `dock_id` so the 22% redundant work
   is not repeated on a re-run.

## The recurring error

This is the third instance of one mistake: **identifying a chemotype by
a pre-reaction or reactive-atom feature when the object being scored is
post-reaction.**

- D0025: `alpha_halo_carbonyl` straddled the core/decoration boundary
  because the warhead was identified by a two-atom reactive SMARTS.
- D0028: decoy warhead classes were assigned by the reactive-atom
  SMARTS, so `[CH2][Cl]` labelled cyclophosphamide a chloroacetamide.
- D0029: nine classes were counted pre-reaction and scored
  post-reaction.

A reactive-atom SMARTS tells the docking engine where to form a bond.
It is not a name for a chemotype, and it is never a name for what
remains after the bond forms.
