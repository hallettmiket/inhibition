---
id: D0058
title: Charge is measured at pH 7.4 with the tool docking used, and phosphate is labelled not filtered
date: 2026-08-05
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/ionisation.py
  - scripts/annotate_ionisation.py
  - tests/test_ionisation.py
evidence:
  - 'descriptors.formal_charge is 0 for all 4,803 T_1 rows, all 1,882 T_2 rows and 5,378 of 5,396 T_3'
  - 'the five T_2 seeds are all formal_charge 0 and span four charge states at pH 7.4: atra -1, du_xu -1, guo_pfizer -2, potter_astex +1, sulfopin 0'
  - 'obabel silently converted 143 of 4,803 T_1 molecules, exit 0, empty stderr; the 144th is OC[P@TB14](O)(O)(O)CO'
  - 'split-retry recovered T_1 from 4,661 unknown to 79; only 41 molecules genuinely fail'
  - 'charge stratification changes the top-25 by 28-64%: t1 18/25 overlap, t2 16/25, t3 9/25, t4 16/25'
  - 'T_1 carries 394 phosphate-containing molecules; T_2/T_3/T_4 carry none'
---

# The charge column had to be built, because the one that existed was the wrong molecule

## Context

#6 items 5 and 7 -- "T_2 phosphate: label, not protect" and charge-stratified
ranking -- were both decided and neither implemented. They land together
because labelling phosphate is only honest once phosphate-free molecules are
actually evaluated rather than filtered out by stratification.

## Why the existing column could not be used

`descriptors.formal_charge` already exists, and it is **0 for essentially every
molecule in the project**. That is not a defect in the descriptor: it is
`Chem.GetFormalCharge()` on the neutral canonical SMILES, exactly what it says.

But docking protonated for pH 7.4, so the molecule that column describes is not
the molecule that was scored. Stratifying on it would have produced ONE stratum
and looked like it had worked -- a populated, plausible column that does not
mean what a reader wants it to mean. The five T_2 seeds make it concrete: all
five are `formal_charge = 0`, and at pH 7.4 they are **-1, -1, -2, +1, 0**.
Four states reported as one.

## Decision

Charge is computed with **`obabel -p 7.4`, the same call
`noncovalent_dock_run` used**, so the stratum and the score describe the same
structure by construction rather than by a pKa model that would need its own
validation. `has_phosphate` is a SMARTS label, never a filter: Pin1 BINDS
phosphate, a rule against phosphorus was proposed and discarded for rejecting
four known binders, and the permeability cost is a chemist's judgement (#12).

`charge_class` is deliberately coarse -- anion / neutral / cation / unknown.
**Vina carries no electrostatic term**, so a stratum is not a claim about
modelled energy; it is the narrower statement that comparing a dianion's score
with a cation's compares two different physical situations. Splitting -1 from
-2 would imply a resolution the argument does not have. `unknown` is its own
class rather than folded into `neutral`, so conversion failures do not land in
the largest stratum.

Stratified ranking needed no new machinery: `rank(group_col="charge_class")`
is the same mechanism T_4 already uses for warhead class.

## Measured effect

Charge stratification moves a large share of every shortlist:

| arm | top-25 overlap with unstratified | charge classes |
|---|---|---|
| T_1 | 18/25 (72%) | 1,574 neutral · 831 cation · 770 anion |
| T_2 | 16/25 (64%) | 1,302 anion · 513 neutral · 67 cation |
| **T_3** | **9/25 (36%)** | 3,830 neutral · 125 cation · 108 anion |
| T_4 | 16/25 (64%) | 1,012 neutral · 671 cation |

T_3 moves most because it is 94% neutral: stratifying lets its small charged
populations surface instead of being buried. **The ranking was partly ordering
by charge state**, which is the artefact this was adopted to remove -- and it
is the same shape as the size sort D0049 removed.

This does NOT make the ranking valid. `rank_validated` stays False.

## A defect found on the way, worth its own line

**obabel silently truncates a batch on one bad molecule.** It reported "143
molecules converted" from a 4,803-line T_1 file, wrote nothing to stderr and
exited 0. The 144th is `OC[P@TB14](O)(O)(O)CO` -- pentavalent phosphorus with
trigonal-bipyramidal stereochemistry, a DiffSBDD artefact -- and everything
after it was lost. 97% of the arm, with no error of any kind.

Docking never hit this because `_prepare_one` converts one molecule per file.
Batch conversion is not immune, so results are matched back **by the id obabel
echoes** and short returns trigger a recursive split until failures are
isolated. Identity matching is what made it detectable at all: by position, a
short return is invisible and every later charge lands on the wrong candidate.
Recovery: 4,661 unknown to 79, of which only 41 genuinely fail.
