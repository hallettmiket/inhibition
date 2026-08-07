---
id: D0074
title: Reactive docking protonates at pH 7.4, like the non-covalent path — so every pose of a titratable candidate so far is the wrong species
date: 2026-08-07
status: accepted
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/nac_screen.py
  - scripts/nac_screen_v2.py
  - shared/ionisation.py
  - scripts/md_residence_3ikd.py
  - docs/outline_2.2.0.md
evidence:
  - 'the two docking paths prepare ligands differently: noncovalent_dock_run uses obabel -p 7.4; nac_screen.prepare_ligand uses largest_fragment -> AddHs -> embed -> MMFF -> meeko, with NO pH protonation'
  - 'charge_ph74 is computed with obabel -p 7.4 (shared/ionisation.py) and therefore describes the NON-COVALENT species only'
  - 'charge_ph74 disagrees with charge(canonical_smiles) in 594 of 1782 T_4 rows (33.3%) and 331 of 5370 T_3 rows (6.2%)'
  - 'T_4 disagreement is strictly one-directional -- charge_ph74 >= SMILES charge, +1 in 534 rows and +2 in 60 -- consistent with obabel adding protons the reactive preparation never added'
  - 'T_3 disagreement is bidirectional (+1: 182, -1: 142, +2: 5, -2: 2) and is NOT assumed to share this cause'
  - 'md_residence reads charge_ph74 regardless of which path produced the pose; its guard refused three of five 100 ns launches rather than parameterise a mismatched species'
  - 'five T_2 seeds all have formal_charge 0 but four distinct charges at pH 7.4: ATRA -1, Du-Xu -1, Guo-Pfizer -2, Potter-Astex +1, sulfopin 0'
runbook: null
---

# Reactive docking protonates at pH 7.4

## The decision

**@tt8804: yes, protonate at 7.4.**

The reactive/covalent preparation adopts `obabel -p 7.4`, matching the
non-covalent path. Warhead-bearing candidates are docked as the species that
exists at physiological pH, not as the SMILES happens to be drawn.

## Why the question arose

Three of five 100 ns launches refused to start, because `md_residence` asked
antechamber for `charge_ph74` (+2) while holding a pose whose formal charges sum
to +1. Investigation showed the guard was right and neither column was stale:
**the two docking paths protonate differently, and `charge_ph74` describes only
one of them.**

`shared/ionisation.py` is explicit that this divergence was deliberate for the
non-covalent arm — `formal_charge` on the neutral SMILES is 0 for essentially
every molecule, so `charge_ph74` was created to record what was actually docked.
The reactive path was simply never brought into that scheme.

## What this costs, stated plainly

**Every reactive pose of a titratable candidate is of the wrong species.** That is
**594 of 1,782 T_4 molecules (33.3%)** and **331 of 5,370 T_3 (6.2%)** on the
newest frames. For those molecules, docking, the geometric criterion, consensus,
the anchoring score and the ranking all describe a molecule that does not exist at
pH 7.4.

A neutral amine and its ammonium are not a small perturbation in a pocket with a
basic cluster (Lys63/Arg68/Arg69) and an anionic-substrate-binding role. This is
not a rounding error to carry forward.

**Consequences that follow, and are not optional:**

1. The affected candidates must be **re-prepared and re-docked** protonated
   before their scores mean anything.
2. **`charge_ph74` becomes the single charge annotation** once both paths use the
   same preparation — the reason for the divergence disappears with it.
3. The two 100 ns runs launched today with `--net-charge +1`
   (`t4_da2e98512d02`, `t4_9a973be6b946`) are simulating the **unprotonated**
   species. They are internally consistent — pose and charge agree — but they are
   the species this decision rules against, and they must be redone from
   protonated poses. Recorded rather than quietly rerun, because the trajectories
   exist and could otherwise be quoted.

## What is NOT decided here

**The T_3 bidirectional disagreement is a separate defect.** T_4's mismatch is
strictly one-directional and fully explained by the missing protonation step.
T_3's runs both ways, including 142 rows where `charge_ph74` is *lower* than the
SMILES charge, which protonation cannot cause. That has its own root cause and
must not be assumed fixed by this change.

## The guard stays

`md_residence`'s refusal to parameterise when the requested charge and the pose's
formal charges disagree is the only thing in the pipeline that noticed this. A
wrong `-nc` does not fail — it produces a force field for the wrong species and a
100 ns trajectory that looks perfectly healthy. **Do not soften it to a warning.**
