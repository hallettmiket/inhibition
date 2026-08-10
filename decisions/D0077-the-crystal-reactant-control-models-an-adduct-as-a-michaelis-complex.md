---
id: D0077
title: The crystal-reactant control models a covalent adduct as a Michaelis complex, and unrestrained MD destroys the construction before production starts
date: 2026-08-10
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/crystal_reactant.py
  - scripts/attack_sweep.py
  - decisions/D0075-the-sweep-rejects-every-known-active.md
  - decisions/D0076-the-dwell-filter-discards-the-approaches-the-observable-was-chosen-to-count.md
evidence:
  - 'rx_6VAJ (crystal reactant): reactive CH2 to SG = 1.98 A, Cl to SG = 3.77 A, S...C-Cl = 180.0 deg'
  - 'the 2.8-4.2 A near-attack window puts 1.98 A BELOW the floor — a formed bond, not an approach'
  - 'the 180.0 deg is constructed, not measured: crystal_reactant.py places the halogen along the S->C vector'
  - 'our docked mode-0 pose: CH2 to SG = 3.36 A (inside the window), S...C-Cl = 156.8 deg (clears the 150 deg bar)'
  - 'docked vs crystal: 5.01 A symmetry-corrected RMSD at only 1.45 A centroid separation — same cavity, pivoted'
  - 'rx_6VAJ enters PRODUCTION at 3.57 A / 100.9 deg, not 1.98 A / 180 deg: equilibration relaxes the construction away'
  - 'rx_6VAJ then spends 87.4% of the run inside the distance window at a median 78.6 deg — right place, wrong direction'
  - 'both controls start production with start_attack_ready = False'
---

# The control was the wrong shape

@tt8804, looking at the viewer: *"sulfopin is in the wrong pose … did our pipeline
select this pose mode over the real one?"* and then, immediately: *"I guess the
real pose for sulf is a covalent pose."*

That second sentence is the finding.

## Our pipeline did not pick the wrong pose

| | reactive CH₂→SG | Cl→SG | S···C–Cl |
|---|---:|---:|---:|
| crystal "reactant" (`rx_6VAJ`) | **1.98 Å** | 3.77 Å | 180.0° |
| our docked mode-0 | **3.36 Å** | 5.04 Å | **156.8°** |

The docked pose is a textbook near-attack conformation: inside the 2.8–4.2 Å
window and over the 150° SN2 bar, as docked. The screen produced a chemically
sensible answer.

The 5.01 Å RMSD between them, at only 1.45 Å centroid separation, is not a flipped
pose. It is the same cavity with the ligand pivoted — which is what bond formation
does.

## The control models the wrong state

6VAJ's deposited ligand is the **covalent adduct**. `crystal_reactant.py` cleaves
the bond and rebuilds the leaving group, and the result is treated as a
pre-reaction complex. It is not one:

* **1.98 Å is a formed bond.** It sits *below* the near-attack window's 2.8 Å
  floor, so the control cannot be attack-ready in its own starting frame — not
  because the geometry is bad but because the criterion is not defined there.
* **The 180° is constructed, not observed.** The halogen is *placed* along the
  S→C vector. Reporting it as if it were a measured approach angle reads the
  builder's assumption back as data.

## And the construction does not survive contact

`rx_6VAJ` does not enter production at 1.98 Å / 180°. It enters at **3.57 Å /
100.9°**. Equilibration relaxes the whole arrangement before the first production
frame, because once the bond is gone nothing restrains it — the 1.98 Å contact is
a strained non-bonded pair and the force field pushes it apart immediately.

It then spends **87.4%** of the run inside the distance window at a median
**78.6°**: in the right place, facing the wrong way, for the entire trajectory.

## What this does to D0075 and D0076

D0075 said the sweep rejects every known active. D0076 found the dwell filter
discarding the brief approaches the observable exists to count, and showed that on
raw visits within-mechanism the controls top their own chemistry.

This adds a third qualification, and it is the most serious: **for the `rx_*`
controls, the input geometry was never a valid Michaelis complex.** Their failure
is partly an artefact of how the control was built, so they cannot carry the
weight D0075 put on them. The sentence "the screen would not have surfaced the
incumbent" is not supported by `rx_6VAJ`.

The **docked** Sulfopin (`ref_Sulfopin__chloroacetamide`) does not share the
defect: it enters at a valid near-attack geometry, and it is the comparable
positive control.

## What follows

1. **Run the docked Sulfopin at 100 ns** as the real positive control. It is the
   one that enters in a state the criterion is defined on.
2. **Do not quote `rx_*` sweep numbers as evidence about the screen.** They are
   evidence about what unrestrained MD does to a cleaved adduct.
3. **If the adduct pose is to be modelled at all, restrain it.** A covalent
   complex held by a bond in the crystal needs the bond, or an equivalent
   restraint, to stay in that geometry. Unrestrained MD from a cleaved adduct
   measures relaxation, not binding.
4. `crystal_reactant.py` should say on its output that the pose it emits is a
   **construction**, and that its 180° is an input rather than a measurement.
