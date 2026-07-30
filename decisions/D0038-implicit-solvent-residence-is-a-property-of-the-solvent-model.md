---
id: D0038
title: The two solvent models disagree; but the dissociation I blamed on water was a single-trajectory artefact
date: 2026-07-30
status: partially_withdrawn
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/gromacs_explicit.py
  - shared/gromacs_analysis.py
  - scripts/run_gromacs_explicit.py
  - scripts/merge_gromacs_results.py
  - decisions/D0036-ensemble-mmgbsa-is-precise-and-still-below-chance.md
evidence:
  - '48 of 48 non-covalent candidates (T_1, T_2) run for 10 ns in explicit TIP3P, 0 failures'
  - 'WITHDRAWN: t1_8a3f4861ac34 9.00 nm was NOT reproducible; implicit re-run gave 1.75 nm (engaged 0.51)'
  - 'WITHDRAWN: t1_bd563e94c862 7.30 nm was NOT reproducible; implicit re-run gave 0.59 nm (engaged 0.91)'
  - 'run-to-run divergence under the SAME model: 5.1x and 12.5x on mean ligand RMSD for those two'
  - 'STANDS: spearman(implicit, explicit) = -0.102 run 1 and -0.144 run 2, i.e. two independent implicit runs both uncorrelated with explicit'
  - 'Spearman(implicit RMSD, explicit RMSD) = -0.102 across 47 paired candidates'
  - 'both models flag a similar COUNT as leaving (4 vs 5 of 47) but they are different candidates'
  - 'candidates stable under GB drift furthest in water: t2_bc8a4b62eb0e 1.78 -> 4.86 nm, t1_c1ec9e35dba7 0.63 -> 4.66 nm'
  - 'GROMACS 2026.3 CUDA sees all 8 A100s at ~740 ns/day; the shared OpenCL build refuses NVIDIA devices entirely'
---

# Two claims, one withdrawn

## Withdrawal notice

This record made two claims. One stands and is now better supported; the other
is withdrawn.

**STANDS — implicit and explicit residence are uncorrelated.** A second,
independent implicit-solvent run gives Spearman **-0.144** against explicit,
beside **-0.102** for the first. Two runs agree that they disagree with
explicit water, so this is not one unlucky trajectory.

**WITHDRAWN — the dissociation was not caused by the water model.** It did not
reproduce under the SAME model. t1_8a3f4861ac34 went 9.00 nm -> **1.75 nm** on
re-run (engaged 0.07 -> 0.51); t1_bd563e94c862 went 7.30 nm -> **0.59 nm**
(0.14 -> 0.91). Both now sit close to their explicit values, 1.52 and 0.47 nm.
The prose below attributes those departures to the absence of explicit water.
That attribution is unsupported: they were a property of a single short
trajectory.

Velocities are drawn afresh each run, so two 2 ns trajectories of one molecule
under one model diverge -- and for candidates that wander they diverge by
**5.1x and 12.5x** on mean ligand RMSD, while stable candidates reproduce to
within 0.7-2.1x. The metric is reproducible for molecules that sit still and
unreproducible for precisely the molecules it was being used to flag.

**The corrected lesson, which is stronger.** A per-candidate residence claim
needs replicate trajectories in either solvent model. One run cannot separate
"this ligand leaves the pocket" from "this trajectory wandered", and the
original reporting did not try to.


## What was claimed

D0036 reported pocket residence from 2 ns of GB implicit-solvent MD, and two
T_1 candidates were singled out as leaving the pocket outright: ligand RMSD
9.00 and 7.30 nm, engaged in 0.07 and 0.14 of frames. That was presented as the
ensemble tier earning its keep -- a behaviour a single minimised structure
could not have shown, because a minimisation cannot dissociate.

## What explicit water says

Both stay bound.

| candidate | implicit RMSD | implicit engaged | explicit RMSD | explicit engaged |
|---|---|---|---|---|
| t1_8a3f4861ac34 | 9.00 nm | 0.07 | **1.52 nm** | **0.977** |
| t1_bd563e94c862 | 7.30 nm | 0.14 | **0.47 nm** | **0.746** |

Ten nanoseconds of TIP3P, the same complex, the same GAFF2/ff19SB
parameterisation, the same starting pose. The dissociation does not reproduce.

The mechanism is not mysterious. GB implicit solvent contains no water
molecules at all; it approximates their averaged effect with a formula. There
is nothing in the model to occupy the space a departing ligand would have to
cross, and nothing to pay a desolvation cost for leaving. A shallow,
solvent-exposed pocket is precisely where that approximation is weakest.

## The stronger result: the two models do not agree at all

Across the 47 candidates run under both:

**Spearman(implicit RMSD, explicit RMSD) = -0.102.**

Not weakly correlated -- uncorrelated, and marginally negative. Each model
flags a similar NUMBER of candidates as leaving (4 of 47 implicit, 5 of 47
explicit), which would look like agreement in a summary table. They are
different candidates. Molecules that sat still under GB drift furthest in
water: t2_bc8a4b62eb0e 1.78 -> 4.86 nm, t1_c1ec9e35dba7 0.63 -> 4.66 nm.

So implicit-solvent residence carries no information about explicit-solvent
residence. It is a measurement of the solvent model.

## What this retracts

Every residence-based statement in D0036 and in the reporting around it is
withdrawn as a claim about molecules. Specifically:

- "t1_8a3f4861ac34 leaves the pocket entirely" -- an artefact.
- "a candidate docking scored well enough to shortlist dissociates" -- an
  artefact.
- The framing of residence as something the ensemble tier revealed that
  minimisation could not. It revealed a property of GB.

D0036's central finding is untouched: it concerns dG and enrichment, not
residence, and none of it rests on these numbers. The residence metric was
always reported as descriptive rather than as a ranking, and was never allowed
into the gate -- which is the only reason this correction is cheap.

## What survives

Residence remains a legitimate quantity in EXPLICIT solvent, where there is
water to be displaced. It still is not a discriminator: 42 of 47 candidates
stay engaged above 0.75, so the metric has almost no dynamic range even in
TIP3P, and for the covalent approaches it cannot vary at all because the
ligand is bonded to Cys113.

The one thing it is good for is catching a candidate whose pose is not
physical, and it should be read that way -- a QC flag, in the solvent model
that can support one.

## Practical notes

**The GPU needed a different GROMACS.** The shared conda-forge build in
`dwi_amber_md` is compiled with OpenCL, and modern GROMACS refuses NVIDIA
devices under OpenCL -- `mdrun -gpu_id 0` reports "incompatible devices",
which reads like a hardware fault and is not one. A CUDA build
(`dwi_gromacs_cuda`, GROMACS 2026.3) sees all 8 A100s and runs a 30k-atom
solvated system at ~740 ns/day against ~85 on 16 CPU cores. The OpenCL binary
is still first on PATH.

**Solvate with tleap, not `gmx solvate`.** GROMACS inserts SPC water
coordinates; an Amber topology expects TIP3P. Solvating on the Amber side and
converting the whole system with parmed keeps the water model consistent and
reuses the existing ligand parameterisation, so the explicit run describes the
same molecule everything else scored.

**Contacts are not comparable between the tiers.** The implicit metric counted
heavy-atom pairs within 0.45 nm; `gmx mindist -on` uses its own definition and
returns several-fold smaller numbers on the same complex (4.1 against ~25
measured). The columns are named differently so nothing can subtract one from
the other and report a solvent effect that is a definition change. Ligand RMSD
IS the same quantity in both and is what the comparison above uses.

## The general lesson

A number produced under one modelling assumption is not evidence about the
world until something varies that assumption. Residence looked like a physical
observation -- it has units, it varies between molecules, it identified
specific candidates -- and it was a property of a formula standing in for
water. The tell was available and unexamined: the metric had almost no dynamic
range (84% of candidates pinned at 1.00), and a measurement that is nearly
constant except for a few extremes is usually reporting on its own machinery.
