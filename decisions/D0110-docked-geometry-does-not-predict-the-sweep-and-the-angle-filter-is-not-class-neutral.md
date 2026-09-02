---
id: D0110
title: Docked geometry does not predict the sweep outcome, post-equilibration geometry does, and an angle filter on a distance-selected set is not class-neutral
date: 2026-09-02
status: accepted
approach: shared
decided_by: '@twu383'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/attack_sweep.py
  - scripts/md_residence_3ikd.py
  - scripts/sweep_supervisor.py
  - config/target.yaml
evidence:
  - 'first five nac_v8 sweeps: DOCKED warhead-SG distance 2.78-2.97 A for all five, frac_attack_ready 0.000 to 0.926 -- the docked geometry has no discriminating power over the range it was selected on'
  - 'the same five ordered by POST-EQUILIBRATION distance (`start_dist_a`, first production frame after 300 ps unrestrained NVT/NPT): 3.58 -> 0.926, 3.76 -> 0.264, 4.40 -> 0.215, 5.73 -> 0.000, 6.46 -> 0.000. Monotonic'
  - 'the docking and MD frames are IDENTICAL: `receptor_cys.pdb` and `3IKD_noligand.pdb` both put Cys113 SG at (13.385, 3.989, -2.040); computing the docked distance against that static sulfur reproduces the worklist values. The drift is physics, not a coordinate mismatch'
  - 'off-normal angle over the WHOLE analysed cloud is near-identical across classes: acrylamide 54.8, bdhi_c4 48.5, bdhi_c5 51.6 deg (isotropic expectation 60)'
  - 'in the 2.8-3.0 A shell it diverges: acrylamide 48.1, bdhi_c4 8.8, bdhi_c5 11.7 deg'
  - 'and it is a distance gradient for BDHI only -- bdhi_c4 by shell: 8.8 (2.8-3.0), 28.2 (3.0-3.5), 47.8 (3.5-4.2), 67.9 (4.2-6.0); acrylamide barely moves: 48.1, 49.7, 58.5, 67.3'
  - 'so a 45 deg cut on the sub-3 A set removes 4,193 of 8,488 modes: acrylamide 7,597 -> 3,404 (-55%), bdhi_c4 649 -> 649 (0%), bdhi_c5 242 -> 242 (0%)'
  - 'BDHI mechanism is `sn2_ring_opening` = perpendicular_to_plane, same as acrylamide''s `michael_addition`, and both carry isotropic_null 0.0816 -- so the angle IS the same physical quantity and the comparison is legitimate; what differs is the freedom each warhead has to adopt it'
runbook: null
---

# Two findings about the geometry criteria, both of which change how the sweep is read

## 1. The docked geometry does not predict the sweep; what happens in the first 300 ps does

The nac_v8 worklist selects modes on the median warhead-to-Cys113 distance of
their member poses. The first five completed sweeps all sat at **2.78–2.97 Å**
docked and returned frac_attack_ready from **0.000 to 0.926**. Over the range it
was selected on, the docked distance carries no signal at all.

What separates them is the distance at the **start of production** — after the
300 ps of unrestrained NVT/NPT every sweep runs first:

| post-equilibration | frac_attack_ready |
|---|---|
| 3.58 Å | 0.926 |
| 3.76 | 0.264 |
| 4.40 | 0.215 |
| 5.73 | 0.000 |
| 6.46 | 0.000 |

Monotonic on five points. That quantity is the **tier-1 readout** — the one
measurement D0071 validated (p = 0.007) and D0108 rested a NO GO on — arriving
here as a free by-product of a sweep that was going to run anyway.

**It is not a frame artefact, and that was checked first.** The obvious
explanation is that docking treats Cys113 as flexible while MD does not, so the
two measure to different sulfurs. They do not: `receptor_cys.pdb` in a sweep
workdir and the docking receptor `3IKD_noligand.pdb` both place Cys113 SG at
(13.385, 3.989, −2.040), and recomputing the docked distance against that static
sulfur reproduces the worklist numbers. The drift is the pose failing to survive
plain dynamics.

### Consequence

An **early give-up** is now available and implemented: measure the warhead after
equilibration and skip the remaining 1,200 ps when it has clearly gone.
Equilibration is 300 of the 1,500 ps a sweep runs, so a departed pose costs a
fifth of a full sweep instead of all of it.

`--abort-above-a`, default **0 (off)** in `attack_sweep` and passed explicitly as
**6.0** by `sweep_supervisor`. The default is off because an abort DISCARDS work
and must be asked for; the threshold lives with the campaign rather than in a
library default another caller could inherit. 6.0 Å is past the 4.2 Å window and
past every pose that has scored above zero, so it catches only the unambiguous
cases. `start_dist_a` accumulates on every completed row, so the cut can be
re-derived from hundreds of sweeps rather than these five.

### What it cost to get one trustworthy number

Three faults, in sequence, and the first was live on the campaign:

1. A naive `.gro` read ignored periodic images and reported **51.18 Å** for a
   ligand in a ~7 nm box — the box length minus the real distance. It gave up on
   one real mode (`t4_b49ffa60a11a_m113`) at a reported 56.6 Å before it was
   caught. That row is invalidated and the mode re-queued.
2. Fixing the wrap by matching `resname == "CYS"` found **two** sulfurs — 3IKD
   has Cys57 as well as Cys113, both reduced in the MD system — and correctly
   refused to choose, so it never measured anything.
3. Selecting Cys113 by residue number through `md_movie.PIN1_OFFSET`, with the
   residue's identity verified rather than trusted, gives **4.56 Å** on the same
   frame that first read 51 Å.

The measurement is **fail-safe by construction**: any uncertainty — missing
frame, triclinic box, atom-count mismatch, ambiguous sulfur, a result beyond half
the box diagonal — returns None and the full sweep runs. Skipping work throws a
molecule away; doing it costs six minutes.

## 2. An angle filter on a distance-selected set is not class-neutral

`sweep_rule` now trims the worklist at **45° off-normal** (@twu383: *"i dont care
about being equal between classes. i just want to trim bad candidates"*). The
effect is not symmetric and the record needs to say why, because the trimmed
worklist otherwise reads as evidence that BDHI aligns better than acrylamide.

| | modes before | after 45° |
|---|---|---|
| acrylamide | 7,597 | 3,404 (−55%) |
| bdhi_c4 | 649 | **649 (−0%)** |
| bdhi_c5 | 242 | **242 (−0%)** |

**Over the whole cloud the three classes are nearly identical** — acrylamide
54.8°, bdhi_c4 48.5°, bdhi_c5 51.6°, against an isotropic expectation of 60°.
The gap appears only inside the shell the worklist already selected:

| median off-normal | 2.8–3.0 Å | 3.0–3.5 | 3.5–4.2 | 4.2–6.0 |
|---|---|---|---|---|
| acrylamide | 48.1 | 49.7 | 58.5 | 67.3 |
| bdhi_c4 | **8.8** | 28.2 | 47.8 | 67.9 |
| bdhi_c5 | **11.7** | 26.1 | 45.8 | 68.6 |

BDHI's angle collapses from 68° to 9° as the sulfur closes in; acrylamide's
barely moves. **That is sterics, not reaction competence.** BDHI's reactive
carbon sits inside a rigid 3-bromo-4,5-dihydroisoxazole ring, so at 2.8–3.0 Å
the only place a sulfur fits without clashing into ring atoms is over the ring
face — perpendicular by necessity. Acrylamide's vinyl is terminal and exposed, so
a sulfur can sit in-plane at close range.

**The angle is the same physical quantity for both** — both mechanisms map to
`perpendicular_to_plane` and both carry `isotropic_null = 0.0816`, so the
comparison is legitimate. What differs is the *freedom each warhead has to adopt
it*. On a set already cut at <3 Å, distance and angle are geometrically coupled
for BDHI and largely independent for acrylamide.

This is D0088's circularity in a new place: filtering on one axis and then
grading on an axis that is not independent of it — for one of the chemistries.

### What must not be concluded from the trimmed worklist

* **Not** that BDHI reaches better attack geometry than acrylamide. Measured
  where neither warhead is constrained, they are the same.
* **Not** that the surviving BDHI modes are better than the trimmed acrylamide
  ones. They survived a filter their geometry cannot fail.
* The campaign is now more BDHI-weighted than the library is, and BDHI still has
  **zero** crystallographic Cys113 positives (D0108). The rebalancing is a
  consequence of the filter, not a finding about the chemistry.

### If class balance ever matters

Apply the cut **within** each class — keep the best N% by angle per warhead — for
the same time saving with the split preserved. Rejected here deliberately, on the
grounds above.
