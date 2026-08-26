---
id: D0091
title: The volume-partition proposal is refuted — the search box entails the bound, and a 3 Å partition rebuilds the bag it was meant to remove
date: 2026-08-17
status: proposed
approach: shared
decided_by: '@tt8804'
origin: adversary
supersedes: []
superseded_by: null
affects:
  - exp/13_pose_volume/run_all.py
  - docs/build_plan_next.md
  - decisions/D0090-the-pose-cloud-does-not-saturate-because-the-energy-landscape-is-flat.md
evidence:
  - '@tt8804: "my intuition is that the use of pose busters and the nature of docking sets a bounded pose space contained within the receptor volume and we can partition this space by 3a and assign poses per partition"'
  - 'THE CLOUD FILLS THE BOX: heavy-atom extent 25.48 x 25.42 x 25.16 A inside a 26 A cube -- 97-98% of every axis (verified independently of the audit)'
  - 'AutoDock hard-clips to the grid, so no pose can place an atom outside the box; the volume question had one possible answer before any pose was docked'
  - 'A 3 A CENTROID PARTITION: 49 occupied cells, largest holds 1,758 of 6,000 poses (29%), with median 6.08 A and max 9.13 A RMSD INSIDE a cell 3.0 A wide (verified)'
  - 'for comparison D0088 condemned the shipped rule at 137 poses spanning 9.3 A -- the 3 A partition is 13x the membership at the same width'
  - 'two poses may share a centroid and be 180 deg flips: a partition of R^3 cannot bound a distance in configuration space'
  - 'THE EXPONENT CONTRAST DISSOLVES: fitting local slopes rather than one OLS line over a curved log-log plot, and correcting rarefaction bias, gives volume b ~ 0.32 against cover b ~ 0.325'
  - 'the rarefaction bias is ASYMMETRIC: subsampling one pooled 6,000-pose cloud inflates the volume ladder low rungs ~19% while leaving the covering number untouched (177 vs 177) -- it manufactured the contrast rather than adding noise'
  - 'the quoted "14% of the box" used POINT occupancy (2,532 A^3, verified 14.4%), the column exp/13 own docstring disowns as boundary-sensitive; sphere occupancy gives ~34%'
  - 'the quoted "1.71x poses -> 1.029x cells" was a single unreplicated draw; over 20 draws the mean is 1.102x'
runbook: null
---

# D0091 — the volume partition is refuted

## What was proposed

Stop grouping poses by similarity; define the 3D volume the poses occupy, partition
it at ~3 Å, and assign each pose to a cell. The motivation was sound: HDBSCAN has
no length scale (D0090) so its mode count grows as n^0.98, whereas a physical
volume is bounded by the pocket and should give a count that stops moving.

## Why the supporting measurement carried no information

**The pose cloud fills the docking box.** Heavy-atom extent is
**25.48 × 25.42 × 25.16 Å** inside a **26 Å cube** — 97–98% of every axis, all six
walls touched. AutoDock hard-clips to the grid; no pose *can* place an atom
outside it.

So "does the occupied volume grow without bound?" had **one possible answer,
fixed by the apparatus before anything was docked.** The experiment could not
have returned anything else.

This is catalogue entry #12 — `resi 101..125`, a sequence window standing in for
the measured pocket shell — in a new coordinate system: **the search box standing
in for the receptor cavity.** The claim's own wording says "bounded pose space
*contained within the receptor volume*", and the envelope measured is the box,
not the receptor.

**And PoseBusters is doing almost none of the bounding it was credited with.** It
removes ~10% of poses (D0089); the box removes everything else.

## Why the partition fails even granting the bound

Assigning a pose to a cell needs one point to stand for it, and the centroid is
the only candidate the proposal offers. Measured over 6,000 poses:

| | 3.0 Å partition |
|---|---:|
| occupied cells | 49 |
| largest cell | **1,758 poses (29% of the cloud)** |
| median RMSD *inside* that cell | **6.08 Å** |
| max RMSD inside that cell | **9.13 Å** |

A cell 3 Å wide containing poses 9 Å apart. **D0088 condemned the shipped rule
for putting 137 poses in a 9.3 Å bag; this is 13× the membership at the same
width** — the proposal rebuilds, worse, the artefact the whole redesign exists to
remove.

The cause is geometric and no cell size fixes it: **two poses can share a centroid
and be 180° flips of each other.** A partition of ℝ³ bounds a distance in ℝ³; it
says nothing about distance in configuration space, which is what "the same pose"
means.

## Why the numbers looked better than they were

Two errors, both in the direction that favoured the hypothesis.

**A single OLS slope was fitted to a curved log-log plot.** Growth is fast early
and flat late, so one exponent averages over systematic curvature — and it
compared volume's late behaviour against the covering number's whole-range
average.

**Every ladder rung was a subsample of one pooled 6,000-pose cloud**, not an
independent docking. A random 500 drawn from 6,000 is already spread across
everything in the pile and looks more complete than a fresh 500-run docking
would. Measured: this inflates the volume ladder's low rungs by ~19% **while
leaving the covering number untouched** (177 vs 177). The bias did not add noise
to the comparison — **it created it.**

Corrected and put on the same footing: **volume b ≈ 0.32, covering number
b ≈ 0.325.** There is no measured difference.

## Two reported numbers were wrong

* **"14% of the box"** used point occupancy — a voxel counted only if an atom
  *centre* fell in it — which `exp/13`'s own docstring names as the
  boundary-sensitive one that overstates growth. The physical (van der Waals
  sphere) figure is ~34%. Catalogue disguise #1, on a live pair the document
  already lists.
* **"1.71× poses → 1.029× cells"** was one unreplicated draw. Over 20 draws the
  mean is 1.102×.

## Why it looked right

The prediction was "bounded", the measurement said "bounded", and the agreement
was produced by the apparatus rather than by the pocket.

`how_this_project_breaks.md` records that the commonest detection route is
*"someone looked at output and it didn't match expectation"* — 9 of 25. **Here
the output matched expectation exactly.** There was no way to notice by looking,
which is precisely the case the document says a positive control is for. None was
used, and one was free: the five replicate clouds of §2.4i have an answer already
on record.

## What survives

> At 6,000 poses the cloud has swept roughly a third of its own docking box and
> is **still adding territory at b ≈ 0.32 by every metric measured.**

A diminishing-returns statement, consistent with D0090. It does **not** license a
fixed partition count, and it is not evidence that the pocket bounds anything.

## Consequences

* **The volume framing is dropped** unless the box test below reverses it.
* **The one measurement that would carry information**: re-dock with the box at
  1.5× and 2×. If the swept envelope expands proportionally, the bound was always
  the box. If it does not, there is a real physical boundary and the intuition
  was right — tested in a way that could show it. Minutes of GPU (D0090 records
  3 × 2,000 runs at 32 s).
* **Any future ladder must be built from independent dockings at each depth**,
  not subsamples of one cloud.
* **Any partition intended to group poses must be built in the pose metric**, and
  must report within-partition RMSD spread as a required output that can fail.
* **`exp/13` needs a positive control** before it is trusted again.
