---
id: D0092
title: Residue-contact space is fixed; the group count climbs because 6,000 poses undersample it, and the groups themselves never move
date: 2026-08-26
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - exp/17_contact_saturation/run_all.py
  - exp/17_contact_saturation/tolerance_sweep.py
  - exp/17_contact_saturation/persistence.py
  - exp/17_contact_saturation/space_growth.py
  - shared/pose_contacts.py
  - docs/build_plan_next.md
evidence:
  - '@tt8804: "lets see the groups as a function of poses generated like before. if it tapers off we just accept the huge number and go to ranking"'
  - '@tt8804: "also can we check if this residue contact space grows with poses"'
  - 'THE COUNT DOES NOT TAPER: 6,000-pose raw cloud at the RMSF tolerance 0.73 A gives 289 groups at n=500 and 1,438 at n=6,000 -- b = +0.693, and the species-accumulation fit implies NO finite plateau'
  - 'nor at any other tolerance that keeps groups tight: b = +0.567 at 1.0 A, +0.446 at 1.5 A, no finite plateau at any of the seven tested'
  - 'LOOSENING THE TOLERANCE REBUILDS THE BAG: largest group 37 poses at 0.73 A, 126 at 1.0 A, 433 at 2.0 A, 2,004 at 3.0 A, 3,147 at the 3.5 A sweep bar (which yields 3 groups total)'
  - 'THE SPACE ITSELF IS FIXED: diameter exponent +0.019, mean pairwise separation exponent +0.001; 60x more poses widened the diameter 1.10x (a max-statistic artefact) and moved the mean 1.00x; the 99th percentile is 3.11 A at every depth from n=500'
  - 'EFFECTIVE DIMENSION ~3.5 out of 420 available coordinates (28 atoms x 15 residues), from the covering-number slope N(eps) ~ eps^-d; still rising with n (2.24 at n=100, 3.54 at n=6,000) so it is itself undersampled'
  - 'GROUPS NEVER MOVE: 100% of shallow (n=500) groups have an n=6,000 counterpart within the tolerance, over 5 draws and 1,460 shallow groups; median centre displacement 0.254 A against a 0.73 A tolerance; 100% for non-singletons and for groups of >= 5'
  - 'GROWTH IS ALL TAIL: of 1,438 deep groups, 26% are singletons and 61% hold <= 3 poses; 65% of the cloud lives in the 431 groups of >= 5'
  - 'consistent across chemistry: exponent b median +0.668 over 12 production molecules, range +0.589 to +0.780, none above 0.9'
  - 'coverage cost at 6,000 poses: ~10,400 poses to cover the visible region at 0.73 A, ~4,500 at 1.0 A, ~1,100 at 1.5 A, ~266 at 2.0 A'
runbook: python exp/17_contact_saturation/run_all.py; python exp/17_contact_saturation/tolerance_sweep.py; python exp/17_contact_saturation/persistence.py; python exp/17_contact_saturation/space_growth.py
---

# D0092 — contact space is fixed, the sampling is not

## The question

exp/16 showed contact-space grouping produces tight groups (median 1.08 Å
within-group Cartesian RMSD, worst 2.71 Å) where both predecessors produced bags.
It also produced ~174 groups per molecule. The proposal on the table was: if the
count tapers with docking depth, accept it and move to ranking.

## What is entailed, and therefore is not a finding

Every coordinate of contact space is a distance capped at `pose_contacts.CAP_A`
(10 Å) and floored at van der Waals contact. **The region is bounded before a
single pose is docked.** Complete linkage at a fixed absolute tolerance likewise
cannot produce more groups than the covering number at that scale. Reporting
either as a result would be D0091 in a new coordinate system — docking into a
26 Å box, measuring that the cloud stayed inside 26 Å, and calling the bound a
finding. The findings below are rates and dimensions, never bounds.

## The count does not taper — at any usable tolerance

| tolerance | n=500 | n=6,000 | growth | exponent b | largest group |
|---:|---:|---:|---:|---:|---:|
| 0.73 Å (RMSF) | 289 | 1,438 | 4.97× | +0.693 | 37 |
| 1.00 Å | 185 | 687 | 3.71× | +0.567 | 126 |
| 1.50 Å | 76 | 209 | 2.75× | +0.446 | 194 |
| 2.00 Å | 32 | 64 | 2.00× | +0.337 | 433 |
| 3.00 Å | 5 | 9 | 1.93× | +0.246 | 2,004 |
| 3.50 Å (sweep bar) | 1 | 3 | 2.25× | +0.300 | 3,147 |

No row admits a finite plateau. And the tolerance is not a free dial: by 2.0 Å the
largest group is 433 poses and by 3.0 Å it is 2,004 — **the bag D0088 and D0091
were both written to condemn.** There is no tolerance at which the count
saturates *and* the groups stay tight.

The exponents for the loosest two rows are fitted to counts of 1→3 and 5→9 and
say nothing about saturation; those tolerances are excluded on group width, not
on growth.

## But the space is fixed, and that changes what the climb means

The occupied region stops changing almost immediately:

* diameter exponent **+0.019**, mean pairwise separation exponent **+0.001**
* 60× more poses widened the diameter by 1.10× — and the maximum of a larger
  sample is larger by construction, so even that is an upper bound on the effect
* the 99th percentile of pairwise distance is **3.11 Å at every depth from 500 to
  6,000**, and the mean is 2.23 Å at every depth

So the climbing group count is **undersampling of a fixed region**, not expansion
of the region. Those are different problems: expansion would mean the search is
unbounded and the pose set is an artefact of runtime; undersampling means the
region is a property of the molecule and we have not covered it.

## The poses use ~3.5 of 420 available dimensions

From the covering-number slope N(ε) ~ ε^−d, the effective dimension is **3.54** at
n=6,000, against 420 coordinates offered (28 heavy atoms × 15 landmark residues).
That is the right order for rigid-body placement, and it is the evidence that the
contact metric tracks pose rather than noise — a metric measuring noise would
report a dimension near 420. The estimate is still rising with n (2.24 → 3.54), so
it too is undersampled and 3.54 is a floor.

## The groups themselves never move — this is the load-bearing result

Because the region is fixed and the tolerance is absolute, a deeper cloud can only
**add** groups; it cannot relocate the ones already found. Measured rather than
assumed, over 5 draws and 1,460 shallow groups:

* **100%** of n=500 groups have an n=6,000 counterpart within the tolerance
* median centre displacement **0.254 Å** against a 0.73 Å tolerance
* 100% holds for non-singletons alone, and for groups of ≥ 5 poses alone

and the growth is entirely in the sparse tail: 26% of deep groups are singletons,
61% hold ≤ 3 poses, while 65% of the cloud sits in the 431 groups of ≥ 5.

This is the property HDBSCAN lacked. D0088's rule lost the MD-validated pose in 3
of 30 replicates and only 1 of 3 modes survived a re-dock, because it had no
length scale and re-derived its clusters from whatever density the sample
happened to show. Here the tolerance is absolute and molecule-owned, so identity
is a property of the region, not of the draw.

## What this decides

**The group count is a sampling-density statistic and must never be reported as a
number of binding modes.** It is a monotone function of docking depth, it has no
plateau, and a shortlist that ranks "modes" invites exactly the reading D0088
found in the shipped pipeline.

Given that, and given the groups are stable, grouping is admissible **as a
cost-saving collapse of the pose cloud** — @tt8804's own reframing — and the
pipeline may proceed to ranking. Two conditions:

1. the artefacts say *groups*, never *modes*, and never report the count as a
   property of the molecule;
2. ranking is compared across a **fixed** docking depth, because the group
   population is depth-dependent even though each group is not.

## Why the alternative looked right

"The count will taper" is the correct intuition for a bounded region and a fixed
length scale — and it is true asymptotically. It fails here for a quantitative
reason that no amount of reasoning about boundedness would surface: at 0.73 Å the
region needs on the order of 10,400 poses to cover what 6,000 have already
revealed, and the true covering number is larger still, since 1,438 is what 6,000
poses found rather than what is there. We are on the rising part of a curve whose
plateau is real and out of reach. Boundedness was never the question; the constant
was.

## Caveats

* **The extent, dimension and persistence results are one molecule**
  (`t4_716800c125a7`), the only 6,000-pose cloud that exists. The exponent is
  replicated across 12 production molecules (median +0.668) but the ladder there
  stops at ~450 poses.
* The coverage estimates use N·ln N, which assumes uniformly weighted cells; the
  cells are plainly not uniform (65% of the cloud in 30% of the groups), so they
  are an order of magnitude, not a spec.
* Every production cloud read here is pre-filtered — see D0093.
