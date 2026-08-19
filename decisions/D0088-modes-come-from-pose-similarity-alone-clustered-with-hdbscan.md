---
id: D0088
title: Modes come from pose similarity alone, clustered with HDBSCAN — the old order clustered on the score and then scored the clusters
date: 2026-08-19
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
amends: D0086
affects:
  - shared/pose_cluster.py
  - shared/pose_modes.py
  - shared/pose_subsplit.py
  - exp/4_election/run_all.py
  - exp/5_mode_saturation/run_all.py
  - exp/6_mode_quality/run_all.py
evidence:
  - '@tt8804: "we first cluster into groups to make modes and then we rank the modes"; "there is only stage 2 -- what you are describing as stage 1 happens later"; "use HDBSCAN"'
  - 'THE CIRCULARITY: pose_modes.split clusters on the reactive atom POSITION and warhead DIRECTION. The viability score is distance(reactive atom, Cys113 SG) in 2.8-4.2 A plus an off-normal angle <= 30 deg. SG is fixed within a run, so position IS the distance term and direction IS the angle term -- the pipeline formed groups along the axis it then graded them on'
  - 'the code comment on that step asserts it clusters "never on the NAC geometry itself, which is the score", which is not what it does'
  - 'MODE QUALITY, 120 replicate x method rows over 3 molecules (exp/6_mode_quality), no scoring involved: widest mode anywhere -- shipped DBSCAN 9.30 A, complete 7.17, 1 A fine cut 11.06, HDBSCAN 3.91'
  - 'largest mode anywhere: shipped 137 poses, fine 208, HDBSCAN 14'
  - 'p90 within-mode width: shipped 7.86 A, fine 3.04, HDBSCAN 2.52'
  - 'THE MODE HOLDING THE VALIDATED POSE: shipped puts it in 108 poses spanning 8.03 A with viable fraction 0.26; HDBSCAN puts it in 8 poses spanning 1.50 A with viable fraction 0.72'
  - 'HDBSCAN settings measured on a real 393-pose cloud: leaf/min_cluster_size 3 gives median width 1.54 A and largest mode 14, against eom/min_samples 10 at 3.68 A and largest 64'
  - 'OPEN: HDBSCAN labels 29% of poses noise, and the validated pose itself was noise in 3 of 30 replicates'
  - 'OPEN: the mode holding the validated pose is still not pure -- viable fraction 0.72, so ~2 of its 8 poses do not reach attack geometry'
  - 'SATURATION (exp/5): mode count under the shipped rule is FLAT in sampling depth because min_samples is 5% OF THE SAMPLE -- a mode must hold 5% of all poses, so the bar rises with depth and no rarer mode is ever found; the count cannot exceed 20 by arithmetic'
  - 'with a FIXED threshold the shipped rule collapses instead: at 2000 poses DBSCAN merges everything into ONE group of 1,173 poses (eps 3.0 chains as the cloud densifies)'
  - 'sweeping eps with a fixed threshold, the log relationship @tt8804 predicted appears at eps 1.5 (R^2 0.98) and nowhere else -- 3.0 collapses, 1.0 and 0.5 grow linearly'
  - 'AutoDock-GPU ceiling on this build, measured: nrun 2000 fine; 5000 exits -6 "stack smashing detected" AND STILL WRITES A .dlg; 10000 exits 0 with no output. 100,000 is unreachable'
runbook: null
---

# D0088 — modes come from pose similarity alone

## The order was wrong, and it was circular

@tt8804: *"we first cluster into groups to make modes and then we rank the
modes"* — and, on being shown the implementation, *"there is only stage 2; what
you are describing as stage 1 happens later."*

The shipped pipeline does the opposite. It clusters on the reactive atom's
**position** and the **direction** the warhead faces, then subdivides those on
whole-molecule RMSD, then scores each group by what fraction of its poses reach
attack geometry.

But attack geometry IS position and direction. The distance term is the distance
from that atom to Cys113's SG, and SG is fixed within a run; the angle term is
that direction. So the groups were formed along the very axis they were then
graded on. The code comment claims the opposite — *"never on the NAC geometry
itself, which is the score"* — and it is wrong.

That circularity explains a day of dead ends: why a mode's viable fraction
behaved like a mixing ratio rather than a property, why tightening the RMSD
sub-split never bought homogeneity (the damage was done above it), and why the
first version of the election benchmark was itself circular in the same way.

## The design

One clustering step, on pose similarity alone — heavy-atom RMSD, no
superposition, nothing about the anchor. Attack geometry is used afterwards, to
rank the groups. `shared/pose_cluster.py`.

**HDBSCAN, not DBSCAN** (@tt8804). DBSCAN needs a radius chosen in advance and
links anything inside it, so a chain of poses each within `eps` of the next
becomes one group however wide it grows. Worse, the right radius moves with
sampling depth: at a fixed threshold and eps 3.0, a 2,000-pose cloud collapses
into a **single group of 1,173 poses**. HDBSCAN asks for a minimum group SIZE
instead, builds the hierarchy over all densities at once, and keeps the groups
that persist — so a tight group inside a diffuse halo survives instead of being
absorbed, and nothing needs retuning when depth changes.

## Measured, on modes only

120 (replicate x method) rows over 3 molecules, no scoring involved:

| | modes | noise | median width | p90 | widest | largest mode |
|---|---|---|---|---|---|---|
| shipped (DBSCAN) | 11.0 | 0% | 5.21 A | 7.86 | 9.30 | **137** |
| complete linkage | 14.4 | 63% | 2.25 A | 5.72 | 7.17 | 32 |
| 1 A fine cut | 38.8 | 0% | 0.90 A | 3.04 | 11.06 | **208** |
| **HDBSCAN** | 59.3 | 29% | 1.35 A | **2.52** | **3.91** | **14** |

HDBSCAN is the only rule that never produces a bag. And on the pose a 100 ns run
validated:

| | size | width | viable fraction |
|---|---|---|---|
| shipped | 108 poses | 8.03 A | 0.26 |
| **HDBSCAN** | **8 poses** | **1.50 A** | 0.72 |

The shipped pipeline puts the pose we know is right into a 108-pose group
spanning 8 A, three quarters of which cannot react. HDBSCAN puts it in an 8-pose
group 1.5 A wide, which is what "essentially the same pose" means.

## What is still open, and why this is `proposed`

- **29% of poses become noise**, and the validated pose was noise in 3 of 30
  replicates. `min_cluster_size` and the leaf/eom choice were tuned on ONE
  molecule's cloud; they need sweeping against noise rate and home-mode purity
  across all three before this becomes the default.
- **The home mode is not pure** (viable fraction 0.72): about 2 of its 8 poses
  do not reach attack geometry.
- **Adopting this invalidates nac_v5.** `consensus` is mode_size / n_poses, so
  changing how modes are cut changes every per-mode number downstream. It means
  a new topic and a full re-screen, which is a budget decision (#79).

## Two findings recorded here because they have nowhere else to live

**Mode count does not saturate, and the reason is a threshold that scales.**
`min_population_frac = 0.05` means a mode must hold 5% of ALL poses to exist —
25 poses at nrun 500, 100 at 2,000, 500 at 10,000. The bar rises with the
sample, so deeper docking never reveals a rarer mode, and the count cannot exceed
20 by arithmetic. @tt8804 predicted a logarithmic relationship; with a FIXED
threshold and eps swept, that curve appears at eps 1.5 (R^2 0.98) and at no other
value — 3.0 collapses to one group, 1.0 and 0.5 grow linearly.

**AutoDock-GPU has an undocumented ceiling on this build.** nrun 2,000 is fine.
At 5,000 it exits -6 with *"stack smashing detected"* **and still writes a
.dlg**. At 10,000 it reports "the job was not successful" and exits **0**. Both
shapes defeat a caller that trusts the exit code or the file's existence, and
`nac_screen.dock` now checks for all three. 100,000 poses is not reachable.
