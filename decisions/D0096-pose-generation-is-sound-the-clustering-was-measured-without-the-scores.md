---
id: D0096
title: Pose generation is sound — the exposed poses are real but score badly, and every clustering result so far was measured without the scores
date: 2026-08-26
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - integration/pose_group_viewer.py
  - scripts/persist_raw_clouds.py
  - scripts/nac_screen_v2.py
  - exp/21_pose_generation_audit/run_all.py
  - exp/21_pose_generation_audit/energy_filtered_grouping.py
  - exp/16_contact_clustering/run_all.py
  - exp/17_contact_saturation/run_all.py
  - exp/20_cloud_concentration/run_all.py
evidence:
  - '@tt8804: "there are literally poses outside of the pocket, how is that possibly the lowest energy. do a full audit of our pose generation step"'
  - 'THEY ARE NOT THE LOWEST ENERGY. On the molecule the viewer opens by default, poses with >30% of their atoms uncontacted are 13 of 500 (2.6%), sit at the 88th ENERGY PERCENTILE, and ZERO of them are in the best decile'
  - 'the energy-geometry relationship has the RIGHT SIGN and is strong: rho(contacts, energy) = -0.590, rho(enclosure, energy) = -0.588, rho(exposed fraction, energy) = +0.446'
  - 'best energy decile: 55 receptor contacts, 4% of atoms uncontacted. Worst decile: 40 contacts, 18% uncontacted. Monotone across all ten deciles'
  - 'PoseBusters passes 99.2-99.4% of poses; no pose is outside the receptor (0% with <10 contacts, 0% beyond 15 A of Cys113 SG across 4,000 poses of 8 molecules)'
  - 'THE CLOUDS CARRIED NO ENERGIES AT ALL: nac_screen_v2 and persist_raw_clouds wrote coordinates with no SD tags, so exp/16, exp/17, exp/19 and exp/20 all weighted the best pose and the 500th equally, and so did the viewer'
  - 'ENERGY-SELECTED POSES AGREE WITH EACH OTHER: keeping the best 25% raises the top-1 group share from 2.8% to 10.4% and drops the effective pose count from 133 to 23'
  - 'AGAINST A SIZE-MATCHED RANDOM CONTROL (the check that makes it a finding rather than an arithmetic artefact of smaller n): best-25% concentrates 2.60x more than a random 25% of the same cloud, best-10% 3.33x, in 21 of 21 molecules, Wilcoxon p = 5.9e-05'
  - 'the saturation exponent falls but does not vanish: b = +0.767 over all poses against +0.566 within the best 25%'
  - 'the energy landscape is flat, confirming exp/12: span 3.5-4.6 kcal/mol, and 76% of poses within 2 kcal/mol of the best on the six-molecule panel'
  - 'THE POSE-ENERGY PAIRING IS SOLVED, NOT ASSUMED: records and conformers matched by an order-invariant signature under Hungarian assignment, returning the identity permutation at 0.00000 A median error'
runbook: python exp/21_pose_generation_audit/run_all.py; python exp/21_pose_generation_audit/energy_filtered_grouping.py
---

# D0096 — the docking was fine; the analysis was blind

## The observation

@tt8804, in the pose viewer: *"there are literally poses outside of the pocket,
how is that possibly the lowest energy."*

## The answer

**It is not the lowest energy.** Those poses exist — they are 2.6% of the cloud —
and the scoring function ranks them near the bottom:

| energy decile | receptor contacts | % of ligand uncontacted |
|---|---:|---:|
| best 10% | 55 | 4% |
| worst 10% | 40 | 18% |

ρ(exposure, energy) = **+0.446**, ρ(contacts, energy) = **−0.590** — the expected
signs, monotone across all ten deciles. Of the 13 most exposed poses, **zero** are
in the best energy decile; their median energy percentile is **88**. PoseBusters
passes 99.2%. Nothing is outside the receptor: across 4,000 poses of 8 molecules,
0% have fewer than 10 receptor contacts and 0% sit beyond 15 Å of Cys113.

**Pose generation is sound.** The defect is that nothing ever said so.

## What was actually wrong

**The persisted clouds carry no energies.** `nac_screen_v2` writes coordinates and
a mode label; `persist_raw_clouds` copied that. So a pose the scorer put 440th of
500 was written, stored, clustered and displayed exactly like the best one.

That has two consequences, and the second is worse than the first:

1. **The viewer invited the wrong conclusion.** It drew 500 poses identically with
   no score anywhere on screen, and defaulted the receptor to a cartoon with the
   surface off — against `pose3d.py`'s own founding note, which says in as many
   words that *"a spectrum cartoon tells a reader where the chain runs; it does
   not tell them where the ligand IS"*. Pin1's site is a shallow surface groove.
2. **Every clustering experiment was energy-blind.** exp/16, exp/17, exp/19 and
   exp/20 all measured the whole cloud, tail included.

## What that cost, measured

Re-persisting the clouds with energies and re-running the grouping:

| kept by energy | poses | groups | top-1 share | effective # poses | singletons |
|---|---:|---:|---:|---:|---:|
| 100% | 500 | 218 | 2.8% | 133 | 48% |
| 50% | 250 | 88 | 6.0% | 51 | 35% |
| 25% | 125 | 41 | 10.4% | 23 | 27% |
| 10% | 50 | 17 | 20.0% | 9 | 30% |

**And this is not an arithmetic artefact of a smaller sample** — which is exactly
what it would look like, and what the first version of this analysis claimed
without testing. Against a **random subset of identical size** from the same
cloud:

| kept | top-1 best | top-1 random | ratio |
|---|---:|---:|---:|
| 50% | 6.0% | 3.6% | 1.67× |
| 25% | 10.4% | 4.0% | **2.60×** |
| 10% | 20.0% | 6.0% | **3.33×** |

21 of 21 molecules at every level, Wilcoxon **p = 5.9e-05**.

**Low-energy poses genuinely agree with each other.** That is a positive result
about the scoring function, and this project has not had many.

## What has to be re-read

* **exp/20's "no molecule has a consensus pose" is too strong.** Measured over the
  whole cloud the top group holds 1.4–7.2%; within the best 10% by energy it holds
  20%, with an effective pose count of 9. Still not one dominant pose — but a very
  different picture, and the honest statement is *"the clouds are far more
  converged among the poses the score favours than the full-cloud numbers show."*
* **exp/17's saturation exponent softens**: b = +0.767 over all poses against
  **+0.566** within the best 25%. The count still climbs, so D0092's conclusion
  (the count is a sampling statistic, never a mode count) stands — but part of the
  climb was the tail.
* **exp/16's group quality is unaffected.** Complete linkage bounds within-group
  distance by the tolerance regardless of which poses enter.

## Fixes

1. `persist_raw_clouds.py` writes `free_energy_kcal` per pose, with the
   pose↔energy pairing **solved rather than assumed** — matched by an
   order-invariant signature under Hungarian assignment, required to return the
   identity permutation at ~0 Å. AutoDock reports a cluster ranking beside the run
   order, so "the order" is genuinely ambiguous here.
2. The viewer shows the energy range, offers a best-N% filter, reports each
   group's best energy, **defaults the pocket surface ON**, and warns loudly on any
   cloud that carries no energies.
3. `--skip-existing` was `store_true` with `default=True` — a flag that could not
   be turned off, and it silently kept every already-persisted cloud energy-less
   after the writer changed. Replaced with `--force`.

## Why it looked right

Every number was real. The poses were real docking output, the geometry was fine,
PoseBusters passed them, and the group statistics were correctly computed over
exactly the set they were given. Nothing was broken — the analysis was simply
answering "what does the whole cloud look like" while everyone read it as "what
does the docking think". The two differ by the 75% of poses the score rejects.

## Guard

`energies_aligned` raises unless the pose↔energy assignment is the identity at
near-zero error. `energy_filtered_grouping.py` reports every filtered statistic
beside a size-matched random control, because a concentration that rises when n
falls is not evidence of anything.
