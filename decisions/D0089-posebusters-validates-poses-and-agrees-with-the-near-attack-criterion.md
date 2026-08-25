---
id: D0089
title: PoseBusters gates pose validity, protects attack-ready poses rather than filtering them, and does not change the ranking
date: 2026-08-17
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/nac_screen_v2.py
  - scripts/nac_screen.py
  - config/target.yaml
  - docs/build_plan_next.md
evidence:
  - '@tt8804: "add posebusters to this step since it is quite lightweight ... run docking checked by posebusters to get up to 500 posebuster stable poses as a quota"'
  - 'posebusters 0.6.5 already installed in dwi_cheminf; no new environment'
  - 'PASS RATE 90.60% over 1,500 poses, 25 molecules from each of acrylamide / bdhi_c4 / bdhi_c5 (91.6 / 88.2 / 92.0); per molecule min 70%, median 90%, max 100%'
  - 'ONLY 2 OF 22 CHECKS EVER FAIL: minimum_distance_to_protein 7.53%, internal_energy 2.00%; the other twenty 0.00%'
  - 'the twenty cannot fail: AutoDock varies rigid-body placement and torsions only, so bond lengths, angles and ring geometry are inherited from the input conformer -- they test the conformer generator, not the docking'
  - 'ATTACK-READY POSES PASS MORE, NOT LESS: 98.57% (n=1,819) against 92.80% (n=4,015), odds ratio 5.35, Fisher exact p = 1.9e-23'
  - 'the clash check fails 0.93% of attack-ready poses against 6.48% of the rest -- it rejects 7x more of what we do not want'
  - 'by distance to SG: <2.8 A 45.2% pass (54.8% clash-fail); 2.8-3.2 A 89.5%; 3.2-3.6 A 96.2%; 3.6-4.2 A 96.5%'
  - 'RANKING UNCHANGED: per-mode viable_fraction recomputed with failures removed over 75 modes >= 12 poses -- median change +0.00 pp, mean +0.67 pp, rho = 0.9989, top-10 overlap 10/10, top-25 25/25, 19 of 75 modes move at all with median move 0 places'
  - 'cost: 232 ms/pose sustained, ~20 CPU-h per screen, ~25 min on the 50-worker cap; docking +10.4% on the mean (552 runs for 500 valid), +43% for the worst molecule sampled (70% pass)'
  - 'PB clash threshold is 0.75 x vdW sum = 2.625 A for reactive C vs Cys113 SG, BELOW the 2.8 A window floor -- which is why it bites the too-close poses rather than the attack-ready ones'
runbook: null
---

# D0089 — PoseBusters as a validity gate

## Context

@tt8804 proposed adding PoseBusters to the docking step and moving the quota
from **500 runs** to **500 PoseBusters-valid poses**. Two questions had to be
answered before it could be adopted: what it costs, and whether it filters out
the poses the pipeline exists to find.

## Decision

**Adopt the gate and the quota — for validity, reproducibility and an equal
denominator, and NOT because it improves the ranking, because it measurably does
not.**

That qualification is the decision. Without it the release notes would claim an
improvement the measurement refuses.

## What was measured

**1. It does not filter out what we want. It does the opposite.**

Attack-ready poses pass at **98.57%** against **92.80%** for the rest — odds
ratio **5.35**, Fisher exact **p = 1.9 × 10⁻²³**. The clash check, the one
expected to do the damage, fails **0.93%** of attack-ready poses against
**6.48%** of the others.

**2. Why the concern was wrong.** PoseBusters flags a clash below 0.75 × the sum
of van der Waals radii, which for the reactive carbon against Cys113 SG is
**2.625 Å** — *below* the near-attack window's 2.8 Å floor. So it rejects poses
that are **too close**: below 2.8 Å only 45% survive, and through the window
89–96% do. PoseBusters and the near-attack criterion agree, independently, about
which poses are physically real.

**3. It does not change the ranking.** ρ = **0.9989** on per-mode
`viable_fraction`, top-25 overlap **25/25**, median rank movement **zero
places**.

## Why it looked like it would help

The proposal is intuitive and the intuition is sound in general: docking
programs do produce physically impossible poses, and 9.4% of ours are. What the
intuition misses is that **the near-attack criterion was already rejecting
almost all of them** — a pose that clashes into the protein is usually not
presenting its warhead at 2.8–4.2 Å and 150°. The two filters overlap far more
than they complement, so removing the invalid poses barely moves a score
computed over the survivors.

The second thing it misses is what PoseBusters can see here. Twenty of its
twenty-two checks test bond lengths, angles and ring geometry — all inherited
unchanged from the input conformer, because AutoDock varies only rigid-body
placement and torsions. **Those checks certify our conformer generator, not our
docking**, and they pass at 100% until ligand prep changes. Describing this as
adding "22 physical checks" would overstate it by twenty.

## Consequences

* **Adopted for three reasons that survive the null result**: the quota gives
  every molecule an equal denominator (molecules are currently ranked against
  each other on 418–486 poses with nothing recording it); the twenty always-
  passing checks are standing insurance against a ligand-prep regression; and
  "PoseBusters-validated" answers a reviewer question #66 will otherwise attract.
* **The equal denominator comes from the QUOTA, not from PoseBusters.** Quotaing
  on raw pose count delivers it for zero CPU. If GPU time is contested, that is
  the cheap version of this change and it should be said out loud rather than
  bundled.
* **Cost**: +10.4% GPU on the library mean, +43% worst case; ~20 CPU-h for the
  gate itself. The per-molecule spread means a flat over-request sized on the
  mean undershoots for roughly half the library — see `build_plan_next.md` §1.8.
* **A molecule that cannot reach the quota must be stamped, not silently short**
  (`quota_met`, `n_runs_requested`, `n_poses_valid`, seeds). Same rule as
  `docked_species_ok` (#58).
* **Reportable in its own right**: 9.4% of raw AutoDock poses on this target are
  physically invalid, and they are disproportionately the non-reactive ones.
* **The sample for the attack-ready test was biased on purpose** — the 13
  molecules with the most viable poses, where the risk would show most sharply.
  It is the wrong sample for a library-wide effect; 90.6% (random, across
  families) is the unbiased pass rate, and a full rerank comparison belongs in
  the re-screen.
* **Today's clouds are not a valid baseline** for that comparison: 105 of 561
  molecules (18.7%) carry a persisted cloud that does not match their own scored
  rows, from an `if not adest.exists()` cache since fixed in `a7cca45`.
