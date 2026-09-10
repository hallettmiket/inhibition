---
id: D0118
title: The best 100 ns result does not reproduce — m50 gives 35.7% and 0.0% on two replicates
date: 2026-09-10
status: accepted
approach: shared
decided_by: '@twu383'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - config/presentation.yaml
  - docs/state_of_the_project.md
evidence:
  - 't4_7b02d1dc4fd2_m50 rep1: frac_attack_ready=0.357285, ligand_pose.sdf md5 73688634b4091a849f3f82b94202257b'
  - 't4_7b02d1dc4fd2_m50 rep2: frac_attack_ready=0.0, same md5, start_dist_a rep1=3.86 rep2=5.39'
runbook: null
related: [D0110, D0113, D0114]
---

## What was decided

**`t4_7b02d1dc4fd2_m50` is not evidence that the molecule engages Cys113.** Two
100 ns replicates from a byte-identical starting pose give **35.7%** and
**0.0%** warhead engagement. The second replicate's warhead never entered the
attack window once in 100 ns.

It must not be presented as the campaign's best result without both numbers.

## The measurement

Same `nac_v8_poses/t4_7b02d1dc4fd2.sdf`, `pose_rank 51`, `mode 50`;
`ligand_pose.sdf` md5-identical between the two work roots; same code path;
both analysed by `attack_sweep.geometry_stats` over a 500-frame series.

| | replicate 1 | replicate 2 |
|---|---|---|
| `frac_attack_ready` | **0.357** | **0.000** |
| `frac_attack_ready_angle` | 0.261 | 0.000 |
| median warhead distance | 3.65 A | 6.01 A |
| MINIMUM over 100 ns | — | 3.83 A |
| `start_dist_a` (production frame 1) | 3.86 A | 5.39 A |
| ligand RMSD mean / max | 7.96 / 9.68 A | 5.51 / 8.85 A |
| left the site | no | no |

## Why the wrong reading looked right

**Three separate disguises, all of which had to be seen through.**

**1. The replicate reports a different quantity under a similar name.** Replicate
2 ran through the elevation pipeline, which writes
`explicit_frac_frames_engaged = 0.9963` -- frames retaining >=25% of starting
protein contacts within 0.45 nm. That is RESIDENCE. Read as "99.6% engaged" it
makes replicate 2 look like a triumphant confirmation of replicate 1, when it is
the opposite. The same confusion is already recorded against the 2.2.0 record's
98.45%/99.80% figures for `t4_716800c125a7`. Catalogue disguise #1: two
populated, plausible columns, and the code names the wrong one.

**2. Both replicates held the pocket**, so every stability readout agrees and
only the warhead geometry disagrees. This is the `inert` tier doing its job --
held, warhead never engaged -- and it is the second time the campaign has caught
a comfortable non-reactive pose (D0114 was the first).

**3. D0113 was read as a statement about SHORT runs.** It measured an 18.9-point
median spread over 1.2 ns and attributed it to a small-sample statistic on a
fast-decorrelating quantity. The natural inference -- that 100 ns would average
the noise away -- is wrong. The spread here is 35.7 points over 100 ns, LARGER
than the short-run median, because the quantity that varies is not sampling
noise within a run: it is which basin the equilibration lands in, and a longer
production run samples that basin longer rather than escaping it.

## What it means

**The divergence is complete at production frame 1** (3.86 A vs 5.39 A). The
100 ns runs did not diverge; they inherited a divergence from a 300 ps
equilibration and then faithfully reported it. This is D0110's finding
(post-equilibration distance predicts engagement, rho = -0.519; docked distance
predicts nothing) seen from the other side: if equilibration sets the outcome
and equilibration is not reproducible, then neither is the outcome, at any
production length.

It also sharpens the chemist's observation, which @twu383 relayed on 2026-09-09:
a short MD moved a mediocre docked pose INTO a reactive one. That is the same
mechanism running favourably. m50 replicate 1 is one draw from it and replicate 2
is another. **A single equilibration is a sample, not a preparation.**

## What follows

- **n = 1 is not a result for this readout.** Any pose proposed on engagement
  needs replicate equilibrations, and the reported figure is the DISTRIBUTION,
  not the best draw.
- **The elevation gate is measuring a lottery.** `elevate_occupancy_min = 0.60`
  applied to a single run selects the favourable tail, which is precisely how
  the nine D0113 re-runs all came back lower.
- **Cheap and worth doing:** several independent equilibrations per pose, then
  score the *fraction of equilibrations* that reach attack geometry. That is a
  property of the molecule; occupancy within one lucky basin is not. 300 ps is
  ~2 min/mode, so 10 replicates cost less than one 100 ns run.
- Do NOT quote m50's 35.7% alone. Both numbers, or neither.
