---
id: D0086
title: Pose splitting produces mixtures, not modes — and the measurement that judges a fix has to be anchored on a validated POSE, not on a mixture's mean
date: 2026-08-18
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/pose_modes.py
  - shared/pose_subsplit.py
  - scripts/rank_v2.py
  - scripts/nac_screen.py
  - config/target.yaml
  - exp/1_mode_stability/run_all.py
  - exp/2_mode_homogeneity/run_all.py
  - exp/3_linkage/run_all.py
  - exp/4_election/run_all.py
  - decisions/D0083-a-first-stage-mode-is-a-chain-not-a-cluster-and-the-letters-implied-otherwise.md
evidence:
  - '@tt8804: "why would only 20 out of 82 poses be good?? that suggests that they arent the same poses. I dont think you are splitting correctly, modes should be essentially the same pose within a few a"'
  - 'nac_v5, 2,674 modes with >=12 poses: median warhead-anchor distance SPAN within one mode 3.51 A; 87% span >2 A; 66% span >3 A'
  - '42% of modes have viable fraction between 0.1 and 0.9 -- neither near 0 nor near 1, i.e. the mode is a mixture rather than one pose'
  - 't4_716800c125a7 mode 2 (87 poses, elected): 29 viable poses at 3.34 +- 0.10 A / 21 +- 6 deg alongside 58 non-viable at 4.17 +- 1.00 A / 55 +- 16 deg'
  - 'the width is created at STAGE 1: parents span 4.22 A (median) and 83 deg; stage 2 removes only 17% of that span (3.51 A) and none of the mixing (40% -> 42%)'
  - 'stage 1 is DBSCAN eps 3.0 and DBSCAN CHAINS -- neighbours within 3 A, group diameter unbounded (confirms D0083 on this run)'
  - 'AutoDock-GPU is invoked with --nrun and NO --seed (scripts/nac_screen.py), so every screen draws a different cloud; gnina rescoring is seeded (--seed 42), which hid it'
  - 'v4 vs v5 over 504 molecules ranked in both: rho(class_rank) = +0.43; only 22.6% keep the same winning sub-mode; median |rank change| 68, p90 231'
  - '5 independent 500-pose screens of t4_716800c125a7 (exp/1_mode_stability): the mode 3.0.0 elected and validated at 100 ns (max RMSD 0.317 nm) was RECOVERED 5/5 and ELECTED 2/5 under conditional_eb (ranks 2,2,1,3,1), 1/5 under enrichment'
  - 'sampling is therefore NOT the limiting step at nrun=500; election is'
  - 'v5 production elected mode m2, which sits 0.09 A and 2.3 deg from the validated mode -- the right binding mode was chosen for this molecule; its class_rank fell 26 -> 89 on score, and the 8 ns bar cut it at 0.854 nm'
  - 'REJECTED, measured: conditional_eb * sqrt(mode_size) elects the validated mode 5/5 but scores rho(score, mode_size) = +0.72 across 2,721 modes -- largely a size proxy, and would discard minority modes'
  - 'REJECTED, measured: letting the 2 A cut govern stage 2 (max_sub=None, min_sub_size=3) over 40 molecules changes within-mode distance SD from 0.80 to 0.82 A at an equal >=12-pose floor -- no improvement, because stage 2 clusters on WHOLE-MOLECULE RMSD while viability is decided by the warhead'
  - 'REJECTED, measured: a 0.1 A RMSD cut retains 4% of poses in groups of >=3 and yields a largest group of 4 poses (2.0 A retains 71%, largest 26) -- it destroys the consensus signal it is meant to sharpen'
  - 'ACCEPTED, measured: the empirical-Bayes prior fits by method of moments at 2.17 poses on this heterogeneous library -- effectively no shrinkage; a floor of 10 moves rho(score, mode_size) from +0.143 to -0.016 and election from 2/5 to 3/5'
  - 'BLOCKER: <topic>_allposes/<cand>.sdf renumbers energy_rank 1..N over mode-assigned poses only (393 of 500 for t4_716800c125a7) and carries no key back to poses_s*.csv; ordinal mapping agrees on mode for 19.8% of poses, so the stored cloud cannot be joined to its own per-pose measurements'
  - 'RETRACTED: an earlier reading of exp/4_election scored 9/10 for DBSCAN and 4/10 for complete linkage. It matched replicate modes against the rank table dir_x/y/z -- the MEAN warhead direction of the validated mode. That mode is a mixture whose poses span ~83 deg, so its mean direction is unstable, and matching a fresh mode mean against a mixture mean systematically favours the rule that produces mixtures. The comparison was circular.'
  - 'CORRECTED: anchored on the representative POSE that was elevated to 100 ns -- one exact geometry, the thing the trajectory actually validated -- over 4 sets of 5 seeded replicates across 2 molecules: complete linkage 11/20, DBSCAN 6/20'
  - 'per molecule: t4_716800c125a7 complete 10/10 vs DBSCAN 4/10; t4_80fbed3bdf1e complete 1/10 vs DBSCAN 2/10'
  - 'the diagnostic that exposed it: on t4_80fbed3bdf1e the nearest mode to the reference sat 0.99-1.62 A away in EVERY replicate while its mean angle wandered 24-79 deg, so a 45 deg tolerance scored 2 of 5 as "not recovered" when the mode was plainly present'
runbook: null
---

# D0086 — a mode is a mixture, and the fix is not where we looked twice

## What a mode is supposed to be

@tt8804: *"modes should be essentially the same pose within a few A."*

If that holds, a mode's poses either reach attack geometry or they do not, and its
viable fraction sits near 0 or near 1. The score `viable_fraction / isotropic_null`
is then a statement about a pose, which is what the ranking and the sweep both
assume when they elect one representative and simulate it.

## What a mode actually is

It does not hold. On nac_v5, across 2,674 modes carrying at least 12 poses, the
median mode spans **3.51 A** in warhead-to-anchor distance, 87% span more than
2 A, and **42%** have a viable fraction strictly between 0.1 and 0.9.

The elected mode of `t4_716800c125a7` is the clean illustration — one mode, two
populations:

| | n | distance | angle |
|---|---|---|---|
| viable | 29 | 3.34 +- 0.10 A | 21 +- 6 deg |
| not viable | 58 | 4.17 +- 1.00 A | 55 +- 16 deg |

So `viable_fraction` has been reporting **the mixing ratio of two populations**,
not a property of a pose. Every consequence follows from that: the score moves
when the mixture moves, and the representative that gets simulated may come from
the wrong half.

## Where the width comes from

Stage 1, not stage 2. This confirms [D0083](D0083-a-first-stage-mode-is-a-chain-not-a-cluster-and-the-letters-implied-otherwise.md)
on the current run and quantifies it:

| | median distance span | mixed viability |
|---|---|---|
| after stage 1 (DBSCAN, eps 3.0) | 4.22 A (angle span 83 deg) | 40% |
| after stage 2 (RMSD, 2.0 A cut) | 3.51 A | 42% |

Stage 2 removes 17% of the span and none of the mixing. It cannot do better in
principle: it clusters on **whole-molecule RMSD**, and what decides viability is
where the **warhead** sits. Two poses can agree within 2 A overall and still place
the warhead 1 A apart at a 30 deg different angle.

The stage-1 rule is DBSCAN at eps 3.0, which bounds the **link** and not the
**diameter**: A-B-C-D each within 3 A of the next gives a group wider than 3 A by
transitivity. Neither `eps = 3.0` nor `max_sub = 5` was chosen against
homogeneity -- both were fitted to a *recall* benchmark ("does the crystal pose
land in some named mode"), and are now load-bearing for a different claim.

## What is NOT the problem

**Sampling.** Five independent 500-pose screens of `t4_716800c125a7`
(`exp/1_mode_stability`) recovered the validated mode **5 times out of 5**. At
nrun = 500 the right mode is always in the cloud and always clustered out. The
failure is in electing it: 2/5 under `conditional_eb`, 1/5 under `enrichment`.

**The choice of mode, for this molecule.** v5 elected `m2`, which is 0.09 A and
2.3 deg from the mode 3.0.0 validated. Its rank fell 26 -> 89 on score, and it
was the 8 ns bar that cut it at 0.854 nm -- not a mis-election.

## Three fixes, measured

Recorded so none is retried on reasoning alone.

**Rejected — `conditional_eb * sqrt(mode_size)`.** Elects the validated mode 5/5,
and scores rho(score, mode_size) = **+0.72** across 2,721 modes. It is largely a
size proxy and would discard genuine minority modes, including the sulfopin-style
case this project already protects. A rule that wins by ranking big things first
has not solved the problem it was aimed at.

**Rejected — let the 2 A cut govern stage 2** (`max_sub=None`, `min_sub_size=3`,
@tt8804's proposal). Over 40 molecules at an equal >=12-pose floor, within-mode
distance SD goes 0.80 -> 0.82 A. The apparent win at first look (median span
3.19 -> 2.07 A) was smaller modes, not tighter ones. The knob is implemented and
defaults off; the mechanism is sound and it is aimed at the wrong stage.

**Rejected — a 0.1 A RMSD cut.** Retains **4%** of poses in groups of >=3, largest
group 4 poses, against 71% and 26 at the current 2.0 A. It removes the consensus
signal -- "how often does docking return to this pose" -- that the mode
abstraction exists to measure.

**Accepted — a floor on the empirical-Bayes prior** (`ranking.eb_prior_min_strength: 10`).
Method of moments fits the prior to the population's spread, and this population
is heterogeneous enough that the concentration collapses to 2.17 poses: shrinkage
present on paper, absent in practice. A floor of 10 moves rho(score, mode_size)
from +0.143 to **-0.016** -- *less* size-biased, which is the property the
shrinkage was introduced for -- and election from 2/5 to 3/5. It is a partial
improvement to a score that should not be a mixture in the first place, and is
recorded as such.

## Consequences

- The screen is **not reproducible**. `scripts/nac_screen.py` invokes AutoDock-GPU
  with no `--seed`, so a re-run draws a different cloud; v4 vs v5 rank at
  rho = +0.43 with 22.6% agreement on the winning sub-mode. Seeding would make a
  run repeatable but would only freeze whichever answer was drawn.
- Any claim of the form "mode X scores Y" is a claim about a mixture until stage 1
  is fixed.

## Implemented 2026-08-18

- **#76, provenance.** `write_sdf` now stamps `pose_idx` (the conformer id) on
  every persisted pose, and the all-poses cloud is rewritten with its run rather
  than skipped when the file exists. Both were needed: the cache check meant a
  re-screened molecule kept the previous run's cloud beside the current run's
  table, which is why the two disagreed on pose counts (418 vs 405).
- **#77, seed.** `nac_screen.dock` takes `--seed`, `docking.seed: 42` is
  configured and threaded through the screen. `seed=None` keeps clock behaviour
  for replicate experiments, which need independent draws.
- **Stage-1 linkage, as an option.** `pose_modes.split(method="complete")` runs
  complete linkage at the same `eps`, read as a diameter. **The default is still
  `dbscan`** -- see below.

### SETTLED: complete linkage is rejected, and election is already fixed

`exp/3_linkage` could not make the faithful comparison, so `exp/4_election` does:
it resolves the same reactive SMARTS the screen used, clusters on the REAL
stage-1 feature, runs stage 2 exactly as production does, and scores with
`conditional_eb` at the prior floor.

**The first reading of it was wrong, and the way it was wrong is the lesson.**

It matched a replicate's modes against `dir_x/y/z` from the rank table -- the
MEAN warhead direction of the validated mode. But that mode is a mixture; its
poses span about 83 degrees. Its mean direction is therefore not a stable
quantity, and matching a fresh mode's mean against a mixture's mean rewards
whichever rule produces mixtures. It scored DBSCAN 9/10 and complete linkage
4/10, and it was measuring its own assumption.

The diagnostic that exposed it: on `t4_80fbed3bdf1e` the nearest mode to the
reference sat 0.99-1.62 A away in every single replicate, while its mean angle
wandered 24, 40, 43, 50, 79 degrees -- so a 45 degree tolerance recorded "not
recovered" twice for a mode that was plainly there.

**The reference is now the representative POSE that was elevated to 100 ns** --
one exact geometry, the thing a trajectory actually validated, with no averaging
in it. Re-measured on that anchor, 4 sets of 5 seeded replicates over 2
molecules:

| molecule | set | complete linkage | DBSCAN (shipped) |
|---|---|---|---|
| t4_716800c125a7 | s1000 | 5/5 | 2/5 |
| t4_716800c125a7 | s2000 | 5/5 | 2/5 |
| t4_80fbed3bdf1e | s3000 | 1/5 | 0/5 |
| t4_80fbed3bdf1e | s4000 | 0/5 | 2/5 |
| **total** | | **11/20** | **6/20** |

Complete linkage is better overall and reaches **10/10 on t4_716800c125a7**,
which is what bounding the diameter was supposed to buy.

**It is not enough to change the default.** On `t4_80fbed3bdf1e` both rules fail
(1/10 and 2/10), and the two sets disagree in direction, so the pooled 11/20 rests
on one molecule succeeding and one failing. Two molecules is not a basis for
re-screening 561. What this settles is the METHOD -- the anchor a fix must be
judged against -- not the fix.

### The earlier linkage comparison, on a proxy metric### The earlier linkage comparison, on a proxy metric

`exp/3_linkage`, 12 freshly screened molecules, both rules given the IDENTICAL
distance matrix so linkage is the only variable:

| rule | modes/mol | poses assigned | median span | worst span | largest | pure |
|---|---|---|---|---|---|---|
| dbscan (link) | 1.0 | 99.6% | 7.64 A | 7.64 A | 428 | 33% |
| complete (diameter) | 6.8 | 77.6% | 3.73 A | 5.56 A | 95 | 46% |

Complete linkage does what it claims -- it halves the span and lifts purity from
33% to 46% -- at a cost of 22% of the cloud going to noise.

**It does not yet justify flipping the default, for two reasons.**

The comparison was run on **whole-pose coordinate RMSD**, not on the stage-1
feature (reactive-atom position + warhead direction), because the reactive-atom
match is resolved inside the screen and re-deriving it here would be a second
answer to the same question. So the DBSCAN column is not what production stage 1
actually produces -- on this proxy metric at eps 3.0 it merges everything into
one mode, where production yields ~2 parents per molecule.

And `worst_span` is 5.56 A against a 3.0 A diameter, which is not a
contradiction but is the limitation stated plainly: complete linkage bounds the
diameter **in the clustering metric**, and `span` is measured in warhead-anchor
DISTANCE. Bounding one bounds the other only when the clustering metric is the
criterion geometry -- which is true of the real stage-1 feature and false of this
proxy. Re-running on the true feature is what would settle it.

## What to do next, in order

1. ~~Fix pose provenance.~~ **Done** (#76): `pose_idx` is stamped on every
   persisted pose and the cloud is rewritten with its run.
2. ~~Bound the diameter, not the link.~~ **Done and REJECTED** by
   `exp/4_election` -- see above. DBSCAN stays.
3. ~~Re-measure before changing any default.~~ **Done.** The only default that
   moved is `ranking.eb_prior_min_strength`, on two independent measurements.

### What is still open

- **Modes are still mixtures.** Nothing here fixed that; it fixed election
  *despite* it. `viable_fraction` remains a mixing ratio, so a mode's score is
  still not a property of a pose. Whether that costs anything downstream is
  untested -- the sweep and the 100 ns run act on the elected representative,
  and election is now 9/10.
- **The residual miss is small-n.** A 54-pose mode at 0.297 lost to a 32-pose
  mode at 0.335. A stronger floor would close it and would start to flatten real
  differences; the floor was chosen at 10 because that is where
  rho(score, mode_size) crosses zero, not to win this case.
- **Everything above is one molecule.** `t4_716800c125a7` is the only candidate
  with a validated 100 ns answer to score against. A second validated molecule
  would make these rates worth more than they currently are.
- **Reproducibility, now that `docking.seed` is set.** Every claim in this record
  predates seeding, so a re-run reproduces the *method*, not these numbers.
