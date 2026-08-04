---
id: D0057
title: The pocket-contact fit score does not recover crystal poses — worse than chance, in both receptor conditions
date: 2026-08-04
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/pose_vector.py
  - tests/test_pose_vector.py
evidence:
  - 'measured on D0046 harness: 82 cases, 1 crystal pose + 9 docked modes each, 34-residue 8A contact vector'
  - 'reference profile built LEAVE-ONE-OUT over the other 81 crystal poses — no leakage'
  - 'cross-docked into 6VAJ: crystal ranked #1 in 0/82 (0.0%), median rank 8.0, top-3 3.7%'
  - 'self-docked into own receptor: crystal ranked #1 in 5/82 (6.1%), median rank 7.0, top-3 22.0%'
  - 'chance with 10 candidates: 10% at rank 1, median rank 5.5'
  - 'the multimodality guard fired on EVERY case — a single median profile describes no mode'
---

# A geometric fit score, measured against ground truth, and it does not work

## What was built and why

Issue #14 proposes ranking on how well a pose sits in the pocket rather than on
the docking score — sound, because D0041 and D0046 have measured that the score
does not discriminate here. It asked specifically for an **agnostic** fit score:
one that does not require knowing in advance that the warhead belongs near
Cys113 (true for T_3, unavailable for T_1 and T_2).

`shared/pose_vector.py` implements one. For each of the 34 residues lining the
site, the minimum distance from any ligand heavy atom, clipped at 8 A — a
fixed-length vector regardless of molecule size, symmetry-invariant for free,
and agnostic by construction: the expected arrangement is **learned** from
poses that were experimentally observed, not declared.

## The test, and the result

Per #13's pre-registration rule, it was scored on D0046's harness before being
allowed to label anything: 82 cases, each a crystal pose plus the 9 docked
modes. The reference profile for each case was built **leave-one-out** over the
other 81 crystal poses, so nothing was fitted on the case being scored.

The question: can the fit score pick the true crystallographic pose out of a
lineup of 10?

| | crystal ranked #1 | median rank | crystal in top 3 |
|---|---|---|---|
| cross-docked into 6VAJ | **0 / 82 (0.0%)** | 8.0 | 3.7% |
| self-docked into own receptor | **5 / 82 (6.1%)** | 7.0 | 22.0% |
| chance | 10% | 5.5 | 30% |

**It is worse than chance in both conditions.** Not weak — inverted. The true
pose is systematically ranked BELOW the docked decoys.

## The receptor mismatch is real but is not the cause

Self-docking beats cross-docking by every measure (6.1% vs 0.0%, top-3 22.0%
vs 3.7%), which is a clean independent confirmation of the induced-fit
criticism in #14: a ligand's crystal pose transplanted into 6VAJ is not
expected to make 6VAJ-like contacts, because 6VAJ's pocket is shaped around
sulfopin. That is worth having.

But self-docking is still worse than chance. The receptor makes it worse; it
does not make it wrong.

## Two hypotheses, neither yet tested

**The reference profile is multimodal, and the code said so on every case.**
`reference_profile` warned for all 82. Pin1's site genuinely binds more than one
way, and a single median over several modes describes none of them — the target
sits in the trough between the modes, which is exactly the hazard
`mean_vector`'s docstring warns about, one level up. A per-mode profile is the
obvious next thing to try.

**More fundamental: the fit score and the docking score may not be
independent.** Vina's function is substantially shape complementarity and
contact terms, and docking *optimises* poses to make good contacts. A crystal
pose is optimised for nothing — it is simply true. So scoring by contact
profile may partly re-score what docking already did, and inherit the failure
D0046 measured rather than providing an orthogonal check. If that is what is
happening, no amount of tuning the vector fixes it, and the orthogonal signal
has to come from somewhere docking does not already look — stability under bias
(BPMD), or an explicit per-mode reference.

## The decision

**The fit score does not label, weight or rank anything.** It failed the
pre-registered test, on ground truth, and it is recorded as failed. Tuning it
until it passes on the same 82 cases would be choosing the answer, which is the
D0045 failure in a new costume.

`shared/pose_vector.py` stays, because two parts of it are correct and useful
independently of the score: the fixed-length contact vector, and
`representative()`, which returns an INDEX into real poses so that a "weighted
average pose" can never be a synthesised conformation (#14's own conclusion).

**What this saves.** #14 sequences BPMD before the fit score. BPMD as published
is ~10 replicas x 10 ns per pose; at 9 poses x 100 molecules that is ~37x every
nanosecond this project has ever simulated. Building it on top of a fit score
that is worse than chance would have spent that first and discovered this
after. Measured here for the cost of an afternoon and no GPU time.
