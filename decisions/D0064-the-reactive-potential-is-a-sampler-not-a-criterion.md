---
id: D0064
title: The reactive pair potential is a sampler, not a criterion — the approach angle must be scored separately
date: 2026-08-05
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/nac_criterion.py
  - scripts/nac_screen.py
  - docs/ranking_rationale.md
evidence:
  - 'reactive docking, published parameterisation (r_eq 1.8 A, eps 2.5, neighbour radii scaled 0.5x), one T_4 chloroacetamide, 3IKD, 20 runs: C1..SG median 1.59 A, 20/20 within 2.5 A'
  - 'the same 20 poses, S-C-Cl angle: median 97.6 deg, max 128.9 deg, 0/20 above 150 deg, 8/20 below 90 deg (leaving group between the sulfur and the carbon)'
  - '40-run repeat of the published parameterisation: angle median 78.8 deg, max 117.2 deg, 0/40 viable'
  - 'NAC parameterisation (r_eq 3.2 A, eps 1.0, neighbour radii unscaled), 40 runs: S..C median 3.39 A, angle median 86.6 deg, max 151.3 deg, 3/40 above 150 deg'
  - 'meeko get_reactive_config: the modified pair term is intnbp_r_eps, a 13-7 function of DISTANCE only — no angular dependence exists in the functional form'
  - 'meeko reactive.py applies r13_scaling / r14_scaling (default 0.5) to atoms 1-3 and 1-4 from the reactive centre, deflating the leaving group''s steric radius'
---

# The reactive potential is a sampler, not a criterion

## What was measured

Reactive docking (D0063) does exactly what it claims, and the claim is narrower
than it first appears. On one T_4 chloroacetamide against 3IKD, **20 of 20 poses
placed the electrophilic carbon 1.55 A from Cys113's sulfur** against a 1.80 A
target. Read as "does the warhead reach the nucleophile", that is a total
success, and against free docking — where a correct pose is *findable* 41.5% of
the time and the score picks one at chance (D0061) — it looks like the problem
solved.

**Every one of those poses was chemically dead.** SN2 at an sp3 carbon requires
the nucleophile *anti* to the leaving group, near 180 deg. The measured S-C-Cl
angle had **median 97.6 deg and maximum 128.9 deg**, and in 8 of 20 the chlorine
sat physically between the sulfur and the carbon it was meant to be attacking. A
40-run repeat gave 0/40 above 150 deg.

## Why, and why no tuning fixes it

The modified term is `intnbp_r_eps`, a 13-7 potential in the interatomic
**distance**. It is isotropic by construction: it rewards the reactive atoms
being close and is indifferent to the direction of approach. **A distance
restraint cannot encode an angle.** This is a property of the functional form,
not of the parameter values, so no choice of `eps_12` or `r_eq_12` changes it.

A second effect compounds it. The published protocol scales 1-3 and 1-4
neighbour radii by 0.5, letting atoms adjacent to the reactive centre
interpenetrate. That is what makes a 1.55 A approach geometrically possible at
all — and it is also what lets the leaving group occupy the attack vector for
free.

## What was tried instead, and how far it got

If deflated sterics let the leaving group sit in the way, real sterics at a real
near-attack distance might push it out — making the angle self-enforcing rather
than needing a separate criterion. Re-parameterised to **3.2 A with full
neighbour radii**, viable backside poses appeared for the first time: **3/40
above 150 deg, versus 0/40**, with a maximum of 151.3 deg.

Better, and not sufficient. The median was still 86.6 deg. The isotropic well
rewards every approach direction equally, so the search settles on whichever is
sterically cheapest, and side-on is always cheapest.

## The decision

**Reactive docking is adopted as a SAMPLER and not as a criterion.** It puts
poses in the right region of the pocket, cheaply — 40 runs in ~2 seconds — and
that is all it is asked to do. Whether a pose is reaction-competent is a separate
question, answered afterwards by `shared/nac_criterion.py` on the poses it
generates.

Two parameter changes follow from the measurement and are adopted with it:

- **`r_eq_12` = 3.2 A, not 1.8 A.** A near-attack conformation is a van der
  Waals contact — the reactant state, bond not yet formed. 1.8 A is a *covalent
  bond* distance, which is past the transition state and a geometry the free
  molecule cannot occupy while still carrying its leaving group. Docking the
  free form to a bonded distance asks it to be somewhere it cannot be.
- **Neighbour radii unscaled.** We want the leaving group's real steric bulk,
  since part of the question is whether it obstructs the approach.

## Consequences

- `docs/ranking_rationale.md` stage 3 stops being optional. It was written as a
  screen applied to poses; it is now the only thing standing between the
  pipeline and a ranking built entirely on chemically impossible geometry.
- The criterion is **mechanism-specific**, and not as a refinement. SNAr proceeds
  by attack along the ring normal with the leaving group staying in-plane, so its
  S-C-LG angle is near **90 deg** — the exact value that is dead for SN2. One
  distance-plus-angle rule applied across mechanisms inverts the verdict for one
  of them.
- The windows are **pre-registered** from textbook stereoelectronics before any
  candidate is scored (D0045), and every result carries the raw measured angle so
  a window can be re-drawn without re-docking.

## What this does not settle

Whether the resulting score **discriminates between molecules** at all. That is
the rationale's second stated failure mode and it is measured separately by
`scripts/nac_screen.py`; a criterion that everything passes ranks nothing.
