---
id: D0075
title: The 10 ns sweep rejects every molecule known to react with Cys113 — Sulfopin included — while passing a quarter of generated candidates
date: 2026-08-10
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/attack_sweep.py
  - scripts/nac_screen_v2.py
  - shared/nac_criterion.py
  - docs/prereg_attack_sweep.md
evidence:
  - 'Sulfopin through the production screen (500 runs, pose splitting): 1 mode, 465 poses in it, 47 reaction-competent — it is NOT rejected at the geometry stage'
  - 'Sulfopin, docked, 10 ns sweep: frac_attack_ready 0.0100, n_visits 0, frac_in_window 0.2236, min_dist 3.22 A'
  - 'Sulfopin, crystal pose with the covalent bond cleaved (rx_6VAJ), 10 ns sweep: frac_attack_ready 0.0080, n_visits 0'
  - 'Liu-2022-ZL-Pin13, crystal pose cleaved (rx_7F0M), 10 ns sweep: frac_attack_ready 0.0379, n_visits 0'
  - 'attack_sweep.py own survivor rule: "A survivor is a molecule with a SUSTAINED episode (n_visits > 0)"'
  - 'of 233 swept candidates, 58 (25%) have n_visits > 0'
  - 'Sulfopin ranks 104 of 234 on attack-ready; the crystal form 106; Liu-2022 76'
---

# The positive control fails our own screen

@tt8804 asked the question directly before the Monday review: *would we have
caught Sulfopin?*

**No.**

## What the screen does and does not reject

Sulfopin is not thrown out early. Through the production protocol — 500 runs,
reactive docking on 3IKD, pose splitting — it yields **one mode, 465 poses, 47 of
them reaction-competent**. The geometry stage finds it perfectly acceptable.

It fails at the **10 ns sweep**, which is the cut that decides what earns a 100 ns
run:

| molecule | provenance | attack-ready | sustained visits | survivor? |
|---|---|---:|---:|---|
| Sulfopin | docked | 0.0100 | **0** | no |
| Sulfopin | crystal pose, bond cleaved | 0.0080 | **0** | no |
| Liu-2022-ZL-Pin13 | crystal pose, bond cleaved | 0.0379 | **0** | no |

`attack_sweep.py` states its own rule in its output: *a survivor is a molecule
with a sustained episode (n_visits > 0); a single 20 ps touch is not evidence of
anything.* By that rule **all three known actives are rejected**, and none would
have been elevated to 100 ns.

Meanwhile **58 of 233 candidates (25%)** do have a sustained visit. The criterion
passes a quarter of our generated matter and none of the chemistry that is known
to work.

## Why this is not a docking artefact

This is the fork #47 left open: either docking mislocates the known actives so the
criterion never sees their real geometry, or the criterion itself rejects real
geometry.

**The two Sulfopin measurements agree.** Docked from SMILES it reads 0.0100; taken
from the crystal structure with the covalent bond cleaved and the leaving group
rebuilt — every other atom at its crystallographic coordinate — it reads 0.0080.
Same answer from a pose that was *determined experimentally*.

Docking is not what is failing here. The criterion rejects the true geometry.

## What it does not establish

- **Not a claim that the criterion is meaningless.** It is a claim that, as a
  gate, it does not admit the molecules it must admit.
- **n = 3 actives.** Small, and two of them are the same molecule by two routes.
- The reactant forms are **rebuilt**, not observed: the leaving group is placed.
  That is one modelled degree of freedom, though the docked/crystal agreement
  argues it is not the cause.
- Whether these molecules would have been recovered *later* by 100 ns MD is a
  different question — they never got that far, which is the point.

## What follows

1. **The sweep cannot be used as a gate in its current form.** Anything it filters
   is filtered by a criterion that rejects Sulfopin. Ordering may survive;
   *cutting* does not.
2. **Re-read every shortlist that the sweep produced** with this in mind. The 58
   survivors are 58 molecules that scored better than the incumbent on a reading
   that puts the incumbent at rank 104.
3. **The 100 ns runs on the controls become the arbiter**, not the confirmation.
   If Sulfopin and Liu engage well at 100 ns having failed the 10 ns sweep, the
   sweep is measuring the wrong thing and the funnel is upside down. Those runs
   are in flight.
4. This is the same shape as D0041 and D0046 one stage further on: a scorer that
   is populated, plausible, and does not separate what it is supposed to separate.
