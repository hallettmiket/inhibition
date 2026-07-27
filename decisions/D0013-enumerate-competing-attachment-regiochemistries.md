---
id: D0013
title: Enumerate competing attachment regiochemistries rather than choosing one
date: 2026-07-27
status: accepted
approach: t4
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/approaches/t4_combinatorial.yaml
  - data/reference/warhead_classes_3.csv
  - shared/warhead_library.py
evidence:
  - 'BDHI: C3 bears the Br and is the Cys attack site, so attachment is limited to C4 or C5'
  - '1,4-naphthoquinone: C2/C3 are the Michael acceptor positions, so attachment is C2 or the benzo ring'
  - 'all four candidates parse and match their mechanism SMARTS'
  - 'the discriminating evidence (5b validity, docking geometry, LUMO) comes from steps already in the pipeline'
runbook: null
---

## Context

Two warhead classes had verified chemotypes but no established attachment
regiochemistry, and the honest answer from the PI was "I do not know". That is a
real state, not an oversight — the literature does not settle it for a sulfolane
core, and guessing would put an arbitrary choice underneath every T_4 result.

The constraint is tighter than it first looks: the reactive atom must stay free.
That reduces an open design question to four concrete candidates.

## Decision

Enumerate **all four** as separate warhead classes and let the gates decide.
T_4 opts into `DESIGNED_UNTESTED` alongside `VERIFIED`.

The discriminating evidence costs nothing extra, because it comes from steps
already in the pipeline:

- **step 5b** — is the attached warhead still a genuine electrophile of its
  class, or did coupling kill it? A regiochemistry that fails this is refuted,
  not merely disfavoured.
- **step 6** — can Cys113 SG reach the reactive atom with the core in the way?
  A blocked approach shows up directly as poor covalent docking geometry.
- **step 7** — is its LUMO inside the window bounded by real actives?

Competing regiochemistries are reported **separately** through to ranking. Which
attachment works is a finding worth carrying, not an implementation detail to
collapse.

## Consequences

Eight enumerable classes instead of four, so the enumerated library doubles.
That is affordable: the alert gate and 5b run before covalent docking, which is
the throughput wall, and MM-GBSA is still capped at per-class representatives.

DESIGNED_UNTESTED classes remain barred from anchoring the reactivity window
(control B5) however well they dock — surviving a gate is not the same as being
a validated Pin1 active.

If both members of a pair survive every gate, that is also an answer: the
attachment is not the discriminating variable, and both should carry forward.
