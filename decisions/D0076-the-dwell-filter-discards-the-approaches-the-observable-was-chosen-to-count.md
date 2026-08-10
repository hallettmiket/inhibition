---
id: D0076
title: The controls fail the sweep on a 100 ps dwell filter that discards exactly the brief approaches n_visits was chosen to count
date: 2026-08-10
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/attack_sweep.py
  - docs/prereg_attack_sweep.md
  - decisions/D0075-the-sweep-rejects-every-known-active.md
evidence:
  - 'rx_6VAJ (Sulfopin, crystal, cleaved) is inside the 2.8-4.2 A window in 87.4% of frames but attack-ready in 0.8% — distance is not what rejects it'
  - 'rx_6VAJ median S-C-LG angle 78.6 deg, rx_7F0M 69.4 deg, against SN2_ANGLE_MIN = 150 deg'
  - 'rx_7F0M made 13 RAW excursions into attack geometry; n_visits after filtering = 0'
  - 'rx_6VAJ 4 raw excursions -> 0; ref_Sulfopin__chloroacetamide 3 raw -> 0'
  - 'every control episode had median length 19.96-39.92 ps against MIN_DWELL_PS = 100.0'
  - 'frame_ps = 19.96, so 100 ps demands >= 5 CONSECUTIVE saved frames'
  - 'of 20 sn2_displacement candidates only 2 survive; survivors are 38 sn2_ring_opening + 18 michael_addition, whose angular criteria are laxer (#47)'
  - "rx_7F0M frac_attack_ready 0.0379 vs the 2 surviving SN2 candidates' median 0.038 — indistinguishable on the fraction; separated only by the dwell filter"
  - 'prereg_attack_sweep.md: "A covalent reaction needs ONE good approach, not sustained occupancy, so visits is the more mechanistically honest observable"'
---

# How a known covalent inhibitor scores zero visits

D0075 recorded that the sweep rejects every known active. @tt8804 asked the
obvious next question: *how is that possible?* It is three findings, and the last
one is a defect.

## 1. It is not the distance

`rx_6VAJ` — Sulfopin's crystal pose with the bond cleaved — sits inside the
2.8–4.2 Å near-attack window in **87.4% of frames**. The molecule is exactly where
it should be, almost all the time. It is scored attack-ready in **0.8%**.

## 2. It is the angle, and the SN2 bar is punishing for everyone

`SN2_ANGLE_MIN = 150°`. The controls' median S–C–LG angles are **78.6°**
(rx_6VAJ) and **69.4°** (rx_7F0M). Nowhere near.

This is not special pleading for the controls: of **20** `sn2_displacement`
candidates only **2** survive, against 58 survivors overall dominated by
`sn2_ring_opening` (38) and `michael_addition` (18), whose angular criteria are
laxer. #47 already flagged that cross-class comparison is biased; this quantifies
it. The SN2 gate rejects 90% of the molecules it judges.

## 3. The defect: the dwell filter discards what the observable exists to count

`rx_7F0M` entered attack geometry **13 separate times** and scored **zero
visits**.

`MIN_DWELL_PS = 100.0` requires each excursion to *last* 100 ps to be counted.
The sweep saves every 19.96 ps, so an episode must span **five consecutive saved
frames**. Every control episode had a median length of one to two frames, so all
of them — all 13 of rx_7F0M's — were filtered away.

The constant has a defensible engineering reason, recorded where it is defined:
counting in picoseconds rather than frames keeps the number comparable between
the 20 ps sweep and the 100 ps validation, since one molecule gave 57/26/14/7/3
visits at 100 ps/200 ps/500 ps/1 ns/2 ns.

**But it contradicts why the observable was chosen.** `prereg_attack_sweep.md`
says, in as many words:

> A covalent reaction needs **one** good approach, not sustained occupancy, so
> visits is the more mechanistically honest observable.

`n_visits` was adopted *because* reactivity does not require persistence — and it
is implemented with a persistence filter. A 20 ps approach is still an approach.
By the pre-registration's own reasoning rx_7F0M should read 13, not 0.

This is the recurring shape: the value taken (`n_visits` **after a dwell filter**)
is not the value the reasoning specified (independent excursions), and both are
populated, plausible, and named the same thing.

## What this does and does not license

- It does **not** rescue the controls' ranking. On `frac_attack_ready` rx_7F0M is
  0.0379 against the surviving SN2 candidates' median 0.038 — indistinguishable.
  The dwell filter is what separated them, not merit.
- It does **not** mean 100 ps is wrong for every purpose. Occupancy questions want
  a dwell floor. Reaction-competence questions, by the prereg's argument, do not.
- It does **not** explain the angle. 70–79° against a 150° bar is a separate
  finding, and the more serious one: an experimentally determined pose, with the
  leaving group constructed anti to the sulfur, does not satisfy our SN2 criterion
  during dynamics.

## What follows

1. **Re-derive `n_visits` with the dwell floor at the save interval** (one frame),
   and report both: `n_visits_raw` and the dwelled count. The raw number already
   exists in every row written so far, so this costs no GPU.
2. **Re-read D0075 against the re-derived number.** If the controls separate from
   the candidates on raw visits, the sweep may be salvageable as a gate; if they
   do not, D0075 stands unchanged and the sweep cannot gate.
3. **The angle question is separate and open.** Why a crystal-derived pre-reaction
   geometry samples 70–79° rather than approaching 180° is not answered here, and
   it bears on whether the SN2 criterion is measurable in unrestrained MD at all.
