---
id: D0082
title: Sulfopin's crystal pose is both generated and carried; the near-attack criterion is what ranks it last
date: 2026-08-12
status: proposed
approach: shared
decided_by: '@twu383'
origin: adversary
supersedes: []
superseded_by: null
affects:
  - scripts/crystal_pose_audit.py
  - tests/test_crystal_pose_audit.py
  - config/target.yaml
  - shared/nac_criterion.py
  - decisions/D0075-the-sweep-rejects-every-known-active.md
  - decisions/D0077-the-crystal-reactant-control-models-an-adduct-as-a-michaelis-complex.md
  - decisions/D0081-only-T4-is-ranked-and-T3-was-why-acrylamide-looked-dominant.md
evidence:
  - 'class_rank 497 of 540 ranked chloroacetamide modes (bottom 8.1%), on the T_4-only ranking of D0081'
  - '455 persisted poses; 99 (21.8%) within 2.0 A of the 6VAJ adduct, best 1.274 A'
  - 'reproduces the 2.2.0 figure #64 said could not be checked: 95/456 = 20.8%, best 1.43 A'
  - 'best exported representative 1.723 A (mode 3 / label 0d) -- a successful recovery by the 2 A criterion, and it is the TOP-RANKED mode'
  - 'of the 99 near-crystal poses, 93 are inside the 2.8-4.2 A window; 85 of those fail on the SN2 angle alone'
  - 'near-crystal median S-C-Cl angle 117.1 deg against a 150 deg bar; whole-cloud median 103.8 deg'
  - 'enrichment by shell: all poses 0.591 | within 2.0 A 1.206 | within 1.5 A 2.764; budget_floor is 4.0'
  - 'only 21 of 455 poses (4.6%) clear 150 deg, against an isotropic null of 6.7%'
  - 'crystal control frame verified: reactive C to Cys113 SG = 1.979 A after superposition (fit 0.524 A over 113 CA)'
  - 'join to the screen rows validated: per-mode counts exact, r = 0.9946, flexible-vs-rigid SG shift 0.351 A'
runbook: null
---

## Context

Sulfopin — the known covalent Pin1 inhibitor, carried through the screen as an
ordinary candidate so the pipeline has to place it without knowing what it is —
ranks in the **bottom 7% of its own warhead class** on the 3.0.0 screen, best
mode scoring `enrichment` 0.71, i.e. below 1.0 and so worse than random
orientation. #64 posed the question as generation versus selection, with
evidence on both sides, and recorded a blocker: no sulfopin crystal reference
existed in the receptor's coordinate frame, so the question could not be
measured.

**The blocker was wrong, and the way it was wrong is the catalogue's shape.**
The reference has existed since 2026-08-09. `crystal_controls.py` superposes
each covalent Pin1 crystal onto the production 3IKD and writes the ligand in
that frame; `crystal_control_poses/xtal_6VAJ.sdf` is sulfopin, 16 atoms, fit
0.524 Å over 113 Cα, reactive carbon 1.979 Å from Cys113 SG. What #64 found
instead was `m2_covalent_smoke/sulfopin.sdf`, an **origin-framed docking
input** — a file with the right molecule, the right name and the wrong frame,
which returned best RMSD 12.07 Å and 0 of 455 poses within 2.5 Å. That reads as
a catastrophic docking failure and is a coordinate-system mismatch.

## Decision

Neither reading in #64 is correct. **The pose is generated, the pose is
carried, and the criterion rejects it.**

1. **Generation is not the problem.** 99 of the 455 persisted poses (21.8%) are
   within 2 Å of the crystal adduct; the best is 1.274 Å. This independently
   reproduces the 2.2.0 measurement #64 said could not be reproduced or refuted
   (95 of 456, best 1.43 Å) on fresh 3.0.0 data — same receptor, same run as the
   scores, pH 7.4 species.

2. **Pose selection is not the problem either.** The best exported
   representative is **1.723 Å**, a successful recovery by the Astex 2 Å
   criterion — and it is the representative of mode 3 (label `0d`), the
   molecule's **top-ranked** mode. The pipeline carried a crystallographic pose
   and ranked it first among sulfopin's four modes.

3. **The criterion is what puts it last.** Of the 99 poses that reproduce the
   crystal, 93 are inside the 2.8–4.2 Å distance window and only 8 are viable:
   **85 fail on the SN2 angle alone**, at a median 117.1° against a 150° bar.
   Across the whole cloud only 4.6% of poses clear 150°, against an isotropic
   null of 6.7% — which is why `enrichment` lands below 1.

**And the decisive number: restricting to poses that ARE the answer does not
rescue it.** Enrichment computed over successively tighter shells around the
crystal pose rises monotonically — 0.591 over all poses, 1.206 within 2 Å,
2.764 within 1.5 Å — and never reaches the sweep's `budget_floor` of 4.0. The
criterion is positively associated with being at the crystallographic pose, so
it is not measuring noise; it is simply too strict to pass the one molecule
known to work, **at any achievable quality of docking**.

So no improvement to pose generation and no improvement to pose selection can
put sulfopin through the current gate. That is a property of the criterion, not
of the molecule or of the search.

## Why it looked right

Three things made "the docking is bad" the natural reading, and each of them is
true in a way that does not support the conclusion.

* **The rank is genuinely terrible.** 496 of 540 ranked chloroacetamide modes
  beat sulfopin — bottom 8.1%, and unchanged by D0081's move to a T_4-only
  ranking, which was expected since chloroacetamide is a T_4 class. Nothing
  about the number is wrong — `enrichment` is computed
  exactly as defined. A metric can be correctly computed, correctly ordered, and
  still not be the quantity you need.
* **The frame mismatch produced a number, not an error.** 12.07 Å and 0 of 455
  is precisely what a failed docking looks like. Nothing raised; the file parsed;
  the molecule was the right molecule.
* **The whole-molecule endpoint hides it.** D0062 already measured that whole-
  molecule RMSD and reactive-region placement correlate at only ρ = +0.433. Here
  the reverse bites: the poses are *right* by RMSD and rejected on an angle, so
  a recovery rate alone would have said "docking works" and stopped, without
  reaching the criterion that does the damage.

The deeper reason is stated in D0077 and holds again here: **the deposited
crystal pose is the covalent adduct, a post-reaction state.** Its reactive
carbon is bonded to SG at 1.78 Å, below the near-attack window's own floor, so
the crystal cannot be scored by the criterion at all. A docked pose that
reproduces it inherits the geometry of a formed bond rather than of an approach
to one — and the 150° anti-periplanar arrangement is reached by motion, which a
rigid docked snapshot is not obliged to show. D0077 measured exactly that on the
crystal reactant: it enters production at 100.9° and spends 87.4% of a 10 ns run
inside the distance window at a median 78.6°.

## Consequences

* **The sweep's `budget_floor` of 4.0 excludes the positive control.** The
  259-mode sweep now running is a defensible experiment — the sweep is what
  tests the criterion — but its shortlist cannot be described as validated,
  and `sweep_rule.floor` remains `null` for the right reason. The stratified
  pilot is now the load-bearing next step, and it should carry sulfopin as a
  known-positive stratum rather than being blind to it.
* **A floor set on `enrichment` cannot be calibrated on this target until the
  criterion admits a known active.** Two routes: relax the angle bar and measure
  what it costs in specificity, or accept that the ground-state pose is the
  wrong observable and let the sweep — not the docked snapshot — supply the
  angle. The second is what D0076 and D0077 were already pointing at.
* **`crystal_controls.py`'s output is the reference of record for pose
  recovery.** `scripts/crystal_pose_audit.py` refuses any structure whose
  reactive atom is not a bond length from the anchor, so the frame mistake in
  #64 raises rather than reporting a recovery rate.
* **The criterion cannot be recomputed from the persisted pose cloud.** The
  screen measures the approach to the **flexible** Cys113 sulfur, taken from the
  docked model (`nac_screen.sg_position`), and that per-pose sulfur position is
  not persisted with the ligand. Recomputing against the rigid receptor's SG
  looked correct — every mode size matched exactly, 2/202/41/210 — and moved
  6 of mode 3's poses across the viability bar (10 → 16). The audit therefore
  joins the screen's own rows and validates the join two ways. If the archive is
  ever to be self-contained, `--all-poses` must persist the flexible sidechain
  beside the ligand.
* **`pose_rank` and `energy_rank` in a persisted SDF are the same write-order
  counter** (`nac_screen_v2.write_sdf` stamps `str(rank)` into both). `pose_rank`
  is used as an in-file identity and is sound for that; `energy_rank` names a
  quantity it is not, is unique and plausible over 1..N, and would be believed by
  anything that read it. Nothing reads it today. Tracked in #13.
