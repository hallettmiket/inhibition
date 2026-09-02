---
id: D0109
title: Every NAC distance and angle was measured against pose 1's Cys113 sulfur, because the function whose docstring forbids that returned the first match
date: 2026-09-02
status: accepted
approach: shared
decided_by: '@twu383'
origin: user
supersedes: []
superseded_by: null
affects:
  - scripts/nac_screen.py
  - scripts/nac_screen_v2.py
  - shared/nac_criterion.py
  - tests/test_per_pose_sulfur.py
  - config/target.yaml
evidence:
  - '`sg_position(dlg)` docstring: "the flexible sidechain moves during docking, so its position must come from the pose being measured" -- and the body `return`s inside the loop on the FIRST CYS SG line, i.e. model 1'
  - '`measure_poses(mol, match, mechanism, sg)` took ligand coordinates PER CONFORMER and one `sg` for all of them, broadcasting pose 1''s sulfur across every pose'
  - 'MEASURED on t4_084b1552486e, 200 poses: SG displacement from pose 1 median 0.18 A, max 1.88 A, 6% of poses beyond 0.5 A'
  - 'same 200 poses, same ligand atom, only the sulfur choice differing: min distance 1.45 A (code) vs 2.88 A (per-pose); 78% of poses shift >0.1 A; 11% change in-window classification'
  - '7 of 200 poses came out below 1.81 A -- shorter than a C-S single bond -- which is the only reason the defect became visible'
  - 'BEFORE/AFTER on the identical cloud (seed 42, n_runs 640, same 208-group split, one molecule): min 2.13 -> 2.94 A; median 4.07 -> 4.51; poses <3.0 A 96 -> 7; in-window 199 -> 238; VIABLE 49 -> 96'
  - 'the error is one-directional and reported distances too SHORT, so poses fell below the 2.8 A window floor and `viable` was UNDER-counted by roughly half'
  - 'PoseBusters passed poses the distance column called 1.22 A: PB checks the real pose geometry and was correct; the two were measuring different sulfur positions'
runbook: null
---

# Every NAC distance was measured against pose 1's sulfur

## What happened

Cys113 is docked as a **flexible** sidechain, so every pose has its own SG
position. `nac_screen.sg_position(dlg)` returned the first `CYS ... SG` record in
the DLG — model 1's sulfur — and `nac_criterion.measure_poses` broadcast that
single value across all 640 conformers while taking ligand coordinates per
conformer. **Poses 2..N were measured against where the sulfur used to be.**

The function's own docstring states the correct rule:

> *"The flexible sidechain moves during docking, so its position must come from
> the pose being measured — reading it from the rigid receptor would measure the
> approach to where the sulfur STARTED."*

It then reads it from pose 1 and hands it to every pose. The docstring is not
wrong about the principle; the body applies it once.

## Why it looked right for so long

**The error is small, systematic and one-directional.** SG moves a median of
0.18 Å between poses, so most distances were plausible and slightly too short.
Nothing was out of range, nothing raised, and a distance of 3.2 Å looks exactly
like a distance of 3.4 Å.

**It only broke the surface at the tail.** Ranking on minimum distance put
molecules at 1.22–1.45 Å at the top, and a C–S single bond is 1.81 Å. An
impossible number is what finally made it visible — and only because a ranking
was asked for that had no lower bound. Every previous use of `distance` went
through the 2.8–4.2 Å window, which *hides* the sub-window tail by design.

**PoseBusters was disagreeing the whole time and was right.** 6,759 poses in
nac_v7 sat below PB's documented 2.625 Å clash floor and PB passed them,
including one the distance column called 1.22 Å. PB checks the actual pose
geometry; there was no clash. The two were measuring different sulfurs, and the
natural reading — "PB's clash check must be broken again, like #31" — was
backwards.

## The direction of the bias matters

Because distances were reported **too short**, poses fell below the window's
2.8 Å floor and were marked not-viable. Fixing it does not lower the scores; it
**raises** them. On the one molecule measured before and after on an identical
cloud, `viable` went from 49 to 96 — nearly double — and in-window from 199 to
238. So `viable_fraction`, `enrichment`, `anchor_quality` and therefore
**`engagement`** were all systematically depressed, and the depression was not
uniform: it depends on how far each molecule's flexible sulfur happened to
wander, which is a property of the docking run and not of the molecule.

That is the part that makes it a ranking defect rather than a calibration
offset.

## Decision

1. **`nac_screen.sg_positions(dlg)`** returns `(n_models, 3)` in MODEL order —
   the same order `rebuild_and_match` and `pose_energies` use.
2. **`measure_poses` refuses a single sulfur for a multi-conformer molecule.**
   This is the fix; the call-site edits are only today's instance. A bare `(3,)`
   raises unless `allow_static_sg=True` is passed, so a static sulfur becomes a
   decision written at the call site instead of a silent broadcast. A wrong-length
   array raises too — pairing pose j with pose k's sulfur would be worse than the
   original defect.
3. **`sg_position` is deprecated, not deleted**, and its docstring names this
   record. Deleting it would turn stale callers into an `AttributeError`
   somewhere that does not explain itself.
4. **nac_v7 is retired as broken**, not superseded. Its tables stay on disk
   because the root is append-only; nothing may read them. The re-screen goes to
   **nac_v8** rather than back into nac_v7, because `rank_v2` concatenates every
   `agg_s*.csv` in a topic and would count every mode twice.

## What this invalidates

Every `distance`, `angle`, `in_range`, `viable`, `viable_fraction`,
`enrichment`, `anchor_quality` and `engagement` produced since Cys113 became
flexible — which includes **nac_v5, nac_v6 and nac_v7**, and the rankings built
on them. The re-screen is cheap (seed 42 makes it reproducible, ~40 min on 4
GPUs for 562 molecules) and reproduces the identical poses with correct
geometry.

It does **not** invalidate:

- **the modes.** Contact linkage splits on weighted residue-contact distance and
  never touches distance-to-SG, so group membership is unaffected. Verified: the
  before/after smoke runs produced the same 208 groups at the same tolerance.
- **the tier-1 MD readout** (D0071, D0108). That measures warhead-to-SG in the
  MD frame via `elevation_run.distance_nm`, a different code path against an
  explicit Cys113 SG in a static receptor. D0108's NO GO stands.
- **PoseBusters validity.** It was independent and correct throughout.

## What else shares the shape

The lesson is narrower and sharper than "read the docstring". It is: **a
function that returns one value where the physics has many will be called once
and broadcast**, and the broadcast is invisible because the shapes are
compatible. Worth auditing every place a per-pose quantity is derived from a
whole-run artefact: `_reactive_xyz` has the identical structure ("coordinates of
the retyped reactive atom, **from the first docked pose**") and is used by other
paths; `pose_energies` correctly returns one per model, which is the shape to
copy.
