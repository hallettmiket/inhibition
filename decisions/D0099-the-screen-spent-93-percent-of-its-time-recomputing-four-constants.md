---
id: D0099
title: The screen spent 93% of its runtime recomputing four constants — 295 s per molecule became 7 s, and the PoseBusters gate D0089 adopted is now built
date: 2026-08-27
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/nac_criterion.py
  - shared/target_config.py
  - shared/pose_contacts.py
  - scripts/nac_screen_v2.py
  - config/target.yaml
  - decisions/D0089-posebusters-validates-poses-and-agrees-with-the-near-attack-criterion.md
  - decisions/D0093-the-file-named-allposes-is-not-all-poses-it-is-dbscan-cleaned.md
evidence:
  - '@tt8804: "can you build it in please. and also lets take a moment to do ultrathink and scrutinize the code and steps, maybe we can save some gpu time"'
  - 'PROFILED: nac_criterion.isotropic_null was 309 s of a 333 s per-molecule budget -- 93% -- called 384 times per molecule at 0.8 s each'
  - 'it is a PURE function of `mechanism`, of which there are four; the perpendicular branch integrates a 2001 x 4000 grid, and its own docstring says the result is "independent of which molecule is being scored"'
  - 'target_config.load re-parsed the YAML 196 times per molecule, 8 s, because `_cfg` reads single keys'
  - 'MEASURED BEFORE AND AFTER, values identical: isotropic_null returns 0.0670 (sn2) and 0.0816 (perpendicular) exactly as documented'
  - 'SCREEN: 295 s per molecule -> 7 s. 42x. The full 561-molecule screen goes from ~21 h to ~35 min on 2 GPUs'
  - 'AutoDock-GPU itself was never the cost: 500 runs take 2.0 s on an idle A100 (100 runs 1.6 s, 1000 runs 2.7 s)'
  - 'the splitting was never the cost either: predict_rmsf 3.2 s median across 8 molecules, contact tensor + pdist + linkage 2.2 s'
  - 'POSEBUSTERS BUILT: 93.0% valid in dock mode (22 checks, protein clashes included) against 99.8% in mol mode (12 ligand-only checks); D0089 measured 92.80% for non-attack-ready poses'
  - 'max_workers 8 takes it from 86 s to 20 s per 500 poses; 16 and 32 buy nothing more'
  - 'THE GATE FAILED CORRECTLY ON ITS FIRST RUN: handed the .pdbqt the docking uses, PoseBusters WARNS rather than raising and every protein-ligand check fails for want of a receptor -- verdict "0 of 500 valid", which reads as chemistry'
runbook: python scripts/nac_screen_v2.py --split-method contact_linkage
---

# D0099 — the screen was 93% overhead

## What the profiler found

Asked to look for GPU savings before committing to a re-screen, the answer was
that **almost none of the time was on the GPU.**

| | per molecule |
|---|---:|
| `isotropic_null` | **309 s (93%)** |
| everything else in `one()` | 16 s |
| — of which `target_config.load` | 8 s |
| AutoDock-GPU, 500 runs | **2.0 s** |
| `predict_rmsf` (50 conformers) | 3.2 s |
| contact tensor + pdist + linkage | 2.2 s |

`isotropic_null` is a pure function of `mechanism`. There are four mechanisms.
The perpendicular branch integrates a 2001 × 4000 grid to get one number, and the
screen called it **384 times per molecule**. Its own docstring already said the
result is *"independent of which molecule is being scored"*.

`target_config.load` re-parsed the YAML **196 times per molecule** because `_cfg`
reads one key at a time.

Both are now cached. **295 s → 7 s per molecule, 42×.** The 561-molecule screen
goes from ~21 hours to **~35 minutes on 2 GPUs**.

## Why the caches are safe, given this project's history

Three defects in the catalogue (#8, #9, #18) are caches keyed on less than their
inputs, so it is worth stating why these are not:

* `isotropic_null` reads `mechanism`, module constants, and nothing else. If a
  window constant ever becomes configurable the cache becomes wrong and must be
  keyed on it — written into the docstring.
* `target_config.load` is keyed on **(path, mtime_ns, size)**, not the path
  alone. A path-only cache would serve a stale config to a long run whose file
  was edited mid-flight, and two configs inside one run is the D0080 defect. It
  returns a **deep copy**, so a caller mutating the result cannot poison later
  readers — asserted in the verification.

## The measurement that stopped a bad conclusion

An earlier comparison had `contact_linkage` at 197 s/molecule against
`warhead_dbscan` at 10 s, which looked like the new splitter being 20× more
expensive. It was not: that run overlapped another session's MD on the same card,
and the splitting is 5 s either way. **The apparent cost of the new method was
another job's contention.** Timing on an idle box is the only reason this was not
recorded as a property of contact grouping.

## PoseBusters, finally built

D0089's Decision section reads *"Adopt the gate and the quota"*. It was never
implemented; `nac_screen_v2` had zero references to it. Now:

```yaml
docking:
  posebusters:
    enabled: true
    config: dock       # 22 checks incl. protein clashes -- what D0089 measured
    max_workers: 8     # 86 s -> 20 s per 500 poses; 16 and 32 buy nothing
```

Measured on three molecules: **91.6–94.8% valid**, against D0089's 92.80% for
non-attack-ready poses. `mol` mode passes 99.8% because its 12 checks cannot see
a protein clash at all.

**Invalid poses are flagged, never deleted.** They keep their per-pose row and
their place in the persisted cloud carrying `pb_valid = False`, and are excluded
only from grouping and the aggregates — `split_poses` takes a `conformers`
subset and labels the rest −1. Verified: 1,500 rows for 3 × 500 poses, 0 invalid
poses grouped, 0 valid poses ungrouped, and the cloud holds all 500.

**This also closes D0093 properly.** The all-poses writer read
`if labels[i] in mode_ids`, which is what dropped 21% of every nac_v5 cloud. It
now writes every pose, ungrouped ones carrying `mode = -1`, so a reader excludes
them by identity rather than by their absence.

## The gate failed correctly on its first run, for the wrong reason

Handed `3IKD_prepared_1.pdbqt` — the receptor the docking actually consumes —
PoseBusters **warns rather than raising** on an unreadable file, every
protein-ligand check then fails for want of a receptor, and the verdict is *"0 of
500 valid"*. That is an infrastructure fault wearing a result's clothes, and the
gate reported it as chemistry.

Two fixes, and the second matters more: the receptor is now required to be a
`.pdb` by suffix, **and a zero-valid verdict now names its own worst failing
checks and says to suspect the receptor.** A guard that reaches the right
conclusion by the wrong route is one input away from reaching the wrong one.

## Cost of the run, restated

| step | per molecule | 561 molecules, 2 shards |
|---|---:|---:|
| dock + split | 7 s | ~35 min |
| **+ PoseBusters** | **29 s** | **~2.3 h** |

PoseBusters is now the dominant cost and it is **CPU, not GPU** — this box has
224 cores, so more shards than GPUs is the right shape: the card is idle 93% of
the time.
