---
id: D0097
title: The contact splitter is built, tested and selectable — and the five pose-grouping frameworks in this repo are mapped rather than merged
date: 2026-08-27
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/pose_contacts.py
  - shared/ligand_flexibility.py
  - shared/pose_modes.py
  - shared/pose_vector.py
  - shared/pose_cluster.py
  - shared/pose_subsplit.py
  - scripts/nac_screen_v2.py
  - config/target.yaml
  - data/reference/pocket_landmarks_1.csv
  - docs/pose_frameworks.md
  - tests/test_pose_contacts_split.py
evidence:
  - '@tt8804: "can we fully build out this pose splitting framework and clean up conflicting old frameworks in this repo. I am happy with the poses now"'
  - 'FIVE modules in this repo group poses: pose_modes, pose_subsplit, pose_contacts, pose_vector, pose_cluster -- and nothing said which one runs'
  - 'pose_vector (2026-08-04) already described a pose by what it touches; pose_contacts (2026-08-26) was written without discovering it'
  - 'they are NOT duplicates: pose_vector reduces to ONE NUMBER PER RESIDUE (orientation lost), pose_contacts keeps one per (atom, residue) (orientation kept)'
  - 'their linkage rationales directly contradict: pose_vector says complete linkage "would split them on the widest pair", pose_contacts says single linkage chains and that is D0088 defect. Both hold, for different n'
  - 'pose_cluster has ZERO production callers -- only exp/4,5,6,7,8,9,10 and their tests'
  - 'END TO END, same molecule and 200 poses: contact_linkage gives 49 groups (largest 29, worst within-group 0.519 A against a 0.52 A bound); warhead_dbscan gives 5 modes'
  - 'the landmark residues are now a versioned reference file resolved BY GLOB, with the water excluded once at write rather than by every caller'
  - 'predict_rmsf moved from exp/15 into shared/ligand_flexibility -- production code was importing it from an experiment directory by file path to dodge two modules both named run_all'
  - 'every frame now stamps split_method, split_tolerance_a and split_landmarks'
  - '17 new tests; suite 1,041 passed'
runbook: nac_screen_v2.py --split-method contact_linkage ; see docs/pose_frameworks.md
---

# D0097 — built, selectable, and mapped

## What was built

`shared/pose_contacts.split_poses` is now a production entry point and a drop-in
for `pose_modes.split` + `pose_subsplit.subdivide` — one integer label per
conformer, so rank_v2, the sweep, the GUI and `mode_key` need no changes.

Four things had to move out of experiment code first:

| | was | now |
|---|---|---|
| landmark residues | an `exp/14` output glob under `append_only/` | `data/reference/pocket_landmarks_1.csv`, resolved by glob |
| RMSF predictor | `exp/15/run_all.py`, imported by file path | `shared/ligand_flexibility.predict_rmsf` |
| tolerance | recomputed in five experiments | `pose_contacts.tolerance_for` |
| method choice | implicit — whatever the caller imported | `config/target.yaml: splitting.method` |

The landmark file drops `A:40:HOH` **once, at write**, rather than making
"remember to exclude water" a rule every reader carries.

## It is selectable, and the default has not moved

```yaml
splitting:
  method: warhead_dbscan     # or contact_linkage
```

`--split-method` overrides it per run, as an **allowlist** — argparse refuses
`contact_linkag` rather than falling through to the shipped rule, which would
produce artefacts naming a method they did not use (#14's denylist defect).

**The default is still `warhead_dbscan`, deliberately.** Every frame, ranking,
sweep result and 100 ns run on disk was produced under it; switching re-groups
every cloud and makes them incomparable in a way nothing in the artefacts would
announce. The switch belongs with the re-screen (#79).

Measured end to end on one molecule and 200 poses:

| method | result |
|---|---|
| `contact_linkage` | **49 groups**, largest 29, worst within-group 0.519 Å against a 0.52 Å bound |
| `warhead_dbscan` | **5 modes** (one mode sub-split into five) |

Ten times the count from the same poses — which is D0092 in one line, and the
reason the count must never be reported as a number of binding modes.

Every frame now stamps `split_method`, `split_tolerance_a` and
`split_landmarks`. Two rules that both emit integer labels are otherwise
indistinguishable downstream, and every mode count, enrichment and rank in a
frame means something different depending on which one ran.

## The cleanup — and the framework I nearly duplicated

There are **five** modules in this repo that group poses. The honest finding is
that `shared/pose_vector.py` (2026-08-04) already described a pose by what it
touches, and `pose_contacts` was written three weeks later without discovering
it.

They are not duplicates, and the differences are the point:

| | `pose_vector` | `pose_contacts` |
|---|---|---|
| granularity | one number **per residue** | one per **(atom, residue)** |
| orientation | **lost** — a flip touching the same residues is identical | **kept** |
| linkage | **single** | **complete** |
| sized for | ~9 Vina modes (O(n³) loop) | 500–6,000 poses (`pdist`) |
| for | a fit score against a reference profile | splitting a cloud |

**Their linkage rationales contradict each other in the source, and both are
right.** `pose_vector` says complete linkage "would split them on the widest
pair"; `pose_contacts` says single linkage chains, which is exactly D0088's
137-pose mode spanning 9.3 Å. Over 9 modes already spread by a minimum-RMSD
floor, chaining cannot run away. Over 500 poses filling a continuous cloud it is
what happens. **The linkage must match the density of the cloud**, and neither
module's reasoning transfers to the other's input.

So the resolution is a map, not a merge: `docs/pose_frameworks.md`, plus a
`STATUS:` banner at the end of each of the five docstrings naming what it is for
and what it is not. `pose_cluster` (HDBSCAN, D0090) is marked **superseded** —
it has no production caller, only the experiments that are the record of why it
was rejected — and is listed in `data/ready_to_delete.md` pending their archival.

## What is still not settled

* **The tolerance** is close to a constant and does not beat writing one number
  down (D0094). `tolerance_from_descriptors` is implemented and unadopted.
* **The count does not saturate** (b = +0.69, no plateau) — a sampling statistic,
  not a mode count.
* **Never tested against the MD-validated pose**, which is the one measurement
  connecting splitting to ground truth, and cheap.
* **SIFt (Deng 2004) has not been run.** It is the published prior art for
  describing a pose by its interactions, and until it is measured here this is a
  re-derivation with an unquantified advance.

## Why the sixth framework would have been easy to write

Because nothing in the repo said the fifth existed. Five modules, five docstrings
each excellent about its own reasoning, and no index. The banner and the map are
the actual fix; `pose_contacts` being better than `pose_vector` at splitting is
incidental to it.
