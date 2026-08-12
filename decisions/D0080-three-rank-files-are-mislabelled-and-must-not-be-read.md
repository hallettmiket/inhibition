---
id: D0080
title: A run's topic is a directory, never a behaviour flag — three mislabelled rank files, and the defect family that produced them
date: 2026-08-12
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/rank_v2.py
  - scripts/nac_screen_v2.py
  - scripts/rebuild_representatives.py
  - scripts/sweep_gap_worklist.py
  - scripts/pose_modes_report.py
  - shared/mode_ranking.py
  - shared/mode_assets.py
  - shared/target_config.py
  - config/target.yaml
evidence:
  - 'rank_v2 resolved input as `V3 if topic == "nac_v3" else V2`, so --topic nac_v4 read the 2.1.0 aggregates and wrote them under a nac_v4 filename, exit 0'
  - 'the wrong output holds 1,683 + 4,082 rows one per MOLECULE; a 3.0.0 ranking holds 34,076 one per MODE, and carries no `mode` or `parent_ident` column at all'
  - 'nac_screen_v2 derived its table topic from --topic but hardcoded POSE_DIR/ALL_POSE_DIR to nac_v3_*, so the append-only guard skipped every write: 5,772 of 5,774 representative files were still the Aug-07 2.1.0 ones'
  - 'that is the same defect the 2.2.0 write-off was diagnosing ("pose counts differ 2 to 15 from the score table") -- treated as a data problem, so it was never fixed'
  - 'recovery needed no re-docking: #44 had persisted 5,420 clouds whose per-pose mode stamps matched the tables for 5,336 of 5,336 molecules and disagreed for none'
  - 'the rebuild reproduces the screen`s own representatives bit-identically -- 0.000000 A across every mode of a molecule re-run through the fixed screen'
  - 'mode_assets, pose_modes_report and sweep_gap_worklist each hardcoded nac_v3_poses, so GUIs and worklists would have drawn one run`s poses beside another run`s numbers'
  - 'per-mode is now detected from the aggregate`s own columns, and run.topic in config/target.yaml names the production run for every reader'
---

# D0080 — three `nac_v4` rank files contain `nac_v2` data and must not be read

## The files

```
rank_v2_T4_nac_v4_weighted_score_1.csv
rank_v2_T3_nac_v4_weighted_score_1.csv
rank_v2_REF_nac_v4_weighted_score_1.csv
```

They are named for the 3.0.0 run and contain the 2.1.0 one. **Nothing may read
them.** They are not deleted because they sit under `append_only/`, and the rule
that keeps a five-hour screen safe is worth more than the tidiness of removing
three files. This decision is the marker.

## How they were produced

`rank_v2` resolved its input as `V3 if args.topic == "nac_v3" else V2`. Every
topic that was not literally `nac_v3` — including `nac_v4`, the whole 3.0.0
screen — read the 2.1.0 aggregates out of `nac_v2` and wrote the result under the
*requested* topic's filename. The run exited 0.

The output is recognisable once you know: it holds 1,683 + 4,082 rows, one per
MOLECULE, where a 3.0.0 ranking has 34,076, one per MODE. It carries no `mode`
and no `parent_ident` column at all.

Two further gates keyed on the same literal, so a topic that was not `nac_v3`
also silently skipped per-mode ranking and the recalibrated consensus gate. All
three are fixed: the topic now resolves to its own directory, and per-mode is
detected from the aggregate's own columns rather than from the topic's name.

## Why it matters more than a bad file

This is the third instance of one defect found on 2026-08-12, and the family is
what to remember, not the individuals:

| where | what was keyed on a topic literal | consequence |
|---|---|---|
| `nac_screen_v2` | the two pose directories | 3.0.0 wrote no representative poses |
| `rank_v2` | the input directory + two behaviour gates | 3.0.0's ranking was 2.1.0's data |
| `mode_assets`, `pose_modes_report`, `sweep_gap_worklist` | the pose directory | GUIs and worklists read the wrong run's poses |

Each was silent, each exited 0, and each produced output that looked exactly
like the thing it was supposed to be. The 2.2.0 write-off — "pose clouds and
score tables came from different runs" — was the first of them, diagnosed as a
data problem and therefore never fixed.

## The rule that replaces them

A run's topic is a **directory name**, resolved once, and every artefact of that
run is derived from it — tables, representative poses, pose clouds, rank files,
GUI assets. There is no privileged un-suffixed filename meaning "the current
screen": `run.topic` in `config/target.yaml` says which run is production, and
the ranking stamps it into every file it writes.

Behaviour must never be selected by a topic's name. Where a stage behaves
differently for per-mode data, it asks the data (`"parent_ident" in columns`),
which is self-describing and cannot go stale.

Guarded by `tests/test_topic_paths.py`.
