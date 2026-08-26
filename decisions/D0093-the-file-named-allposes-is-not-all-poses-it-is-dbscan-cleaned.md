---
id: D0093
title: The file named `allposes` is not all poses — it is the DBSCAN-cleaned subset, so every replacement for DBSCAN has been measured on clouds DBSCAN already cleaned
date: 2026-08-26
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/nac_screen_v2.py
  - exp/14_residue_selection/run_all.py
  - exp/15_rmsf_predictor/run_all.py
  - exp/16_contact_clustering/run_all.py
  - exp/17_contact_saturation/run_all.py
  - docs/how_this_project_breaks.md
evidence:
  - 'nac_screen_v2.py:501 writes `order = [i for i in argsort(labels) if labels[i] in mode_ids]`, and nac_screen_v2.py:350 sets `mode_ids = sorted(set(labels) - {-1})` -- DBSCAN noise is excluded from the file'
  - 'docking.n_runs = 500 against 393 poses in the production cloud for t4_716800c125a7: 21% of the cloud is absent'
  - 'MEASURED CONSEQUENCE: centroid extent 19.19 A in the raw 6,000-pose cloud against 7.14/7.09/7.22 A in three independent 500-run production clouds -- dropping 21% of poses removes 62% of the spatial spread, because the dropped poses are the scattered ones'
  - 'at n ~ 400 the raw cloud gives 241-254 contact-space groups and the filtered clouds give 109-118 -- a factor of 2.2, entirely from poses that would have been singletons'
  - 'exp/5 deep_dock knew: "Not the screen `--all-poses` file: that one drops DBSCAN noise, and a saturation curve needs the raw cloud" -- the caveat existed in one docstring and travelled with nothing'
  - 'RULED OUT FIRST: both paths call ns.build_reactive_receptor(ns.RX_RECEPTOR), one cached 26 A box and one receptor, so the difference is not the box (checked before the filter was found)'
runbook: null
---

# D0093 — `allposes` is not all poses

## What was taken, and what it should have been

`shared/run_paths.allposes_dir()` resolves `<topic>_allposes/`, and every
experiment that needed "the pose cloud" read it. It is not the pose cloud. At
[`nac_screen_v2.py:501`](../scripts/nac_screen_v2.py#L501):

```python
order = [int(i) for i in np.argsort(labels, kind="stable")
         if labels[i] in mode_ids]
```

with `mode_ids = sorted(set(labels) - {-1})` thirty lines earlier. **Poses DBSCAN
labelled noise are not written.** For `t4_716800c125a7` that is 107 of 500 — 21%
of the cloud, and specifically the 21% that failed to join a dense region.

## How it surfaced

Not by audit. exp/17 compared five independent 500-run dockings against
subsamples of the raw 6,000-pose deep cloud at the same n, as a check that
subsampling a pooled cloud is a fair stand-in for an independent one. The
independent clouds returned 109–118 groups where the subsamples returned 241–254
— a factor of 2.2 in the wrong direction for any sampling artefact.

The first hypothesis was the docking box, since `config/receptor.yaml` gives T_4 a
20 Å covalent box while `nac_screen.py:130` hardcodes 26 Å. **That hypothesis was
wrong and was checked before being reported:** both paths call
`ns.build_reactive_receptor(ns.RX_RECEPTOR)`, which is one cached receptor and one
box. The centroid extents then made the real cause plain — 19.19 Å raw against
7.1 Å filtered. Removing a fifth of the poses removes nearly two-thirds of the
spread, because the poses removed are the outlying ones.

## Why it looked right

The name. `allposes` is a promise, `--all-poses` is a flag that reads as "keep
everything", and #44 — *every docked pose is persisted, always* — is a standing
rule in `CLAUDE.md` that the directory appears to implement. Nothing about reading
that path suggests a filter has already run. The counts are plausible (393 poses
is a believable cloud), and no downstream number looks wrong.

This is disguise #1, **selection by name**: a column — here a directory — that
*describes* what you want rather than being derived from the thing that defines
it. It is entry #1 of the catalogue (`shortlist` vs `shortlist_synth`) with a
directory in place of a column.

## Why it matters more than the usual instance

The filter is the clustering we are trying to replace. exp/14 (which residues to
use), exp/15 (RMSF calibration) and exp/16 (does contact grouping produce tight
groups?) all read this path. **Every candidate replacement for DBSCAN has been
evaluated on clouds DBSCAN had already cleaned** — the poses hardest to group are
exactly the ones absent from the test. exp/16's headline numbers (median 174
groups, 47% singletons, worst within-group RMSD 2.71 Å) are therefore measured on
the easy case; the raw cloud has roughly twice the groups.

The knowledge existed. `exp/5_mode_saturation/run_all.py:139` says so in a
docstring — *"Not the screen's `--all-poses` file: that one drops DBSCAN noise,
and a saturation curve needs the raw cloud"*. It travelled with nothing: no guard,
no name change, no note on the reader.

## Fix the class, not the case

1. **Name the file for what it holds.** `<topic>_modeposes/`, with `allposes`
   retained only if it is genuinely written unfiltered.
2. **Persist the raw cloud too**, which is what #44 actually asks for, and stamp
   each pose with its DBSCAN label so the filter is a query rather than a
   deletion. A pose set that has been filtered cannot announce it; a label column
   can.
3. **Make the reader state its assumption.** `run_paths` should expose
   `mode_assigned_cloud()` and `raw_cloud()` as separate calls that raise when the
   requested one is absent, so no caller gets the other by default.
4. **Re-run exp/14–16 against raw clouds** before any of their numbers are cited.
   exp/17's conclusions (D0092) are unaffected: they are measured on the raw
   6,000-pose cloud throughout, which is why the discrepancy was visible at all.

## Guard

A test asserting that the persisted cloud for a molecule holds `docking.n_runs`
poses, or that its per-pose label column contains at least one `-1`. Both fail
loudly on a silently filtered file, and neither can pass by the thing it inspects
being absent.
