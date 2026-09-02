---
id: D0105
title: BPMD ran pose_rank 1 for every molecule because no caller ever passed the argument its own reader supports
date: 2026-08-31
status: accepted
approach: shared
decided_by: '@twu383'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/bpmd_run.py
  - tests/test_bpmd_pose_identity.py
evidence:
  - '`read_pose(sdf, pose_rank=1)` has selected by the `pose_rank` PROPERTY since it was written, and its docstring says why: "never by position ... BPMD on the wrong pose completes normally and reports a perfectly plausible stability"'
  - '`run_pose(..., pose_rank: int = 1)` accepts and forwards it, and `prepare_pose` uses it; the workdir is even named `<stem>__p<rank>` for rank != 1'
  - 'BOTH call sites in `main()` omitted it: `run_pose(cand, replicates=..., production_ps=..., gpu=..., threads=..., nrun=..., dock_gpu=..., allow_redock=..., on_row=chunk.add)` -- no `pose_rank`, so every production and convergence run took the default 1'
  - 'there was no `--pose-rank` flag on the CLI at all, so no operator could have passed one'
  - 'the mode representatives are multi-pose: `elev_sulfone_vs_cf2_poses/t4_80fbed3bdf1e.sdf` holds 165 poses, one per contact-linkage mode'
  - '`already_done()` keyed on `(ident, replicate)`; `load_results()` de-duplicated on `["ident", "replicate"]`'
  - 'the 100 ns run launched the same day used pose_rank 11 (mode 10, the top mode by engagement, 0.7247); an unfixed BPMD would have measured pose_rank 1 (mode 0, engagement 0.0000) and filed it under the same molecule name'
runbook: null
---

# BPMD ran pose_rank 1 for every molecule

## What happened

Elevating `t4_80fbed3bdf1e` meant running a 100 ns trajectory and BPMD **on the
same pose**, so the two readouts describe one binding mode. `md_residence_3ikd.py`
takes `--pose-rank` and `--mode`. `bpmd_run.py` did not.

Not because it could not: `read_pose` selects by the `pose_rank` property,
`run_pose` takes the argument, `prepare_pose` honours it, and the workdir is
already named `<stem>__p<rank>` for ranks other than 1. The whole mechanism was
built. **Both call sites in `main()` simply did not pass it, and no CLI flag
existed to supply one.** Every BPMD replicate this project has ever run was on
pose_rank 1.

## Why it looked right

Three things kept it invisible, and they are the usual three.

**The default is a legal value.** `pose_rank=1` is a real pose of the right
molecule. It parameterises, solvates, biases and reports a perfectly ordinary
stability score. Nothing is malformed and nothing raises.

**The evidence of intent was in the reader, not the caller.** `read_pose`'s
docstring argues at length that selecting by position "would begin biasing a
different pose the day that changes", and warns that "BPMD on the wrong pose
completes normally and reports a perfectly plausible stability". That warning is
about `mols[0]`. It did not notice that its own caller was pinning the rank to a
constant, which has the identical effect by a different route.

**The pose set changed underneath it.** The contract used to be one pose per
file — `read_pose` still carries the comment recording that "one pose per file
used to be the export's contract and is not any more". When the export became
multi-pose, `read_pose` was fixed to cope. The call sites were not revisited,
because they still worked.

## Why it matters here specifically

The top mode of `t4_80fbed3bdf1e` by the config's own score is **mode 10**
(engagement 0.7247) at pose_rank 11. Mode 0 — pose_rank 1 — has engagement
0.0000. Had this not been caught, the run would have produced a BPMD number for
mode 0 and a 100 ns number for mode 10, written them under one molecule name,
and invited exactly the comparison they cannot support. That is
`how_this_project_breaks` **#23** almost verbatim: an asset matched on the
MOLECULE where the runner writes one per `(molecule, pose_rank)`.

## The second half: the resume key

`already_done()` returned `(ident, replicate)` pairs. With a `--pose-rank` flag
and no change here, a completed run of pose_rank 1 would mark pose_rank 11 as
already done — the molecule matches, the replicate matches — and the second pose
would **never be simulated while the table said it had been**. That is
`how_this_project_breaks` **#22**, which is this same function, keyed on
trajectory length instead. The fix that adds the flag has to fix the key in the
same change, or it installs a worse bug than it removes.

`load_results()` de-duplicated on `["ident", "replicate"]` for the same reason
and had the same blind spot: two poses of one molecule are two results, and one
of them was being dropped.

## Decision

1. `--pose-rank` on the CLI, defaulting to 1, threaded to **both** `run_pose`
   call sites.
2. `already_done()` returns `(ident, pose_rank, replicate)`; rows written before
   `pose_rank` was recorded read as rank 1, which is what they were — no caller
   could produce anything else.
3. `load_results()` de-duplicates on `["ident", "pose_rank", "replicate"]`.
4. A test asserts the argument actually reaches the runner and that the resume
   key distinguishes two poses of one molecule — not that the flag parses.

## What else shares this shape

The question this raises is not "where else is `pose_rank` missing" but **which
other functions accept an identity argument that no caller passes**. A default
that is a legal value, on a parameter the code went to the trouble of
supporting, is invisible in exactly the way a pinned version literal is: it was
right when written, it cannot announce that it is not, and the output looks
normal. Worth a sweep for keyword arguments with defaults that are never
overridden anywhere in the repo.
