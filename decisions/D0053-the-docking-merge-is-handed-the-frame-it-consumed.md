---
id: D0053
title: The docking merge is handed the frame it consumed, never re-resolves it
date: 2026-08-04
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/noncovalent_dock_run.py
  - tests/test_dock_merge_provenance.py
evidence:
  - '2a22970 split merge_poses_onto_frame out of run(); the frame_path reference moved, the binding did not'
  - 'frame_path is assigned only inside the `if df is None:` branch; the manifest call reads it unconditionally'
  - 'run() passes df, so the T_1/T_2 03_dock.py route raised UnboundLocalError at the manifest write'
  - 'dock_chunked.py passes no df, binds it, and works -- the only pool docked since the refactor was liu_2024_c3, which requires chunking'
  - 'reproduced with symtable (frame_path is_local=True, is_parameter=False) and at runtime'
---

# The merge cannot name a frame it was not given

## What was wrong

`merge_poses_onto_frame` was split out of `run()` in `2a22970` so a chunked run
and a single-GPU run would share one merge — correctly, because two code paths
that both "merge the docking results" is how the covalent and GROMACS frames
acquired suffixed columns nobody noticed.

The reference to `frame_path` came across. The **binding** did not:

```python
if df is None:
    frame_path = dio.latest(...)      # only bound here
    df = dio.read_frame(frame_path)
...
inputs={"frame": frame_path, ...}     # read unconditionally
```

So the chunked path (no `df`) binds it and works, and the `run()` path (passes
`df`) raises `UnboundLocalError` — **at the manifest write, after the entire
GPU run is already spent.** ~1.3 h for atra, ~10 h for du_xu.

## Why it looked right

Three things hid it, and all three are worth naming.

**It was tested by the path that works.** The only pool docked after the
refactor was `liu_2024_c3` at 16,806 molecules — large enough to require
chunking. The chunked driver is precisely the caller that does *not* pass `df`.
The refactor was exercised end to end, successfully, on the one branch that
binds the variable.

**It fails at the end, not the start.** Nothing about the run looks wrong until
the last statement. A guard that fired at launch would have cost seconds.

**It is a crash, not a wrong number** — which makes it the *unusual* case for
this project, and is why it survived a code review that was looking for the
familiar shape. `how_this_project_breaks.md` catalogues twenty-one silent
defects; this one is loud, and loudness made it feel like something that would
have been noticed already.

## The decision

**The caller that supplies `df` must also supply `frame_path`.** Absent, the
merge raises rather than proceeding.

**It is explicitly NOT fixed by calling `dio.latest` when `df` was supplied.**
That is the obvious one-line fix, it makes the crash disappear, and it is
wrong: it would record whichever frame is newest *at merge time* rather than
the one the caller actually read. A manifest exists to record the SHA-256 of
every input a run consumed, so a manifest naming the wrong frame is worse than
a crash — it is the silent-provenance failure this project already has five
instances of. The caller knows which frame it read. It passes it.

## The guard

`tests/test_dock_merge_provenance.py`:

* a runtime test of the `run()` call shape, asserting `write_full_frame` was
  **actually reached** and recorded the frame passed in — not merely that
  nothing raised, because a merge returning early would otherwise pass for
  free (the `test_stale_guard` lesson);
* the refusal path, asserting no frame is written;
* an **AST guard** over `shared/ scripts/ approaches/ integration/ tests/`
  failing for any *future* call passing `df` without `frame_path`;
* a test that the AST matcher registers the defect when fed it, so the class
  guard cannot pass vacuously.

Verified to fail against the pre-fix module and pass after.
