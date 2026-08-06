---
id: D0060
title: Vina-GPU segfaults when it cannot write its kernel cache, and it fails silently
date: 2026-08-05
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/noncovalent_dock_run.py
  - scripts/redock_3ikd_benchmark.py
evidence:
  - 'every Vina-GPU invocation as @tt8804 segfaulted (exit 139 / -11): 6VAJ with its own box, 3IKD, and 3IKD with waters stripped, identically'
  - 'the governed wrapper cds into /data/lab_vm/envs/_src/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1, owned by mhallet, not writable by other members'
  - 'the log shows "Build kernel 1 from source" — it recompiles OpenCL kernels each run and writes Kernel1_Opt.bin / Kernel2_Opt.bin into that directory'
  - 'the cached kernels are READABLE (test -r succeeds); only the write fails'
  - 'ATRA degree-2 docking reported driver exit 0 while all six chunks failed: 127 of 30,000 molecules docked'
  - 'byte-identical copy in a writable directory: exit 0, 3 of 3 poses, 17.1 s'
---

# A permissions failure that presented as a docking result

## What happened

Every Vina-GPU run failed with a segmentation fault. The obvious readings were
all wrong: it was not the new 3IKD receptor (6VAJ failed identically with its own
box), not the retained waters (stripping them changed nothing), and not the box
geometry (verified inside the receptor's extent, all atom types standard).

The governed wrapper `cd`s into a source tree owned by `@mhallet` that other lab
members cannot write. Vina-GPU **recompiles its OpenCL kernels at startup** and
writes `Kernel1_Opt.bin` / `Kernel2_Opt.bin` into its working directory. When
that write fails it segfaults, rather than falling back to the cached kernels it
can **read perfectly well**.

## Why it cost something

**It failed silently at the level anyone was watching.** The ATRA degree-2
docking run reported driver **exit 0** while all six chunks failed:

```
docked successfully 127 / 30000
CHUNKS FAILED       6
```

and the pose directory still held the **previous** run's 15,653 files, so it
looked populated and plausible. That was reported as "docking finished, exit 0"
and believed, because the exit code of the wrapper was checked instead of the
count of what it produced. Catalogue shape: a value taken from the wrong place,
failing silently because both the right and the wrong candidate were populated.

The lesson is narrow and repeatable: **a driver's exit code is not evidence that
the work happened.** `dock_chunked` already prints the count — it was printed and
not read.

## The decision

A **byte-identical** copy of the install lives in a directory we can write, and
`shared/noncovalent_dock_run` resolves between them:

```python
if os.access(workdir, os.W_OK):  governed wrapper
else:                            the writable copy, with a warning
```

Resolved at import, by capability rather than by a hardcoded switch, so the
governed path returns automatically once `/data/lab_vm/envs/` is writable and
nobody has to remember to change it back.

**The binary is the same file** — sha256 `a53d33554320d41b0ac22d9d...`, verified
at copy time. Only the working directory differs, so this is a permissions
workaround and not a different tool version. That distinction is what keeps the
results comparable to everything docked before.

## Proper fix

Write access to `/data/lab_vm/envs/` for `ssmd-ud-vmlab`. It is the one tree
John's 2026-08-05 permission change did not cover, and the same fix that
unblocked `immutable/` and `append_only/` applies. Retire the copy when it lands
— the resolver will pick the governed wrapper up on its own.

## What this invalidates

Any docking run by a non-owner since the permissions changed. Known: the ATRA
degree-2 campaign (127 of 30,000). The 3IKD pose-recovery benchmark had not
produced a number at all, so nothing was concluded from it.
