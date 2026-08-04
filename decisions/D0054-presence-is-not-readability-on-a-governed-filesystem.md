---
id: D0054
title: Presence is not readability — guards on the governed filesystem must test the access they depend on
date: 2026-08-04
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - integration/app/pose3d.py
  - tests/test_subpockets.py
  - tests/test_pose_modes.py
evidence:
  - 'measured 2026-08-04 for @tt8804: RECEPTOR.is_file() True, os.access(RECEPTOR, R_OK) False'
  - '13 tests failed with assertions about Arg-loop residues and pocket size; the cause was an ACL deny'
  - 'pocket_resi() fell back from the measured 8 A shell (>25 residues) to 11 declared ones, silently'
  - '0 of 166 sampled files under append_only/inhibition/ were readable at the time'
---

# A file you can stat is not a file you can read

## What was wrong

`needs_receptor` was `pytest.mark.skipif(not RECEPTOR.is_file())`. On this data
root a file can be present, `stat`-able and listed while every read of it
raises `PermissionError` — the Isilon ACL is invisible to an NFSv3 client and
denies per user, independently of the POSIX mode.

So the skip never fired. Thirteen tests ran against a receptor they could not
open and failed with assertions about **Arg-loop residues and pocket size** —
which reads as a broken pocket definition, not as an access problem. The same
defect sat in `_first_pose_file`, whose caller then skipped on "no docked poses
on disk" while a `PermissionError` was raised from inside the parser.

Worse, and quieter: `pocket_resi()` catches every exception and falls back to
the union of the declared sub-pockets. That fallback is right for a *missing*
shell file — a blank surface is worse than a small one — but it made the grey
surface silently shrink from the measured 8 Å shell (>25 residues) to **11**,
and render as a perfectly ordinary picture. A populated, plausible value
computed from the wrong thing.

## Why it looked right

`is_file()` is the idiomatic "is it there" check and it is correct on an
ordinary filesystem, where presence and readability differ only for root-owned
files nobody was going to touch. The whole class of bug only exists because
this project's data lives behind an ACL the client cannot query — and the
project already knew that, in prose, in §8 of `state_of_the_project.md`. The
knowledge existed; the guards predated it.

It also produced a **confident wrong diagnosis**. Anyone reading the failures
would have started debugging `SUBPOCKETS`, because that is what the assertions
name.

## The decision

Guards on the governed filesystem test **the access they depend on**, not
presence:

* `pose3d.receptor_readable()` — one predicate, owned by the module that owns
  the path, used by every caller.
* `_first_pose_file` returns the first pose file that can actually be opened.
* `pocket_resi()` keeps the fallback but **warns when the file exists and
  cannot be read**, because missing and unreadable are different facts and only
  the first is benign.

Skip reasons name the real cause ("absent, or denied by the data root's ACL")
rather than "not on this machine", which was untrue and sent the reader the
wrong way.

## Result

13 failures → 0. The suite reports 44 honest skips instead of 13 misleading
failures, and the skip text tells the reader what to fix.
