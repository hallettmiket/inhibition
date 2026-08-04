---
id: D0055
title: The orientation document's measured numbers are generated, and drift fails the suite
date: 2026-08-04
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - docs/state_of_the_project.md
  - scripts/refresh_orientation.py
  - tests/test_orientation_current.py
evidence:
  - 'the doc says of itself that it drifted badly within 24 h of being written and a new maintainer read it as fact'
  - 'on the day this landed the doc claimed 59,323 molecules across the six T_2 variants; the measured total was 60,123'
  - 'issue #11 asks for a device to keep the two orientation docs current'
---

# Tier-II memory has to be true or it is worse than nothing

## Context

Issue #11: *"these two files should give CC a quick way to get the context it
needs to move forward. (A type of tier II memory but in the repo not vault)."*

That framing is what makes the drift serious. These documents are not
documentation in the usual sense — they are the context a fresh Claude Code
session loads **before** touching anything, imported directly into `CLAUDE.md`.
A stale number is not an untidy doc; it is bad context handed to whoever works
next, with the authority of a file that says "start here".

The doc already knew: *"it drifted badly within 24 h of being written and a new
maintainer read it as fact."*

## Decision

**The counts are generated; the prose is not.**

`scripts/refresh_orientation.py` regenerates only the numbers that drift —
latest frame per experiment, rows, docked, ranked, shortlisted, per-seed pool
sizes, the decision count — into `<!-- AUTO:key:BEGIN/END -->` fences.

**"What is established", "what is ruled out" and "what to do next" are
deliberately untouched.** They are judgements. A generator that rewrote them
would produce a confident document nobody decided, which is a worse failure
than a stale number because it would not look stale.

A missing marker is a hard error, not a skip: a refresh that silently updates
nothing is indistinguishable from never running it.

## Why the script alone would not have worked

**A generator only helps if somebody runs it, and nobody runs it exactly when
the numbers are moving fastest** — which is when the drift happens. That is the
whole history of this problem.

So the load-bearing half is `tests/test_orientation_current.py`, which shells
out to `--check` and **fails the suite** when the document no longer matches
the frames. The script makes the fix a one-liner; the test makes skipping it
impossible. This is the same reasoning as the stale-pin guard: a pin cannot
announce that it is out of date, only a comparison against the directory can.

## It caught something immediately

First run: the hand-maintained line claimed **59,323** molecules across the six
T_2 variants. The measured total was **60,123** — 800 molecules of drift, in
the document whose stated purpose is to stop people being misled, four days
after it was written.

## Checking the guard can fail

Per rule 3 of `how_this_project_breaks.md`, the test corrupts a copy of the doc
and asserts `--check` rejects it, and separately asserts the expected AUTO keys
are present — so a refactor that renamed the markers cannot leave the check
inspecting nothing and passing for free.

Writing that test surfaced a real bug in the script: `_shown()` called
`Path.relative_to(REPO)`, which **raises** for a path outside the repo, so the
script crashed with a pathlib `ValueError` while formatting its own "this doc
is stale" message. Found only because the test drove it with a tmp path.
