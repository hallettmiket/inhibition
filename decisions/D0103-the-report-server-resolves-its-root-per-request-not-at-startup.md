---
id: D0103
title: The report server resolves its root per request, because startup resolution is a pin with a 14-day half-life
date: 2026-08-31
status: accepted
approach: integration
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/serve_reports.py
  - tests/test_report_server_follows_topic.py
evidence:
  - 'MEASURED 2026-08-31: the server on :8931, up 14 days, returned mdprio_reports_nac_v5/index.html byte-for-byte (13,423 bytes) while run.topic was nac_v6 and mdprio_reports_nac_v6/index.html (13,448 bytes) was being rebuilt beside it every few minutes'
  - ':8931 is the DEFAULT port and the one the docs tell you to forward -- docs/gui_spec.md:166 `ssh -L 8931:127.0.0.1:8931`, docs/pipeline_in_plain_language.md:301'
  - 'audit of every running report server: :8031 nac_v6 (up 1d18h), :8931 nac_v5 (up 14d), :8950 nac_v6 (up 3d14h), :8933 galena_3.0.0 archive (by design), so three of four live servers were right only because they happened to start after the topic bumped'
  - 'serve_reports.py module docstring claimed immunity: "the root is resolved from run_paths rather than typed, so the server cannot end up serving a superseded topic''s pages -- which it did once, for hours, from a literal path"'
  - 'target_config.load is cached on (path, mtime_ns, size), so per-request resolution costs a dict deepcopy and picks up a topic change on the next request'
runbook: null
---

## Context

Asked to launch the current ranking GUI, the first thing to check was whether
the one already running was showing the current run. It was not. The server on
**:8931 had been up 14 days and was serving `nac_v5`** while `run.topic` had
been `nac_v6` for days and `mdprio_reports_nac_v6/` was being rebuilt beside it
every few minutes by a live `pipeline.write_status()` loop.

Nothing about it looks wrong. The pages render, the title says *"DWI covalent
screen — Home"*, the nav works, and the numbers are a real screen's numbers —
just the previous screen's. This is catalogue #25 (*"the report server had the
same defect, serving the old directory from a literal path"*) recurring in the
component written to close it.

## Decision

`LiveRun.translate_path` re-resolves `rp.reports_dir()` on **every request**.
The root is no longer captured in `main()`.

`--archive` keeps the pinned handler, deliberately: an archived GUI is a frozen
snapshot and a released run has to stay browsable after the topic moves on. The
two modes want opposite things, and the test asserts both.

## Why the wrong thing looked right

Because the fix that preceded it was real, and its docstring says so:

> *"the root is resolved from `run_paths` rather than typed, so the server
> cannot end up serving a superseded topic's pages — which it did once, for
> hours, from a literal path."*

That is true of the **literal**, and it is the sentence that stopped anyone
looking further. What it does not cover is the **lifetime**. Resolving from
`run_paths` once, at startup, is still a pin — it is just a pin whose value was
correct when it was taken. The class this project keeps rediscovering is
"a value that was right when written and cannot announce that it is not any
more" (disguise #3), and a startup-resolved root is exactly that, with a
half-life set by how long the process happens to live. These live for weeks.

The three servers that were serving the right topic were not evidence the design
worked: they were started after the last topic bump. Correct by luck, and the
luck runs out at the next bump rather than at a code change — so there is no
edit to blame it on and no moment at which anyone would think to look.

## Consequences

`:8951` now serves the live run and follows a topic change. The stale `:8931`
is still up and still serving nac_v5 — it needs restarting, and until it is,
**anyone following the documented `ssh -L 8931` instruction is reading nac_v5's
ranking under nac_v6's name.** That is the most likely way this bites someone
next.

Per-request resolution costs one cached-config deepcopy per request, on a
loopback server for a handful of readers. It is not a performance question.
