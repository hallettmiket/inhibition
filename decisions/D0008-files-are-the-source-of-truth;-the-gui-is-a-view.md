---
id: D0008
title: Files are the source of truth; the GUI is a view
date: 2026-07-27
status: accepted
approach: integration
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - decisions/
  - shared/decisions.py
  - integration/app/
evidence:
  - 'provenance was scattered across 7 places: commit messages, .provenance.md, prep_log.json, manifest.json, runbooks, config comments, ready_to_delete.md'
  - 'a published result must be reproducible from the repo alone, without running a web app'
  - 'the GUI is built last (M5); M0-M4 need provenance now'
runbook: null
---

## Context
Searching seven formats to answer 'why is the box 26 A' is untenable, and a single interface holding everything is genuinely more discoverable. But making the GUI the STORE would make reproduction depend on standing up Streamlit, and would lose git diff, code review, and bisect.

## Decision
One store, many views, with files as the store. Decisions live in decisions/ as git-versioned records with machine-readable frontmatter. The GUI aggregates decisions + manifests + runbooks into a single Decisions pane, per approach and for the choreography.

## Consequences
The GUI gains a real job beyond candidate display. If lab members ever need to AUTHOR decisions rather than read them, the GUI should write THROUGH to these files rather than into a database, so the repo stays authoritative.
