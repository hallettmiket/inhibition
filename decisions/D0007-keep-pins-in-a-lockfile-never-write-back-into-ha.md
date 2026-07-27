---
id: D0007
title: Keep pins in a lockfile, never write back into hand-authored config
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - config/sources.yaml
  - config/sources.lock.json
  - shared/sources.py
evidence:
  - 'the first sources.py round-tripped sources.yaml through yaml.safe_dump'
  - 'this erased all 22 explanatory comments in the file'
  - 'the file had not been committed, so git could not restore it'
  - 'pin enforcement verified by tampering: mismatch raises SourceError'
runbook: null
---

## Context
Acquisition needs to record observed hashes so they can be enforced on later runs. The obvious implementation - write them back into the source config - destroys anything a human wrote there.

## Decision
config/sources.yaml is hand-authored and never written by code. Observed hashes and resolved commits go to config/sources.lock.json. Config is written by humans; pins are written by code; the two never share a file.

## Consequences
Standard dependency-lockfile pattern. Adding a source means editing YAML; pinning happens automatically on first acquisition and is enforced thereafter.
