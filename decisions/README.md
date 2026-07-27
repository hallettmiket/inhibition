# Decision records

One file per consequential decision. Git-versioned, reviewable in a PR,
machine-readable, and aggregated into the GUI's **Decisions** tab so nobody has
to grep seven places to learn why the docking box is 26 Å.

**Files are the source of truth; the GUI is a view over them.** A published
result must be reproducible by reading this repo, without standing up a web
app. The GUI reads these records — it does not own them.

## When to write one

Write a record when a choice would be **expensive or confusing to reverse**, or
when someone reading a number six months from now would reasonably ask "why?":

- a shared artifact everything depends on (the receptor, the boxes, the seeds);
- rejecting an available option (not reusing a prior library, not going wide);
- a control's interpretation (what the reactivity window is *for*);
- anything where the spec and reality diverged (a stale install note, a
  requirement that turned out unsatisfiable).

Do **not** write one for routine implementation choices. A record that says
"used pandas" is noise, and noise is what makes people stop reading.

## Format

YAML frontmatter (machine-readable, drives the GUI) then prose:

```yaml
---
id: D0007                      # sequential, never reused
title: short imperative phrase
date: 2026-07-27
status: accepted               # proposed | accepted | superseded | rejected
approach: shared               # shared | t1 | t2 | t3 | t4 | integration
decided_by: '@mhallet'
origin: spec                   # spec | adversary | implementation | user
supersedes: []                 # ids this replaces
superseded_by: null            # id that replaced this
affects:                       # files/artifacts this governs
  - config/receptor.yaml
evidence:                      # what backs it — numbers, not adjectives
  - 'LINK record: SG CYS A 113 - C10 QT7 A 201, 1.78 A'
runbook: docs/runbooks/receptor_selection.md   # or null
---

## Context
What forced a choice. Include the constraint, not just the topic.

## Decision
What was decided, stated so it can be checked.

## Consequences
What follows — including what is now harder, and what would need to change if
this is revisited.
```

## Superseding

Never edit an accepted record to change its meaning. Write a new one, set its
`supersedes`, and set the old one's `status: superseded` and `superseded_by`.
The history of *why the answer changed* is usually more informative than the
current answer.

Correcting a typo or adding evidence to an existing record is fine.

## Validation

`shared/decisions.py` validates the frontmatter and refuses duplicate ids,
dangling `supersedes` references, and unknown status values:

```bash
python -m shared.decisions check
python -m shared.decisions list --approach t4
```

## Relationship to the other provenance layers

| Layer | Answers |
|---|---|
| **decision record** | *why* we chose this |
| [runbook](../docs/runbooks/) | *how* to make this kind of choice again |
| `manifest.json` | *what* a specific run actually consumed |
| `.provenance.md` | where the reference data came from |

A decision that keeps recurring should get a runbook. A runbook applied to a
specific case should produce a decision record.
