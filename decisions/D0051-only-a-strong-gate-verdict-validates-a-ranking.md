---
id: D0051
title: Only a STRONG gate verdict validates a ranking — the check is an allowlist, not a denylist
date: 2026-08-01
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/rank_shortlist.py
  - tests/test_rank_within_class.py
evidence:
  - 'the check was `verdict not in (UNDERPOWERED, UNGATED, FAIL)` — permissive by default for any unanticipated verdict'
  - 'the non-covalent gate verdict moved UNDERPOWERED -> WEAK during 2026-08-01'
  - 'WEAK is defined in enrichment_gate.py as "enriches, but within noise of not enriching"'
  - 'the WEAK verdict carries ROC-AUC 0.599, CI [0.311, 0.874], EF1% 0.0 — the D0041 numbers'
  - 'D2_26 (written 09:17, before any change today) already carried gate_verdict=WEAK with rank_validated=True'
  - 'D1_25/D2_27 reproduced it; D3_21/D4_33 were unaffected because the covalent gate is still UNDERPOWERED'
  - 'graded vocabulary is STRONG | WEAK | UNDERPOWERED | FAIL (D0012), plus UNGATED for a missing token'
---

# A ranking is validated only when the gate says STRONG

## Context

`attach_gate` decided `rank_validated` by exclusion:

```python
out["rank_validated"] = str(g.get("verdict", "UNGATED")) not in (
    "UNDERPOWERED", "UNGATED", "FAIL")
```

Every verdict the list did not name validated the ranking. That was invisible
while the only verdicts ever produced were on the list.

On 2026-08-01 the non-covalent enrichment gate moved from `UNDERPOWERED` to
`WEAK`. `WEAK` is not a pass — `enrichment_gate.py` defines it as *"enriches,
but within noise of not enriching"*, and the verdict carries ROC-AUC 0.599 with
a confidence interval of [0.311, 0.874] and EF1% of 0.0, which are exactly the
D0041 numbers that established docking does not demonstrably enrich here.

So T_1 and T_2 began shipping frames stamped **`rank_validated = True`**,
contradicting D0041, `docs/state_of_the_project.md`, the README, and
`rank_shortlist`'s own docstring, which says in as many words that nothing here
currently clears the gate. `D2_26` carried it before any of today's changes.

## Decision

**`rank_validated` is true only for `STRONG`.** The check names what validates
rather than what does not:

```python
VALIDATING_VERDICTS = {"STRONG"}
```

An unrecognised verdict logs a warning and does **not** validate. The graded
vocabulary is listed as `GATE_VERDICTS` so a new grade is noticed rather than
silently bucketed.

## Why this is the catalogued failure and not a typo

`docs/how_this_project_breaks.md`, disguise 4: *a guard that is scoped out,
mis-ordered, or vacuous.* The defence it prescribes is to ask **what would make
this pass when it should fail** — and here the answer was "any verdict string
the author did not think of", including every verdict that might be added
later. A denylist cannot fail closed. That is the class; `WEAK` was merely the
first instance to arrive.

Note also how it surfaced: not by a test and not by an audit, but by reading
the log line of a re-rank run for an unrelated reason. That is route 1 in the
catalogue's own scorecard ("someone looked at output and it didn't match
expectation"), which remains the way most of these are found.

## Consequences

The four arms were re-ranked. T_1 and T_2 return to `rank_validated = False`;
T_3 and T_4 never changed, because the covalent gate is still `UNDERPOWERED`
and that was already on the denylist.

`gate_verdict` is unchanged and still carried on every row — the *grade* is
information, and D0012's position is that the gate reports evidence strength
rather than adjudicating. What changes is only the boolean that claims the
ranking is validated.

**Only `FAIL` demotes `dock_score` to a displayed label** (D0012), and that is
untouched: a WEAK ranking still ranks, still shortlists, and still carries its
verdict and interval. This decision does not make the pipeline more
conservative about *using* the score — only about *claiming* the score has been
validated.
