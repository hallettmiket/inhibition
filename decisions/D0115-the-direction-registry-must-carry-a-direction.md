---
id: D0115
title: The rank-direction registry must carry a direction, not just a name
date: 2026-09-09
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/rank_shortlist.py
  - shared/pharmacophore.py
  - tests/test_rank_direction.py
evidence:
  - 'LOWER_IS_BETTER was a set of 3 names; rank() checked membership then hardcoded ascending=True at 4 sites'
  - 'the ValueError said "add it to LOWER_IS_BETTER with its direction" — the structure had nowhere to put one'
  - 'simulated on the old code path: a higher-is-better metric ranks similarity 0.10 first and 0.30 last'
  - 'ph2d_max_similarity is the first higher-is-better rank metric in the project'
runbook: null
---

## Context

`shared/rank_shortlist.rank()` refuses any metric not named in
`LOWER_IS_BETTER`. That guard is the one this project points to as a success —
it caught catalogue #4, an analysis ranking on `cnn_affinity` when the ranking
column was `affinity_kcal`, and `how_this_project_breaks.md` records it as one
of only three defects a guard ever caught.

It was a bare **set of names**:

```python
LOWER_IS_BETTER = {"affinity_kcal", "vina_affinity", "size_decorrelated_score"}
```

All three members are kcal/mol, so all three are lower-is-better, and every
sort underneath simply hardcoded `ascending=True` — four sites, plus two
`.min()` aggregations for the identity-collapsed path. The direction was never
read; it was assumed, and the assumption happened to be right for every metric
that existed.

Adding a pharmacophore similarity forced the issue. It is a Tanimoto: **higher
is better**.

## Why it looked right

Because the guard *worked*, and had a track record. It refused an unregistered
metric, loudly, with a message naming the fix.

The message was the trap:

> `"...is not a known rank metric; add it to LOWER_IS_BETTER with its direction
> rather than assuming one"`

There is no way to add a direction to a set. The only action the message
admits is `LOWER_IS_BETTER.add("ph2d_max_similarity")` — which registers a
higher-is-better metric in a container whose name asserts the opposite, and
every sort then runs `ascending=True` against it. Simulated on the old code
path: the molecule with similarity **0.10 ranks first** and the one at **0.30
ranks last**. The shortlist is the 25 molecules *least* like a known Pin1
binder, every value populated, every value plausible, the ordering exactly
inverted, and nothing raises.

The guard checked the **name** and assumed the **direction** — which is the
project's signature defect one level up from the defect it was built to catch.
It sat this way for six weeks because a set with three same-signed members
cannot express the bug.

The alternative workaround is worse and was rejected: storing `-similarity` so
that lower is better. That produces a column whose name says one thing and
whose sign means another, which is disguise #1 with the evidence hidden inside
an arithmetic operation instead of a column name.

## Decision

The registry carries the direction per metric and every sort **reads** it:

```python
RANK_DIRECTION = {
    "affinity_kcal": True,              # kcal/mol, lower binds better
    "vina_affinity": True,
    "size_decorrelated_score": True,    # residual of a lower-is-better metric
    "ph2d_max_similarity": False,       # a Tanimoto, higher is better
}
LOWER_IS_BETTER = {m for m, lower in RANK_DIRECTION.items() if lower}
```

`LOWER_IS_BETTER` is kept as a **derived** view so existing readers
(`composite_rank`, `pose_vector`, `test_composite_rank`) keep working and the
two cannot drift. `direction_of(metric)` raises on an unregistered name — the
original guard is unchanged in strength.

**The residual inherits the direction of its base metric.**
`size_decorrelated_score` is registered lower-is-better because it has only
ever been a residual of kcal/mol, but `rank()` resolves the direction from the
metric it was *asked* for, not from the residual's own entry. Otherwise
enabling decorrelation on a higher-is-better metric would invert the ranking —
the same bug, reachable through a keyword argument.

## What was checked

`tests/test_rank_direction.py`, seven tests. The load-bearing one is
`test_higher_is_better_is_not_ranked_backwards`, and it was **verified to fail
against the old behaviour** rather than assumed to: forcing `direction_of` to
return `True` (what the old hardcoded path did) makes the least-similar
molecule rank first. Two of this project's own tests have passed for free by
having nothing to check (`how_this_project_breaks.md`, disguise #4), so a new
guard is not trusted until it has been seen to fail.

`test_the_two_directions_produce_opposite_orderings` fixes the mechanism rather
than today's data: both fixture metrics order the molecules identically by
value, so a direction-blind implementation would return the same ranking for
both, and the test asserts the raw orderings really are opposite so the
agreement is the direction doing work rather than a coincidence.

## What else shares this shape

Worth a sweep, in the spirit of catalogue #9 (caught *before* it bit by asking
what a cache key omits right after #8 taught the question):

- **Any allowlist whose entries carry an implicit shared property.** This
  registry's members shared a sign. `GATE_VERDICTS` and `VALIDATING_VERDICTS`
  are allowlists of bare strings — they are fine today because the property
  being tested *is* membership, but the pattern is the one to watch.
- **`shared/engagement_rank.LOWER_IS_BETTER = False`** is a module-level
  direction sitting outside this registry, and `pose_consensus` documents a
  higher-is-better metric that is "NOT A RANK METRIC YET". Both are directions
  declared somewhere other than the place `rank()` looks. Neither is wrong
  today; both would be found by the same question.
- **`pose_vector`'s docstring claims its metric is "registered in
  `rank_shortlist.LOWER_IS_BETTER` before anything sorts on it".** It is not in
  the registry. Nothing sorts on it yet, so nothing is broken — but the
  docstring describes a state that does not hold, which is how a reader
  concludes the guard covers a path it does not.
