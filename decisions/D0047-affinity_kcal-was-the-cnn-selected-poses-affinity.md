---
id: D0047
title: affinity_kcal was the CNN-selected pose's affinity, not the best affinity — and a quarter of T_4 was ranked on clashing poses
date: 2026-07-31
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/covalent_protocol.py
  - approaches/t3_reinvent/04_rank.py
  - approaches/t4_combinatorial/04_rank_within_class.py
  - decisions/D0011-cnn-affinity-is-uncalibrated-for-covalent-docking.md
  - decisions/D0043-the-scoring-functions-rank-on-molecular-size.md
evidence:
  - "gnina's results table is ordered by CNN pose score; covalent_protocol.py read affinity off rows[0]"
  - 'rows[0] was NOT the affinity-best mode for 3631/4080 T_3 (89.0%) and 1503/1683 T_4 (89.3%)'
  - 'median affinity change on correction: T_3 -1.29, T_4 -1.57 kcal/mol'
  - 'POSITIVE (clashing) affinity under rows[0]: 69/4080 T_3, 423/1683 T_4 (25.1%), max +159.7'
  - 'after correction: T_3 has 0 positive affinities; T_4 has 181'
  - 'the residual 181 T_4 candidates have NO non-clashing pose in any of 9 modes — they do not fit'
  - 'shortlist churn T_3: 25 -> 25, 12 kept, 13 replaced'
  - 'shortlist churn T_4: 27 -> 27, 10 kept, 17 replaced (within-class quotas amplify a re-order)'
  - 'no re-docking: every pose already existed; only the choice of which to read was wrong'
  - 'D0046 context: the scoring function being corrected recovers known Pin1 poses 5% of the time in production'
---

# We rejected the CNN as a ranking signal and then let it choose the pose

## The defect

gnina returns nine modes and prints a results table **ordered by CNN pose
score**. `covalent_protocol.py` did:

```python
best = rows[0]
...
"affinity_kcal": best.get("affinity"),
```

So `affinity_kcal` meant *"the affinity of the pose gnina's CNN liked best"*,
not *"the best affinity"*.

That is incoherent with a decision already on the books. **D0011 demoted
`cnn_affinity` to advisory** because gnina itself reports that CNN scoring is
not calibrated for covalent docking, and T_3 and T_4 rank on `affinity_kcal`
instead. The project rejected the CNN as a *ranking* signal and then accepted
its *pose selection* without noticing that the two are the same choice wearing
different clothes.

## How wrong

| | T_3 | T_4 |
|---|---|---|
| candidates with poses | 4,080 | 1,683 |
| rows[0] **not** affinity-best | 3,631 (**89.0%**) | 1,503 (**89.3%**) |
| median penalty | 1.29 kcal/mol | 1.57 kcal/mol |

## The part that is worse than "suboptimal"

A positive Vina-style affinity is a **sterically clashing pose**. Under
`rows[0]`:

* T_3: **69** candidates carried a positive affinity, up to **+159.7**
* T_4: **423 of 1,683 — 25.1%** — carried a positive affinity **as their
  ranking value**

The CNN was preferring poses the affinity function scores as physically
impossible, and those scores were ranking the covalent shortlists.

After correction T_3 has **zero** positive affinities. **T_4 still has 181**,
and that is a different and real signal: those candidates have no non-clashing
pose in *any* of their nine modes. They do not fit the pocket. That should be
treated as a failure flag rather than a score, and it is recorded here as open
rather than quietly folded in.

## What changed downstream

| | old | kept | replaced |
|---|---|---|---|
| T_3 shortlist | 25 | 12 | **13** |
| T_4 shortlist | 27 | 10 | **17** |

T_4 churns harder than the raw re-order implies because it ranks *within
warhead class* with a fixed quota, so a re-order inside one class cascades.

**Consequences not yet addressed:**

* `shortlist_synth` is stale — it was built on the old ranking and must be
  rebuilt.
* MM-GBSA and explicit MD ran on the old shortlist members. 13 T_3 and 17 T_4
  of those are no longer shortlisted, and the newly promoted candidates have no
  physics.

## How it was fixed

`min(rows, key=affinity)` rather than `rows[0]`, and **all three scores now come
from the same mode**. Reading affinity from the affinity-best pose while leaving
`cnn_score`/`cnn_affinity` on the CNN-best pose would describe two different
geometries in one row — a subtler instance of the same defect. `selected_mode`
records which mode won, so the disagreement is visible rather than inferred.

`scripts/reextract_covalent_affinity.py` re-derives the columns from the
**existing** poses. No re-docking: every pose already existed and only the
choice of which to read was wrong. It preserves the previous values as
`*_rows0`, because a correction that erases the thing it corrected cannot be
audited, and it deliberately does **not** rewrite ranks — re-ranking is a
separate explicit step, so the script cannot silently reorder a shortlist as a
side effect.

## How it was found, which is the uncomfortable part

Not by auditing the ranking metric. It surfaced while building a **multi-pose
viewer for the GUI** (issue #3) — the viewer needed to show all nine modes, and
showing them made the disagreement visible.

Nobody was looking. The metric had a plausible name, produced plausible
negative numbers for 75% of rows, and had already survived the D0043 audit,
which corrected *which column* was analysed without ever asking which *pose*
the column came from.

This is the pattern in `docs/how_this_project_breaks.md` — a value taken by
**position** rather than by identity, failing silently. `rows[0]` is only
correct if the row order encodes the property you are reading, and here it
encoded a different score entirely.

## What this does NOT mean

The corrected shortlists are **not** more trustworthy in absolute terms. D0046
measured this same scoring function recovering known Pin1 poses **5% of the
time** on the production receptor. We have corrected which pose a barely-working
function is read from. The metric now means what its name says; it does not
follow that the ranking is now informative.

Both are true and both should be quoted together.
