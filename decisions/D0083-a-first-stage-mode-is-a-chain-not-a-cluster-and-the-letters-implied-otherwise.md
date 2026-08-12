---
id: D0083
title: A first-stage mode is a chain, not a cluster — the lettered sub-mode labels claimed a similarity the geometry does not have
date: 2026-08-12
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/mode_ranking.py
  - shared/pose_modes.py
  - shared/pose_subsplit.py
  - scripts/mode_homogeneity.py
  - decisions/D0079-pose-splitting-and-per-mode-ranking-are-accepted.md
evidence:
  - '@tt8804 on the ranking view: "how are the sub modes so diff from eachother, these should be individual modes not submodes ... they should be numbered individually not sub"'
  - '7,968 first-stage clusters were sub-split on this run'
  - 'spread of sub-mode median distance-to-sulfur WITHIN one cluster: median 0.63 A, 90th pct 2.25 A, max 7.68 A'
  - '1,763 of 7,968 (22.1%) span MORE than the criterion window (2.8-4.2 A) is wide'
  - '1,445 of 7,968 (18.1%) hold sub-modes typically INSIDE the window and sub-modes typically OUTSIDE it'
  - 'worst case t4_152e517aa867 cluster 0: sub-modes at 3.17 A and 9.73 A, enrichment 0.00 to 12.25 (the arithmetic ceiling) in ONE first-stage mode'
  - 't4_e0b03662d460 cluster 0: sub-mode median distances 3.57 / 3.82 / 5.93 / 5.96 / 6.83 A, enrichment 0.00 to 2.66'
  - 'stage 1 is DBSCAN, eps 3.0, on positional_separation + 2.0 * angular_difference_rad: co-located poses are neighbours up to 86 deg apart, and DBSCAN chains'
  - 'the identity was already flat: t4_e0b03662d460_m0..._m8, nine modes; the letters existed only in the GUI'
runbook: null
---

## Context

The ranking view grouped a molecule's modes under the first-stage mode they were
cut from, named them `0a`…`0e` / `1a`…`1d`, gave a group one hue with lightness
steps, and told the reader they were "sub-splits of ONE first-stage mode".

@tt8804, looking at `t4_e0b03662d460`: *"this also makes no sense?? how are the
sub modes so diff from eachother, these should be individual modes not
submodes."*

The observation has force because of what the first stage clusters on. It groups
poses by the reactive atom's **position** and the **direction** the warhead
faces (`pose_modes.features`) — and the criterion scores a pose on the reactive
atom's **distance** to Cys113 SG and its **approach angle**. Those are the same
two quantities. So sub-modes of one first-stage mode ought to score alike, and
should differ only away from the reactive end. That is precisely the claim #61
made when the second stage was introduced: *"two poses that place the warhead
identically and hang the scaffold differently."*

Measured, that claim fails for a large minority of the library
(`scripts/mode_homogeneity.py`):

| | |
|---|---:|
| first-stage clusters sub-split | 7,968 |
| median spread of sub-mode median distance-to-sulfur | 0.63 Å |
| 90th percentile | 2.25 Å |
| **spanning more than the whole 1.4 Å criterion window** | **1,763 (22.1%)** |
| **holding sub-modes both inside and outside the window** | **1,445 (18.1%)** |
| worst | 7.68 Å |

`t4_152e517aa867`'s first-stage mode 0 contains sub-modes at 3.17 Å and 9.73 Å,
scoring enrichment 0.00 and 12.25 — the floor and the arithmetic ceiling, inside
one "mode".

## Why it happens

Stage 1 is DBSCAN at `eps = 3.0` over
`positional_separation + ANGSTROM_PER_RADIAN * angular_difference`, with
`ANGSTROM_PER_RADIAN = 2.0`. Two poses in the same spot are neighbours up to
1.5 rad — **86°** — apart. And DBSCAN is a **chaining** method: A joins B, B
joins C, and A and C share a mode having never been within `eps` of one another.
A first-stage mode is a **connected component, not a ball**, and its diameter is
unbounded.

The second stage then cuts the chain apart at a 2 Å whole-molecule diameter. It
separates poses by reactive-atom distance as a *side effect*, because that
distance is part of whole-molecule geometry. **The second stage has been doing
the first stage's job**, which is why sub-modes of one cluster can score
anything.

## Decision

**Two things, and only the first is done here.**

1. **The ranking view names every mode by its own index** — `m0`…`m8`, which is
   its identity in the rank table, the sweep table and the pose file. The
   lettered names are gone. The first-stage cluster is still shown, as
   provenance in the group header and the detail caption, and still drives the
   shared hue so a cluster can be drawn at once — with the header now saying
   these are separate modes and the note warning that a shared hue does not mean
   the rows are alike.

   **This is display-only and costs nothing.** The identity was already flat:
   `t4_e0b03662d460_m0` … `_m8`, nine modes, each with its own `class_rank`, its
   own enrichment and its own slot in the sweep worklist. Every number on the
   screen was already computed treating them as independent. Only the label said
   otherwise, so no frame, no score and no join changes.

2. **The clustering itself is NOT changed here.** Tightening stage 1 moves every
   mode boundary in the library, which moves `consensus` (mode_size / n_poses)
   and everything derived from it — a full re-screen and another MAJOR version
   by the CHANGELOG's own rule. That is a decision to take deliberately, with
   the pilot, not as a side effect of fixing a label. Tracked in #65.

## Why it looked right

The lettering was not careless; it was a considered fix to a real problem, and
its own rationale is still in the code it replaced. Sub-modes had been coloured
by their renumbered index, so one cluster came out in five unrelated hues and
read as five unrelated modes. The fix asserted the opposite relationship — one
hue, lightness steps, lettered names — and that assertion was **true of the
design and false of the data**. The design intends stage 1 to bound warhead
geometry; nothing measured whether it does.

This is the catalogue's shape at the level of a claim rather than a value: a
label describing what the code was *meant* to produce, sitting beside numbers
computed from what it *did* produce, with both populated and plausible. It took
a reader noticing that two rows in one group scored 0.00 and 2.66.

## Consequences

* **The ranking view no longer implies a similarity it cannot support.** A
  reader comparing `m3` and `m4` of one molecule now sees two modes, which is
  what every downstream stage already treated them as.
* **`mode_label` stays on the frames.** It is still written by
  `nac_screen_v2` and still records which first-stage cluster a mode came from —
  useful provenance, and removing it would invalidate frames for a naming
  preference. It is simply no longer used as a mode's name.
* **The 22% figure bounds how much of the ranking is affected by chained
  clusters**, and `scripts/mode_homogeneity.py` re-measures it on any run, so a
  future change to `eps` can be judged rather than argued.
* **The consensus gate compounds this.** Splitting divides `consensus` among the
  pieces, so a chained cluster that splits into five leaves five small consensus
  values; 1,447 of the 1,772 modes clearing enrichment 4.0 are then cut by the
  `consensus >= 0.05` gate. That gate was already flagged as suspect when the
  sweep depth was measured — its median had landed on its own threshold — and
  this is a second reason to revisit it. It can be revisited **without**
  re-docking, unlike the clustering.
* **If stage 1 is tightened, the second stage's measured benefit must be
  re-derived.** The k=1→4 recovery gain (22.0% → 39.0%, p = 1.2×10⁻⁴) was
  measured against the *current* first stage. Part of that gain is the second
  stage repairing chained clusters, and a tighter first stage would leave it
  less to repair.
