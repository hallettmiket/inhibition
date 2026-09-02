---
id: D0104
title: The Ranking page displayed conditional_eb while sorting on engagement, so it contradicted itself on a third of its rows
date: 2026-08-31
status: accepted
approach: integration
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/mode_ranking.py
  - tests/test_ranking_page_shows_its_own_score.py
evidence:
  - '@tt8804, reading modes.html: "what is that number and why does it seem like they are not being ranked by it, the 2.94 over the 4."'
  - 'the rail printed `fmt(x.eb, 2)` where `eb` = `conditional_eb`, named as a literal in _rows_json; the list sorts on class_rank, which rank_v2 computes from ranking.score_by_tier (engagement for T_4, D0098)'
  - 'MEASURED on nac_v6 acrylamide (34,888 modes): rho(engagement, conditional_eb) = 0.320'
  - 'the displayed number rose as the rank fell at 9,077 of 24,927 adjacent pairs (36.4%)'
  - 'the screenshotted rows: class_rank 7-18 have engagement 0.9129, 0.9114, 0.9080, 0.9041, 0.9023, 0.9017, 0.9010, 0.9004, 0.8932, 0.8848, 0.8833, 0.8759 -- STRICTLY DESCENDING, i.e. the ordering was correct throughout'
  - 'their conditional_eb: 2.935 x6, then 4.315 (rank 13), then 2.483 (rank 14), then 2.935 x4'
  - 'conditional_eb collapses to its prior at small n: 15,988 of 34,888 acrylamide modes hold ONE pose and take just 2 distinct conditional_eb values, while their engagement spans 0.0000 to 0.9557'
  - 'the 2.94 rather than 2.93 is DOUBLE ROUNDING: 2.934896 -> round(.,3) = 2.935 in Python -> toFixed(2) = "2.94" in JS'
  - 'SECOND DEFECT: global_rank was computed from conditional_eb while class_rank comes from engagement, so the scope selector''s two settings ordered by different quantities'
  - 'AFTER THE FIX: 0 of 34,887 adjacent pairs contradict the order'
runbook: null
---

## Context

@tt8804 looked at the Ranking page and asked why rows 7–12 and 15–18 all showed
**2.94** while row 13 showed **4.32** — the list plainly not ordered by the
number printed on it.

**The ordering was correct.** Ranks 7–18 have `engagement` 0.9129 → 0.8759,
strictly descending. The number on screen was a different column.

## Decision

The rail displays the score the list is ordered by. `gather()` stamps `_score`
and `_score_col` onto every row from the column `ranking.score_by_tier` names
for that tier, and the payload sends it as `sc`. `conditional_eb` stays in the
detail panel as context, where it is genuinely useful, and is no longer the
headline. The score's NAME is rendered with it, so the page states what it is
showing.

`global_rank` now ranks the same `_score`. It ranked `conditional_eb` while
`class_rank` arrives from `rank_v2` computed on `engagement`, so switching the
scope selector from "within class" to "global" silently changed the score as
well as the comparison set.

`gather()` raises if the configured score is absent from the frame rather than
falling back to another column — falling back is the defect.

## Why the wrong thing looked right

**Because `conditional_eb` is a perfectly good number.** It is a real quantity,
computed correctly, in a plausible range, formatted to two decimals in the slot
where a score belongs. Nothing about it looks like a bug. This is disguise #1 in
[`how_this_project_breaks.md`](../docs/how_this_project_breaks.md) — two columns
exist, both populated, both plausible, and the code names one.

**Because the same fix had already been made one step away.** `gather()` reads
the score name from config, with a comment explaining that a hardcoded
`conditional_eb` had made all 327,167 of nac_v6's modes invisible to the GUI —
*"which reads as 'the ranking has not run' rather than 'this reader is looking
for a different filename'."* That fix was applied to choosing the FILE and not
to the number the reader looks at. The literal survived in the place with the
most eyes on it.

**And the symptom pointed away from the cause.** `conditional_eb` is an
empirical-Bayes estimate that shrinks to its prior at small n, and 46% of these
modes hold a single pose — so the displayed column takes two values across
16,000 rows while `engagement` spans the full 0–0.96 underneath. The page
therefore looked like it *had no ordering*, which invites "the ranking is
broken" rather than "the wrong column is printed". The ranking was fine; the
label was lying about it. That is the more dangerous failure of the two, because
a broken ranking gets fixed and a mistrusted correct one gets abandoned.

A smaller wrongness rides along: the printed **2.94** is not even
`conditional_eb` to two places. `_rows_json` rounds to 3 (2.935), then the
client's `toFixed(2)` rounds again to 2.94, where the true value 2.934896 is
2.93. Rounding twice is its own small instance of a value that stopped being
what it says it is.

## Consequences

0 of 34,887 adjacent acrylamide pairs now contradict the ordering, against 9,077
of 24,927 before.

Nothing about the SCIENCE changed — no number moved, no molecule reordered, and
the ranking this page shows is the same ranking it was showing. What changed is
that the page now agrees with itself. Every ranking read off `modes.html` before
2026-08-31 was ordered correctly; anyone who compared the printed figures
against that order and concluded the ranking was broken was reading a display
bug.
