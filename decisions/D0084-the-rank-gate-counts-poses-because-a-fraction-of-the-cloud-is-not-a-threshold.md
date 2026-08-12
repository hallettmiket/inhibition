---
id: D0084
title: The rank gate counts poses — a fraction of the cloud is a threshold that moves when the workload does
date: 2026-08-12
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/target.yaml
  - shared/target_config.py
  - scripts/rank_v2.py
  - tests/test_rank_gate.py
  - decisions/D0083-a-first-stage-mode-is-a-chain-not-a-cluster-and-the-letters-implied-otherwise.md
evidence:
  - '@tt8804: "that gate can be revisited without re-docking - it is the cheapest real improvement available right now"'
  - 'consensus == n_poses_mode / n_poses exactly, for all 34,059 rows (checked, 0 exceptions)'
  - 'every cloud on this run holds 500 poses, so `consensus >= 0.05` IS `n_poses_mode >= 25`'
  - 'the smallest mode passing the old gate holds exactly 25 poses'
  - 'nothing measured 25; 12 is measured (sweep_rule.min_mode_poses == splitting.stage2.min_mode_size)'
  - 'T_4 production ranking: gate admits 5,998 -> 8,334 modes (+2,336); ranked 5,132 -> 6,338 (+1,206); 0 lost'
  - 'among the 998 modes clearing enrichment 4.0: ranked 289 -> 434 (+145)'
  - 'of the strong modes the old gate excluded, 589 held 0-3 poses and 452 held 4-11 -- correctly excluded; the 157 holding 12-24 were not'
  - '1,996 modes now admitted remain unranked because conditional_eb is NaN -- not measured, so not ordered'
  - 'no re-dock, no re-screen: every enrichment, viable_fraction and pose is unchanged'
runbook: null
---

## Context

`rank_v2` decided which modes may hold a `class_rank` with `consensus >= 0.05`,
where `consensus` is the mode's share of its molecule's pose cloud. The stated
purpose was *"the mode is real"*.

`consensus` is **exactly** `n_poses_mode / n_poses` — checked against all 34,059
rows, no exceptions — and every cloud on this run holds exactly 500 poses. So
the gate was, precisely, **`n_poses_mode >= 25`**. The smallest mode that passed
it holds 25 poses.

Nothing measured 25. It is what 0.05 of 500 comes to.

## The defect is the division, not the number

A size floor written as a fraction of the workload is a threshold that **moves
when the workload moves**. Raise `docking.n_runs` from 500 to 1,000 — a change
nobody would think of as touching the ranking — and the same `0.05` silently
becomes "at least 50 poses", roughly halving what is ranked, with nothing in any
output saying the gate changed. Lower it and the ranking fills with two-pose
modes scoring the arithmetic maximum.

That is `how_this_project_breaks.md` disguise #3 — *a constant sized against a
workload that has since grown, which cannot announce that it is out of date* —
in the one place that decides what the shortlist is drawn from. It is the same
shape as the `timeout=86400` that killed a 24-hour run after the pool grew from
1,882 to 16,806.

`config/target.yaml` had already written down half of this when the sweep depth
was measured: *"NOTE this is NOT the `consensus >= 0.05` gate in rank_v2, which
is a fraction OF THE CLOUD and therefore conflates a small mode with a large
cloud."* The sweep was moved onto absolute size then. The ranking was not.

## Decision

**The rank gate counts poses.** `ranking.mode_gate` in `config/target.yaml`:
`n_poses_mode >= 12`, read through `target_config.rank_min_mode_poses()`, which
refuses rather than defaulting.

**12 is the number already measured on this target.** It is
`sweep_rule.min_mode_poses`, itself tied to `splitting.stage2.min_mode_size` —
the size below which a cluster is not treated as a mode at all. The three are
kept as separate keys so one can be changed without silently dragging the
others, and a test asserts they agree until someone decides otherwise.

**The old rule stays reachable**: `rank_v2 --gate consensus_fraction --floor
0.05` reproduces 3.0.0's earlier tables. Which gate ran is stamped on every row
(`rank_gate`, `rank_gate_min`), so a file cannot be read under the wrong rule.

**The gate refuses a frame it cannot evaluate.** A 2.1.0 aggregate has no
`n_poses_mode`; `rank_v2` exits with the flag to use rather than falling back to
the fraction. A silent substitution here is D0080's defect exactly — two topics
ranked by different rules under one filename.

## What it changed

T_4, the production ranking, `conditional_eb`:

| | old (`consensus >= 0.05`) | new (`>= 12 poses`) |
|---|---:|---:|
| modes admitted by the gate | 5,998 | **8,334** |
| of those, actually ranked | 5,132 | **6,338** |
| lost a rank | — | **0** |
| ranked, among the 998 clearing enrichment 4.0 | 289 | **434** |

**Nothing was re-docked and no score moved.** Every `enrichment`,
`viable_fraction` and pose is the value it was; only which modes are allowed to
hold a position changed.

The gain is specifically the **12–24 pose band**. Of the strong modes the old
gate excluded, 589 held 0–3 poses and 452 held 4–11 — those remain excluded, and
correctly: a two-pose mode scores the arithmetic ceiling of 12.25 whenever both
poses are viable, which is noise wearing a maximum. Only the 157 holding 12–24
were being thrown away for having been measured on a cloud of 500 rather than
for anything about the mode.

**1,996 admitted modes are still unranked**, because `conditional_eb` is NaN —
no pose reached the distance window, so the quantity was never measured. That is
the right treatment and it is a different statement from "failed": not measured
is not the same as measured and found absent.

## Why it looked right

`consensus` is a genuinely meaningful quantity — the share of a molecule's poses
that chose this mode — and gating on it reads as "is this mode the molecule's
real answer or a corner of the cloud". That question is worth asking, and it is
*still* asked: `conditional_eb` multiplies the score by consensus, so weight of
evidence remains part of the ranking rather than a threshold bolted beside it.

What the threshold form added was a hidden dependence on `n_runs`. The gate
looked like a statement about modes and was a statement about modes **divided by
a screening parameter**, and the two are indistinguishable as long as that
parameter never changes. It has not changed yet. D0083's sub-splitting is what
made it visible, by pushing many real modes just under it.

## Consequences

* **The shortlist is drawn from a larger, honestly-gated pool.** 1,206 more
  modes carry a rank, 145 of them clearing the enrichment floor of 4.0.
* **The running sweep is not affected.** Its worklist (`sweep_gaps_6.csv`,
  13:46) is a snapshot taken under the old gate, and the workers hold it;
  `attack_sweep` reads poses from `run.topic` and never reads a rank file.
  **Regenerating the worklist would change what a restarted worker does, and is
  deliberately NOT done here** — that is a scheduling decision for whoever owns
  the run.
* **Class ranks moved because the denominators grew**, not because any mode got
  better. Sulfopin's chloroacetamide class went from 545 ranked modes to 613, so
  its position reads bottom 9.3% rather than bottom 7.5%. D0082's conclusion is
  untouched — it never depended on the exact percentile.
* **`consensus` is still written and still scored on.** Nothing was removed; a
  threshold was replaced by the measured one.
* **This does not fix D0083.** The clustering still chains, and sub-splitting
  still divides consensus among the pieces. This removes the arbitrary part of
  what that cost, at zero compute. #65 keeps the clustering question open.
