---
id: D0106
title: The screen samples to a fixed number of VALID poses, because an equal denominator over an unequal eligible set is not an equal comparison
date: 2026-09-02
status: accepted
approach: shared
decided_by: '@twu383'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/target.yaml
  - scripts/nac_screen_v2.py
evidence:
  - '`consensus` = mode_size / n_poses, and `n_poses` was `len(res)` = the number DOCKED, equal at 500 for every molecule in nac_v6'
  - 'but `labels` is -1 for PoseBusters-invalid poses, so only valid poses can join a mode -- the numerator was drawn from a smaller and UNEQUAL set'
  - 'MEASURED over nac_v6''s 561 acrylamide/bdhi molecules: pass rate 0.812 to 0.982, median 0.920, 1st pct 0.834'
  - 'so a molecule at 0.812 valid had a consensus CEILING of 0.812 while one at 0.982 could reach 0.982 -- a 21% difference in achievable headroom, on the axis the ranking reads'
  - 'n_runs 640 chosen from that distribution: 500/0.812 = 616 is the worst case observed, 640 clears it (0.812 x 640 = 520 valid)'
  - 'smoke test on t4_80fbed3bdf1e: 640 docked, 578 valid, 500 kept, 208 modes, consensus sums to 1.0000 across modes, n_poses_mode sums to 500'
  - 'per-pose table retains all 640 rows; 78 poses are valid-but-not-kept, a third state that `pb_valid` alone cannot express'
  - 'exp/21: filtering a cloud to its best-scoring 25% concentrates attack-ready poses 2.60x over a random 25% (21 of 21 molecules, p = 6e-05) -- which is why the truncation is NOT by energy'
runbook: null
---

# Sample to 500 valid poses, not 500 attempts

## The defect this closes

`consensus` is a mode's share of the pose cloud, and every downstream quantity
that depends on mode population depends on it. It was computed as

    consensus = mode_size / n_poses,   n_poses = len(res) = docking.n_runs

which is 500 for every molecule. That reads as an equal denominator and the
config said so in as many words: *"`consensus` is exactly mode_size / n_poses
(checked: true for all 34,059 rows), and every cloud on this run holds 500
poses."*

**The denominator was equal and the numerator's ceiling was not.** `labels` is
`-1` for every PoseBusters-invalid pose, so only valid poses can be counted into
a mode. Measured over the 561 acrylamide/bdhi molecules of nac_v6, the valid
fraction runs from **0.812 to 0.982**. A molecule whose cloud was 81% valid could
not produce a consensus above 0.81 however tightly its poses agreed; one at 98%
could reach 0.98. Two molecules were being ordered on a quantity with 21%
different room to grow, and nothing in any artefact said so.

## Why it looked right

Because the fix that introduced it was itself correct, and it was checked.

D0089 adopted PoseBusters and chose to **flag rather than delete** — the right
call, with D0093 as the cautionary record of what a deleting filter costs. The
denominator was then verified to be equal across all 34,059 rows, and it was.
Nobody asked the next question: equal *out of what population*, and is the
population a mode can actually be drawn from the same one?

This is the catalogue's disguise #4 in its subtlest form. The guard ran, the
verification ran, the number it reported was true, and it answered a slightly
different question than the one its result was used for — the same shape as #27,
where a real calibration was measured for the wrong quantity.

## The decision

`docking.target_pb_valid: 500`. Dock `n_runs: 640`, keep the **first 500 valid
poses in docking order**, and compute every aggregate over those 500.

**640 is derived, not padded.** From the pass-rate distribution of the same 561
molecules this run screens:

| quantile | pass rate | n_runs for 500 valid |
|---|---|---|
| median | 0.920 | 544 |
| 5th pct | 0.864 | 579 |
| 1st pct | 0.834 | 600 |
| **minimum** | **0.812** | **616** |

At the worst rate observed anywhere in that set, 640 runs yield 520 valid. Well
below the ~2,000 where AutoDock-GPU begins failing silently (#77).

**The FIRST 500, in docking order, and the reason is the whole argument.**
AutoDock-GPU's runs are independent GA replicates, so the first N valid poses
are an unbiased sample of the valid population and the cloud keeps meaning *what
docking produces*. Position-based selection is defensible here and **only**
because the runs are i.i.d. — which is exactly the question disguise #2 says to
ask of any index ("if you index, write down what guarantees the order").

Truncating by **energy** was rejected. It would make the cloud "the best-scoring
500", a different population: exp/21 measured that a cloud's best 25%
concentrates attack-ready poses **2.60x** over a random 25% (21 of 21 molecules,
p = 6e-05). `engagement` would rise for reasons that are not chemistry, and the
number would not be comparable with nac_v6's.

## Nothing is deleted, and the record gains a third state

Every one of the 640 docked poses keeps its per-pose row and its place in
`<topic>_allposes`. `pb_valid` says whether PoseBusters passed it. **`pb_kept`
says whether it is in the analysed 500** — because "valid" and "analysed" are
now different, and a column that cannot express the difference would hide it.
Verified on the smoke test: 640 rows, 578 valid, 500 kept, 78 valid-but-not-kept,
0 kept-but-invalid.

`n_poses_docked`, `n_poses_pb_valid`, `n_poses_kept`, `target_pb_valid` and
`target_met` all travel on every aggregate row. A molecule that cannot reach the
target is **logged as a warning and flagged `target_met = False`** rather than
quietly carried: it is not comparable with the rest of the run, and a table that
mixes the two silently reintroduces this defect.

## What this invalidates

**nac_v6's `consensus` and everything derived from it are not comparable with
nac_v7's.** The run moves to its own topic for that reason. `target_pb_valid: 0`
restores the previous behaviour exactly, so nac_v6 remains reproducible rather
than only described.

Note that `engagement` — the score `ranking.score_by_tier` names for T_4 — is a
pose property and does not depend on the denominator, so it is affected only
through which poses are in the analysed set. `conditional_eb` and anything built
on `consensus` are affected directly.

## What is still not fixed

The mode COUNT still grows with sampling depth (D0092, n^0.69 with no plateau),
and this run docks 640 where nac_v6 docked 500. The analysed cloud is held at 500
so the count is comparable at fixed analysed depth, but the *deeper* docking
means the 500 kept poses are drawn from a wider sample. Whether that changes the
group structure is measurable against nac_v6 on the shared molecules and has not
been measured yet.
