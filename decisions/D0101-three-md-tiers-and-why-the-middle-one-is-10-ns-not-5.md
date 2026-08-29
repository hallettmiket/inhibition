---
id: D0101
title: Three MD tiers — and the middle one is 10 ns, not 5, because the bar is nested and culling harder is free
date: 2026-08-29
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/target.yaml
  - shared/target_config.py
  - tests/test_spec_has_one_home.py
  - decisions/D0100-the-triage-sweep-drops-to-1-2-ns-as-standalone-triage-not-as-the-optimum.md
evidence:
  - '@tt8804: "okay now we can do 10 ns md on the 90 and then 100 ns , so three tiers" and "maybe only need 5 ns for 2nd tier screen. do an optimization consideration for the 3 MD tiers"'
  - 'TIER 1 RESULT: 299 of 300 modes completed at 1.2 ns, 0 errors, 9.4 h on 5 cards; 90 held at the 0.35 nm max-ligand-RMSD bar (30.1%)'
  - 'the 90 are 54 acrylamide + 36 bdhi_c5, max-RMSD 0.150-0.349 nm'
  - 'THE BAR IS NESTED: max ligand RMSD over a longer window is >= the max over a shorter one, so anything holding 100 ns necessarily held 10 ns and 5 ns. A longer tier 2 CANNOT lose a tier 3 survivor'
  - 'so tier-2 length is a PURE COST TRADE. Over the 90, using exp_sweep_length pass rates: 3 ns -> 331 GPU-h total, 5 ns -> 212, 8 ns -> 170, 10 ns -> 122'
  - 'at 10 ns the two stages balance almost exactly (62 GPU-h tier 2 against 61 GPU-h tier 3), which is the optimum'
  - '10 ns is the EDGE of what exp_sweep_length measured (0.1-10 ns at 0.1 ns resolution), so longer is not defensible without new data'
  - 'AGAINST TWO TIERS (1.2 ns straight to 100 ns): 615 GPU-h for the 90. The middle tier saves 493 GPU-h, 80%'
  - 'THE PASS RATES ARE TRANSFERRED, not measured on this cohort: exp_sweep_length UNCONDITIONAL curve over the nac_v5 ranked list, while these 90 are pre-selected by both 1.2 ns and engagement'
  - 'a 100 ns run is 83x a 1.2 ns one at the measured 4.10 min/ns, which is why an hour spent culling buys more than an hour saved'
runbook: bash scratchpad/tier2.sh
---

# D0101 — three tiers, and why the middle one is 10 ns

## The cascade

| tier | length | population |
|---|---:|---|
| 1 — triage | **1,200 ps** | 300 modes → **90 held** |
| 2 — screen | **10,000 ps** | the cull that pays for tier 3 |
| 3 — production | **100,000 ps** | the result |

All three apply the **same** bar: `md.sweep_survivor_rmsd_nm`, 0.35 nm max ligand
RMSD. That is what makes the cascade coherent rather than three unrelated filters.

## Why a middle tier exists

Straight from 1.2 ns to 100 ns costs **615 GPU-h** for the 90 survivors. Culling
at 10 ns first costs **122** — an **80% saving**. A 100 ns run is **83×** a 1.2 ns
one at the measured 4.10 min/ns, so an hour spent culling buys far more than an
hour saved.

## Why 10 ns and not 5

**The bar is nested.** Max ligand RMSD over a longer window is ≥ the max over a
shorter one, so anything that holds 100 ns necessarily held 10 ns and 5 ns.
**A longer tier 2 cannot lose a tier 3 survivor.** It is a pure cost trade with no
recall risk, and culling harder therefore always wins:

| tier 2 | → tier 3 | tier 2 GPU-h | tier 3 GPU-h | total |
|---:|---:|---:|---:|---:|
| 3 ns | 46 | 18 | 313 | 331 |
| 5 ns | 27 | 31 | 182 | 212 |
| 8 ns | 18 | 49 | 121 | 170 |
| **10 ns** | **9** | **62** | **61** | **122** |

At 10 ns the two stages balance almost exactly — 62 GPU-h against 61 — which is
the optimum of a cascade like this. The instinct that a shorter middle tier is
cheaper is right about tier 2 in isolation and wrong end to end.

**10 ns is also the edge of the evidence.** `exp_sweep_length` measured 0.1–10 ns
at 0.1 ns resolution; going longer would be extrapolation, and the balance point
being at the boundary is a reason to stop there rather than a reason to push on.

## What this does NOT establish

* **The pass rates are transferred, not measured on this cohort.** They come from
  `exp_sweep_length`'s **unconditional** curve over the nac_v5 ranked list. Our 90
  are pre-selected by *both* the 1.2 ns triage and engagement, so their true rate
  at 10 ns is probably higher and tier 3 will receive more than 9. The nesting
  argument is untouched by this; the arithmetic is an estimate and the run will
  replace it.
* **Nothing here says a survivor binds.** The cascade measures whether a pose
  stays where it was put. `rank_validated` remains False.
* **Tier 1's own length is still a triage compromise** (D0100): 1.2 ns sits
  outside D0085's 4.3–9.5 ns CI and passes 36% where 5 ns passes 11%, so some of
  the 90 are false positives. Tier 2 is what culls them — which is the second
  job it is doing, and the reason its bar was set deliberately rather than
  inherited.

## Guards

`md_tier_ps()` is an **allowlist of three** and raises on a fourth name, so a
caller cannot invent a tier and receive a default. `test_the_three_md_tiers_are_a_cascade_at_one_bar`
asserts the lengths strictly increase and that tiers 1 and 3 agree with their
existing readers — a middle tier no longer than the triage cannot cull anything,
and one no shorter than production cannot save anything.
