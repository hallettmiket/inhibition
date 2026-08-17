---
id: D0085
title: The triage sweep is 8 ns and the survivor bar is 0.35 nm max ligand RMSD
date: 2026-08-16
status: accepted
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/target.yaml
  - scripts/exp_sweep_length.py
  - docs/sweep_length.md
  - scripts/sweep_combine.py
evidence:
  - 'the design (@tt8804): poses are ranked on docked geometry -- the best case IF the pose held -- and the short MD asks only whether the pose is stable; survivors earn 100 ns'
  - 'truncation is one-sided: max@x <= max@10, so a shorter sweep can never drop a genuine survivor from the set, only add extras'
  - 'exchange rate: one 100 ns run (4.5 h observed) costs as much as 66 ns of sweeping one mode (4.10 min/ns over 251 completed sweeps)'
  - 'measured at 0.1 ns resolution over 168 finished sweeps; at the 0.35 nm bar the optimum is 8.3 ns, bootstrap 95% CI 4.3-9.5 ns'
  - 'the optimum is a PLATEAU: pass rate is a step function moving 0.60% per crossing, total variation across 8-10 ns is 3.7%, and between 9.0 and 9.6 ns no trace crosses at all'
  - 'the earlier "9.0 ns on the dot" was spurious precision -- the argmax lands wherever one of 168 molecules happens to cross'
  - 'at 0.35 nm, 12 of 168 modes survive (7.1%) against 39 at 0.5 nm'
  - 'CONTRARY: of the first three 100 ns runs the only one that held (0.83 nm) had a sweep RMSD of 0.483, which this bar rejects; the two that flew out at ~6 nm were at 0.465 and 0.473'
  - 'the 10 ns failure hazard was still running at 9 ns (0.015/ns, down from 0.045), so 10 ns is not a converged verdict'
---

# D0085 — the triage sweep is 8 ns, the survivor bar is 0.35 nm

## The design this serves

@tt8804: *"we order by best case scenario and incrementally see how good the
reality is"*, and *"we run x ns MD down a ranked list given time constraints and
then only run 100 ns MD on mols that stayed under 0.5 rmsd for the x ns"*.

Three stages, each cheaper than the next and each asking a different question:

1. **Docking + NAC ranking** — how good would this pose be *if it held*? Best
   case, and the order everything downstream is spent in.
2. **Short MD (this decision)** — does it hold at all? A stability gate, nothing
   more. Engagement is not re-asked here; it was the ranking.
3. **100 ns MD** — the real filter.

Because stage 2 only gates, its cost must be small relative to stage 3, and its
job is to be *cheap and enriching*, not to be a good predictor of stage 3.

## Sweep length: 8 ns

Measured at 0.1 ns resolution over the 168 finished sweeps
(`scripts/exp_sweep_length.py`, written up in `docs/sweep_length.md`).

The estimand is **genuine survivors found per GPU-hour**:

    cost(x) = x · 4.10 min  +  pass_rate(x) · 4.5 h

The yield term is constant in `x` because truncation is one-sided: `max@x ≤
max@10`, so every genuine survivor passes at every `x`. Shortening cannot lose
one from the survivor set — it only lets extras through, and each extra costs a
full 100 ns run.

**The exchange rate decides it.** One 100 ns run costs as much as **66 ns of
sweeping one mode**. The gate is very cheap relative to what it gates, so buying
coverage by weakening it is a bad trade. At the 0.5 nm bar, 5 ns passes 38% of
modes while only 23% are genuine — 18% less efficient than the optimum.

### The optimum is a plateau, and the first version of this overstated it

An earlier draft recorded "9.0 ns" as the answer. That was **spurious
precision**. Pass rate is a *step function* — it moves only when one of 168
traces crosses the bar, 0.60% per step — so the argmax lands wherever a molecule
happens to cross. Between 9.0 and 9.6 ns **no trace crosses at all**, and the
apparent decline there is pure sweep-cost arithmetic. Total variation across
8–10 ns is **3.7%**, well inside noise.

What the data does support:

* below ~7 ns the curve falls away steeply and is outside the bootstrap interval
* past 10 ns it declines on sweep cost alone
* anywhere in **8–10 ns** the choice is free

**8 ns** is chosen as the short end of that plateau: it is inside the interval
(optimum 8.3 at this bar) and saves 2 ns × 4.10 min on every mode swept.

## Survivor bar: 0.35 nm max ligand RMSD

@tt8804. A strict reading of "minimally moving", and deliberately not the
**1.2 nm** used for the 100 ns *held* verdict — that is a different question
asked over a window ten times longer.

It is severe: **12 of 168 modes pass (7.1%)**, against 39 at 0.5 nm. The screen
is explicitly trading recall for cost. Only poses that barely move earn 100 ns.

The bar was checked physically before being tightened. As it tightens, the
warhead stays on the anchor:

| bar (nm) | kept | median warhead–Sγ | within 5 Å |
|---|---|---|---|
| 0.3 | 4 | 3.96 Å | 100% |
| 0.5 | 39 | 4.22 Å | 67% |
| 1.0 | 124 | 5.47 Å | 41% |
| all | 168 | 6.39 Å | 32% |

At 0.5 the median warhead is already inside the 2.8–4.2 Å attack window; at 1.2
it has drifted to 5.6 Å. So the loose bar kept poses that stayed in the *site*
while losing the *presentation*, which is the thing being screened for.

## The evidence that points the other way

Recorded here rather than lost, because it is the only direct test we have.

Of the first three completed 100 ns runs:

| mode | sweep max RMSD | 100 ns max RMSD |
|---|---|---|
| t4_0c53d06ce193_m3 | **0.483** | 0.83 — **held** |
| t4_35344d777bbe_m1 | 0.465 | 6.01 — gone |
| t4_966a63abe1d6_m2 | 0.473 | 5.78 — gone |

**The only one that held had the loosest sweep RMSD of the three, and a 0.35 nm
bar rejects all three.** n = 3 decides nothing, but if that pattern survives more
data the bar is cutting into real candidates and this decision should be
revisited.

Related: the 10 ns failure hazard was still running at 9 ns (0.015/ns, down from
0.045 early), so "held for 8 ns" is a weak claim about 100 ns behaviour. Under
this design that is acceptable — stage 3 is the filter — but it is the reason the
sweep should not be trusted as a predictor.

## What this cannot settle

`exp_sweep_length` defines "genuine" against the **10 ns** verdict, because too
few 100 ns runs exist to define it against 100 ns. So it finds the cheapest `x`
that reproduces the 10 ns answer, not the best `x` for predicting 100 ns. The
script takes `--truth-md100` for re-running once enough rows land, and it should
be re-run then.

## The larger lever

At the optimum, the 100 ns stage is **64% of total cost per mode**. Sweep length
spans 113 → 138 candidates per 1000 GPU-hours across its entire useful range.
Halving the 100 ns stage — shorter production, or a cheaper intermediate before
committing to it — would move throughput far more than any choice of `x`.

## Recorded

```yaml
md:
  sweep_ps: 8000
  sweep_survivor_rmsd_nm: 0.35
  bound_rmsd_nm: 1.2        # the 100 ns held verdict, unchanged
```

Current queue: 12 survivors, 5 already run at 100 ns, **7 queued** on 3 GPUs.
