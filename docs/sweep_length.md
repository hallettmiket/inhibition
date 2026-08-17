# How long should the triage sweep be?

**Experiment:** `scripts/exp_sweep_length.py`
**Data:** 168 finished 10 ns sweeps (campaign `sweep_gaps_6`, 3.0.0 Galena)
**Result:** a **plateau at 8–10 ns**, not a point. Adopted: **8 ns** with a
**0.35 nm** survivor bar (D0085). At the 0.35 bar the argmax is 8.3 ns with a
bootstrap 95% CI of 4.3–9.5; at the 0.5 bar it is 9.0 with a CI of 7.1–9.7.
The decimal is not meaningful — see "Is the peak real?" below.

## The design being optimised

Poses are ranked on docked geometry — the best case, *if the pose held*. A short
MD then asks one question: is the pose **stable**? Every mode that stays under a
max-RMSD bar for `x` ns earns a 100 ns run. Budget is spent going down the ranked
list, so the question is not "is `x` ns accurate" but **"what `x` finds the most
real candidates per GPU-hour"**.

## Why `x` cuts both ways

A shorter sweep gets you further down the ranked list for the same money — more
coverage, more candidates. But max RMSD only grows with time, so a shorter sweep
also passes poses a longer one would reject, and **every one of those costs a
full 100 ns run**.

Truncation is one-sided: `max@x ≤ max@10`, so every genuine survivor passes at
every `x`. Shortening cannot drop one from the survivor **set** — it only changes
how many extras come with them. So the yield term is constant in `x` and the
whole question is the cost term.

    cost(x) = x · 4.10 min  +  pass_rate(x) · 4.5 h
    metric  = genuine_rate / cost(x)

## The exchange rate, which is what decides it

One 100 ns run (4.5 h) costs as much as **66 ns of sweeping one mode** (4.10
min/ns, from 251 completed sweeps at a 41 min median). The sweep is very cheap
relative to what it gates, so buying coverage by weakening the gate is a bad
trade almost everywhere.

| x (ns) | pass rate | sweep/mode | 100 ns/mode | total/mode | genuine per 1000 h |
|---|---|---|---|---|---|
| 1.0 | 71.4% | 0.068 h | 3.214 h | 3.283 h | 70.7 |
| 3.0 | 45.2% | 0.205 h | 2.036 h | 2.241 h | 103.6 |
| 5.0 | 38.1% | 0.342 h | 1.714 h | 2.056 h | 112.9 |
| 7.0 | 29.2% | 0.478 h | 1.312 h | 1.791 h | 129.6 |
| 8.5 | 25.6% | 0.581 h | 1.152 h | 1.733 h | 134.0 |
| **9.0** | **23.8%** | 0.615 h | 1.071 h | 1.686 h | **137.7** |
| 10.0 | 23.2% | 0.683 h | 1.045 h | 1.728 h | 134.3 |

At 5 ns, **38% of modes pass but only 23% are genuine** — you would run 100 ns on
15% of everything you sweep for nothing, and each of those costs more than
sweeping another mode to completion. 5 ns is 18% less efficient than 9 ns.

## Is the peak real? No — it is a step artefact

Bootstrap over 2,000 resamples: **95% CI 7.1 – 9.7 ns** at the 0.5 bar
(4.3 – 9.5 at 0.35), median resampled optimum 8.8 ns.

More decisive than the CI is the fine structure. **Pass rate is a step
function** — it moves only when one of 168 traces crosses the bar, 0.60% per
step — so the argmax lands wherever a molecule happens to cross:

| x | pass | yield | % of peak |
|---|---|---|---|
| 8.8 | 24.4% | 136.6 | 99.2% |
| **9.0** | 23.8% | **137.7** | **100%** |
| 9.2 | 23.8% | 136.5 | 99.2% |
| 9.6 | 23.8% | 134.4 | 97.6% |
| 9.7 | 23.2% | 136.0 | 98.8% |

Between 9.0 and 9.6 **no trace crosses** — the decline there is pure sweep-cost
arithmetic, not a property of the molecules. Total variation across 8–10 ns is
**3.7%**.

So reporting "9.0 ns" as the optimum was spurious precision. The defensible
statement is **8–10 ns, choice free inside it**. What the data *does* exclude is
the short end: below 7 ns falls away steeply and is outside the interval.

## Sensitivity

| 100 ns cost | optimum | | RMSD bar | optimum |
|---|---|---|---|---|
| 2.0 h | 7.3 ns | | 0.4 nm | 8.6 ns |
| 3.0 h | 7.3 ns | | 0.5 nm | 9.0 ns |
| 4.5 h | 9.0 ns | | 0.6 nm | 8.8 ns |
| 6.0 h | 9.0 ns | | 0.8 nm | 9.9 ns |
| 8.0 h | 9.0 ns | | 1.0 nm | 7.3 ns |

Stable. Only a much cheaper 100 ns stage (≤3 h) pulls the optimum below 8 ns, and
that is the lever worth pulling — see below.

## What this cannot settle

"Genuine" is defined against the **10 ns** verdict, not against 100 ns, because
only a handful of 100 ns runs exist. The 10 ns failure hazard was still running
at 9 ns (0.015/ns, down from 0.045 early), so some 10 ns survivors will fail at
100 ns. Under this design that is *not a problem* — the 100 ns stage is itself
the filter — but it does mean this experiment finds the cheapest `x` that
reproduces the 10 ns answer, not the best `x` for predicting 100 ns.

Early 100 ns results are consistent with that caveat: of the first three
completed, one held (0.83 nm) and two went to ~6 nm, and the one that held had
the *loosest* sweep RMSD of the three. Re-run this experiment against 100 ns
truth once enough rows exist.

## The larger point

Sweep length is a **second-order lever**: the whole range from 5 to 10 ns spans
113 → 138 candidates per 1000 GPU-hours. At the optimum the 100 ns stage is
**64% of total cost per mode** (1.071 h of 1.686 h). Halving *that* — shorter
production, or a cheaper intermediate stage before committing to 100 ns — moves
throughput far more than any choice of `x`.

## Decision

**8 ns, 0.35 nm** — D0085. Recorded in `config/target.yaml` as
`md.sweep_ps: 8000` and `md.sweep_survivor_rmsd_nm: 0.35`.

8 ns is the short end of the plateau. 0.35 nm is a strict reading of "minimally
moving" and passes 12 of 168 modes (7.1%) against 39 at 0.5 nm — an explicit
trade of recall for cost.

**The evidence against the tight bar, recorded rather than lost:** of the first
three 100 ns runs, the only one that held (0.83 nm) had a sweep RMSD of 0.483,
which this bar rejects; the two that flew out to ~6 nm sat at 0.465 and 0.473.
n = 3 decides nothing, but if it holds up the bar is cutting into real
candidates.
