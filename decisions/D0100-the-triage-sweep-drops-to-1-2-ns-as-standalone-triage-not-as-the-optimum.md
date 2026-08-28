---
id: D0100
title: The triage sweep drops to 1.2 ns as STANDALONE triage — outside the measured plateau, safe only because truncation is one-sided
date: 2026-08-27
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/target.yaml
  - decisions/D0087-the-triage-sweep-drops-to-5-ns.md
  - scripts/engagement_curve_plot.py
evidence:
  - '@tt8804: "first off lets trim down to 1.2ns sweep", then "run the 1.2ns sweep on bdhi and acrylamide using 3 cards", then "BDHIC5 AND ACRYLAMIDE" and "run only on the top pose for each mol"'
  - 'D0085 measured the yield PEAK at 8.3 ns (45.7 per 1000 GPU-h), bootstrap 95% CI 4.3-9.5 ns. 1.2 ns is WELL OUTSIDE that interval, so D0087 "safe because it is inside the CI" argument does not transfer'
  - 'MEASURED TRADE-OFF at 1.2 ns vs 5 ns: pass rate 36.3% vs 10.7%; sweep 0.082 vs 0.342 h/mode; the 100 ns stage it feeds 1.634 vs 0.482 h/mode; yield per 1000 GPU-h 20.8 vs 43.4'
  - 'so the sweep gets 4.2x cheaper and the stage it feeds gets 3.4x dearer, and end to end the yield per GPU-hour HALVES'
  - 'WHAT MAKES IT SAFE: the bar is one-sided, so anything clearing 5 ns also clears 1.2 ns. A shorter sweep cannot lose a survivor; it admits false positives instead'
  - 'ACRYLAMIDE DOES NOT NEED A LAXER ENGAGEMENT CRITERION, checked on 22,486 valid poses: decomposing anchor_quality into its distance and angle terms, acrylamide has the BEST distance term of the three families (median 0.658, only 20% outside the window against 43% for bdhi_c4) and is within a few degrees on angle (22% inside 30 deg against 26%)'
  - 'the term that binds is the 30 deg off-normal window, FOR ALL THREE FAMILIES: median a_term 0.11-0.22 against d_term 0.38-0.66, and median off-normal angle ~50 deg'
  - 'acrylamide and BDHI already share a mechanism class (perpendicular_to_plane) and an isotropic_null (0.0816), so their enrichment is already comparable'
  - 'THE ENGAGEMENT CURVE HAS NO CLIFF: rank 1 scores 0.96, rank 150 about 0.78, rank 300 about 0.76 -- 0.2 of a 0-1 scale across the whole selectable range'
  - 'THE CAP SELECTS, NOT THE FLOOR: 43,559 in-scope groups clear the 0.05 budget floor and the 150/family cap admits 450; the lowest SELECTED engagement is 0.760, 15x the floor'
runbook: bash scratchpad/sweep_top1.sh
---

# D0100 — 1.2 ns, and why that number is not a measurement

## The decision

`md.sweep_ps` goes from 5,000 to **1,200**, for a standalone triage pass over
300 molecules — 150 acrylamide and 150 bdhi_c5, one mode each, on 3 cards.

**It is not the measured optimum and must never be read as one.** D0085 put the
yield peak at **8.3 ns**, bootstrap 95% CI **4.3–9.5 ns**. D0087 could drop to 5 ns
and call it safe *because 5 sits inside that interval*. 1.2 ns does not, so that
argument does not carry over and a different one is needed.

## What actually makes it safe

**The bar is one-sided.** A mode passes the sweep by staying under a max-RMSD
bar; a shorter window can only be a subset of a longer one, so **anything that
would clear 5 ns also clears 1.2 ns.** A short sweep cannot lose a survivor. It
admits false positives instead, and a false positive costs GPU time while a miss
costs a candidate — the asymmetry D0085 chose deliberately.

## What it costs, measured

| | 1.2 ns | 5 ns |
|---|---:|---:|
| pass rate | **36.3%** | 10.7% |
| sweep, h/mode | 0.082 | 0.342 |
| **100 ns stage it feeds, h/mode** | **1.634** | 0.482 |
| total h/mode | 1.716 | 0.824 |
| **yield per 1000 GPU-h** | **20.8** | **43.4** |

The sweep gets **4.2× cheaper** and the stage it feeds gets **3.4× dearer**. End
to end the yield per GPU-hour **halves**. Over the 413-mode worklist that is
34 GPU-h of sweep against 141 — but 675 GPU-h of 100 ns against 199.

**So this is a triage-only setting.** Anyone going on to 100 ns should put it
back to 5 ns first, and the config says so where the value lives.

## Acrylamide does not need a laxer engagement criterion

Asked directly, and checked rather than reasoned about. `anchor_quality` is
`d_term × a_term`; over 22,486 PoseBusters-valid poses:

| family | d_term | a_term | outside distance window | median off-normal |
|---|---:|---:|---:|---:|
| acrylamide | **0.658** | 0.111 | **20%** | 53.3° |
| bdhi_c4 | 0.383 | 0.220 | 43% | 46.8° |
| bdhi_c5 | 0.555 | 0.168 | 38% | 49.9° |

**Acrylamide is the best-placed family on distance**, not the penalised one —
median 3.42 Å against an ideal of 3.5, with less than half bdhi_c4's rate of
falling outside the window. On angle it is within a few degrees of the others
(22% inside 30° against 26%). The two also already share a mechanism class
(`perpendicular_to_plane`) and therefore an `isotropic_null` of 0.0816, so their
enrichment is comparable by construction.

Loosening acrylamide alone would hand it an advantage on the one axis where it is
already competitive, against families held to the original rule.

**The real finding is that the 30° off-normal window binds for everything.**
Median `a_term` is 0.11–0.22 against `d_term` 0.38–0.66, and the median pose sits
~50° off normal — so the angle term is near zero almost everywhere and the score
is dominated by a gate almost nothing clears. That is a chemistry judgement about
what a perpendicular attack tolerates, it applies to all three families equally,
and it belongs in #12 with the Lu lab rather than as a per-family adjustment.

## The engagement curve has no cliff, and the cap is what selects

`scripts/engagement_curve_plot.py`, on nac_v6's 327,167 groups: rank 1 scores
0.96, rank 150 about **0.78**, rank 300 about **0.76**. Across the entire
selectable range the score moves **0.2 on a 0–1 scale**, and all nine warhead
classes trace nearly the same curve.

**43,559 in-scope groups clear the 0.05 budget floor; the 150/family cap admits
450.** The lowest selected engagement is 0.760 — **15× the floor**. The floor is
doing no work at all, and because the curve is flat through that region the cap
is close to arbitrary: mode #150 is barely distinguishable from #300 or #1,000.

That is the same compression the docking energy showed (76% of poses within
2 kcal/mol of the best). Engagement **predicts the MD outcome** where the
incumbent did not (D0098, ρ = +0.652 against −0.015) — but it does not separate
the top of the library from itself. It says a mode is worth simulating; it does
not say which of 43,559 qualifying modes to pick.

**The open question that follows:** is the flat band flat in OUTCOME too? A
stratified sweep spanning engagement 0.95 → 0.60 would answer it for roughly a
third of the cost of sweeping the top band, and it decides whether the cap
matters at all.
