---
id: D0061
title: No cheap pose-selection rule beats random, and the bar is random rather than the docking score
date: 2026-08-05
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/pose_selection_bench.py
  - docs/ranking_rationale.md
evidence:
  - '3IKD, 82 crystal cases: ceiling (a <=2A pose exists in the nine) 41.5%, random pick 19.8%, Vina score top-1 18.3%'
  - 'centroid_closest 24.4% (z=+1.57, p=0.12), contact_medoid 22.0% (z=+0.74, p=0.46), largest_cluster_medoid 19.5%, vina_score_top1 18.3% — none significant'
  - '6VAJ, same 82 cases: ceiling 15.9%, random 5.3%; every rule within noise'
  - 'three of four rules FLIP SIGN between the two receptors'
---

# The docking score is at chance, so chance is the bar

## The reframe

The obvious target for a pose-selection rule is "beat the docking score". On
3IKD that target is **worthless**: the score picks a correct pose 18.3% of the
time and picking **uniformly at random** among the nine gets **19.8%**. The score
is not a weak selector, it is *indistinguishable from a coin flip* — and
marginally below it.

So the bar is **random**, and it is slightly *harder* than the status quo. Any
rule that beats the score but not random has achieved nothing.

## What was measured

Four rules, each mapping a pose ensemble to one index, none permitted to see the
crystal reference:

| rule | 3IKD | vs random | z | p |
|---|---:|---:|---:|---:|
| `centroid_closest` | 24.4% | +4.6% | +1.57 | 0.12 |
| `contact_medoid` | 22.0% | +2.2% | +0.74 | 0.46 |
| `largest_cluster_medoid` | 19.5% | −0.3% | −0.09 | 0.93 |
| `vina_score_top1` | 18.3% | −1.5% | −0.51 | 0.61 |

**None is distinguishable from random.** The best, `centroid_closest`, is 3.8
cases of difference over 82.

**And three of the four flip sign between receptors** — `centroid_closest` is
+4.6% on 3IKD and −1.6% on 6VAJ. A rule capturing something real about pose
quality should help on both. This is noise with a direction.

## Why it was nearly reported as a result

`centroid_closest` beating random by **+4.6%** reads as a finding. It survives
until you attach the number of cases behind it. The null here is a
**Poisson-binomial** — each case has its own success probability, since cases
differ in how many of their nine poses are correct — and pooling to a single rate
would understate the variance and manufacture significance.

`significance()` is now part of the harness and the report will not describe a
rule as beating random without z and p attached.

## What this means for the pipeline

`docs/ranking_rationale.md` stage 2 said: if nothing free clears the floor, that
is the strongest available argument that **poses inside one Vina ensemble are not
separable by cheap geometry**, and BPMD becomes the only remaining option rather
than one of several. That is where we are.

It also sharpens what BPMD must do. The headroom is 19.8% → 41.5%, and stability
under bias is information the docking search does not already contain, which
every rule tested here was — in one form or another — a re-reading of.

**Not yet tested and genuinely different: replicate consensus.** All four rules
above read a *single* ensemble, whose nine modes Vina deliberately spreads with a
minimum-RMSD floor. Consensus across *independent* runs with different seeds is
new information, is what #10 identifies as the honest version, and costs minutes.
It is the last cheap candidate; if it also fails, the free options are exhausted.
