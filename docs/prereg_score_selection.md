# Pre-registration — choosing the 2.2.0 ranking score

*Written before any score is evaluated. @tt8804, 2026-08-07: "no rush we want to
do our best with ranking."*

---

## Why this document exists

I am about to compare six candidate scores against a validation set. **Picking
the winner after seeing the results is how a score gets chosen by noise** — and
this project has already shipped one ranking (`topn_viable_frac`) that placed the
crystallographically-confirmed parent compound **dead last of 5,765**. So the
criteria, the thresholds, and what disqualifies a score are fixed here first.

## What we are choosing between

| score | definition |
|---|---|
| `viable_fraction` | fraction of all poses in a mode that are reaction-competent |
| `enrichment_joint` | viable_fraction ÷ isotropic null — **2.0.0's score** |
| `enrichment_conditional` | P(viable \| in distance window) ÷ null |
| `conditional_x_consensus` | the above × mode population — *currently wired* |
| **`conditional_lcb`** | **Wilson 95% lower bound** on P(viable \| in window) ÷ null |
| `anchor_quality_max` | best continuous anchoring geometry in the mode |

## The three tests

### T1 — convergence *(disqualifying)*

D0068 has demanded this since the beginning and **no score in this project has
ever passed one.** Measured against the two things that should not change a
molecule's score:

- **sampling depth**: ρ between the score at 500 runs and at 2,000 runs, on the
  same molecules
- **re-run**: ρ between two independent 500-run dockings of the same molecules

Already known and alarming: between 200 and 500 runs the *current* score gives
**ρ = +0.395**, a median rank move of 853 places, and only 3 of the old top-10
surviving into the new top-10.

### T2 — reference separation *(primary)*

The 22 reference molecules carry potency annotations, and **this has never been
used to validate a ranking.** Ten carry a warhead this criterion can score.

- **leads** — Sulfopin (nM covalent), BJP-06-005-3, Liu-2022-ZL-Pin13 (IC50
  67 nM), Reddi-2023-4d/4g, Tian-6a, Ieda-(S)-2
- **historical promiscuous** — KPT-6566, Juglone, ATRA

**Readout:** AUC separating leads from promiscuous, and Sulfopin's percentile
against the library.

### T3 — independence from docking energy *(disqualifying)*

A score rank-correlated with AutoDock energy has re-imported the defect this
version exists to remove (#23/#30: the correct pose sits at a uniformly random
energy rank, KS *p* = 0.666). **Readout:** |ρ(score, mean_energy)|.

## Readings, fixed now

| criterion | disqualifies | acceptable | good |
|---|---|---|---|
| **T1 sampling** ρ(500, 2000) | < 0.5 | 0.5–0.7 | ≥ 0.7 |
| **T1 re-run** ρ(500a, 500b) | < 0.6 | 0.6–0.8 | ≥ 0.8 |
| **T2 Sulfopin percentile** | < 25th | 25–60th | ≥ 60th |
| **T2 lead-vs-promiscuous AUC** | < 0.5 | 0.5–0.7 | ≥ 0.7 |
| **T3** \|ρ(score, energy)\| | > 0.4 | 0.2–0.4 | ≤ 0.2 |

**Selection rule, fixed before the numbers:** any score failing a *disqualifying*
criterion is out regardless of the others. Among survivors, rank on **T1 sampling
convergence first**, then T2 AUC. Convergence comes first because a score that
does not reproduce cannot be validated by anything — a good AUC on an
irreproducible score is a good measurement of that day's noise.

**If every score is disqualified, that is the result.** The honest conclusion
would be that pose geometry at this sampling depth does not support a ranking,
and the next move is more sampling or a different observable — not a seventh
score.

## What this cannot settle

- **n = 10 references, 3 of them promiscuous.** An AUC on 7 vs 3 is a weak
  instrument; it can catch a catastrophic score, not distinguish two good ones.
- **The promiscuous three are all naphthoquinone-type Michael acceptors** and
  most leads are not, so the comparison is confounded with warhead class. Any AUC
  above chance may be measuring chemistry rather than quality. Reported alongside
  a within-class version wherever the class has both.
- **Non-covalent leads cannot be scored at all** — 12 of 22 references have no
  warhead this criterion applies to, including the most potent (Wildemann, Ki
  1.2 nM). The validation set is therefore biased toward covalent chemistry by
  construction.
- **Convergence is necessary, not sufficient.** A score can be perfectly
  reproducible and measure the wrong thing; `enrichment_joint` was reproducible
  and dimensionally inconsistent.
