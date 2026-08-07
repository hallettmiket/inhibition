# Class stratification — why composition effects cannot drive selection

*Recorded 2026-08-07 at @tt8804's direction, closing the open question in
[#21](https://github.com/hallettmiket/inhibition/issues/21) / D0073.*

---

## The concern

D0073 measured that pose consensus **depletes** validated chemistry rather than
enriching for it:

| | validated-mechanism share |
|---|---:|
| the library | **90.3%** |
| the consensus-surviving pool | **77.8%** |

Fisher odds ratio **0.34**, p = 2.3×10⁻¹⁴. Mechanisms with no positives and no
validation — BDHI above all — pass consensus at more than twice the base rate
(`sn2_ring_opening` 16.6%, lift 2.41×).

#29 found the same thing in a second metric. Against T₃ acrylamide,
`anchor_quality` separates by class with large effect:

| class | median `anchor_quality` | Cliff's d vs T₃ acrylamide |
|---|---:|---:|
| bdhi_c5 | 0.3723 | **+0.632** |
| bdhi_c4 | 0.3042 | **+0.574** |
| chloroacetamide | 0.1466 | +0.003 |
| acrylamide | 0.1314 | −0.011 |
| naphthoquinone_benzo | 0.0658 | −0.343 |
| sulfamate_acetamide | 0.0255 | **−0.740** |

So on two independent scores, the unvalidated chemistry outranks the validated
chemistry. A **globally** ranked list would hand elevation to BDHI.

## Why it does not

**Ranking and selection are both stratified by warhead class.** This was already
the design; it is recorded here because the D0073 finding is easy to read as an
unaddressed threat when it is in fact already controlled.

`scripts/rank_v2.py` computes rank **within** class, never across:

```python
for cls, g in df.groupby("warhead_class"):
    s["class_rank"] = s[score].rank(ascending=False, method="min")
```

`scripts/select_elevate.py` then takes the top *n per class*:

```python
for (tier, cls), g in surv.groupby(["tier", "warhead_class"]):
    for r in g.nsmallest(args.per_class, "class_rank").itertuples():
```

The consequence is the point: **a molecule competes only against others carrying
the same warhead.** BDHI's higher scores buy it a better position among BDHIs and
buy it nothing against acrylamides. The share of each mechanism in the elevated
set is fixed by `--per-class`, not by the score — so a composition effect, however
large, cannot change what gets elevated.

This also explains why T₄ being **class-balanced by construction** (187 molecules
in each of 9 classes) is not the problem it looks like. A balanced library would
badly distort a global ranking. Under stratified selection it is simply the
design.

## What this does and does not settle

**Settled:** the depletion cannot leak into selection. That was the operational
risk, and it is controlled by construction rather than by monitoring.

**Not settled, and deliberately left open:** *why* unvalidated chemistry scores
higher. Two live explanations, not distinguished by any measurement we have:

1. **BDHI really does form better near-attack geometry**, and the validated
   chemistry is validated for reasons — potency, selectivity, tractability — that
   our geometric scores do not measure.
2. **The scores are biased toward BDHI's shape.** `anchor_quality` multiplies a
   distance and an angle factor; a rigid strained ring reaches a narrow angular
   window more reliably than a flexible acrylamide, independent of whether it
   would react.

Explanation 2 would mean the scores are partly measuring conformational rigidity
under a different name. Distinguishing them needs measured activity for BDHI, of
which there is none — which is exactly why the class is unvalidated, and the
circularity is why stratification rather than re-weighting is the right response
for now.

**A caution on re-measuring.** [#30](https://github.com/hallettmiket/inhibition/issues/30)
has since shown that the top-10-by-energy window these scores are computed over
contains the crystal pose only **33.3%** of the time. Both `consensus_*` and
`topn_viable_frac` read that window, so the class differences above are measured
through a mostly-wrong lens. When #23's mode-diverse keep rule lands, **the class
comparison should be re-run** — the sign is not guaranteed to survive, and
nothing should be re-weighted on the current numbers.

## Related

- `decisions/D0073-consensus-selects-away-from-the-validated-chemistry.md`
- [#21](https://github.com/hallettmiket/inhibition/issues/21) — the review thread
- [#23](https://github.com/hallettmiket/inhibition/issues/23) — the pose window
- [#29](https://github.com/hallettmiket/inhibition/issues/29) — the second measurement
