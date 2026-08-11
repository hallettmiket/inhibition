# Sweep depth: how a target-agnostic screen decides when to stop

*Method note for the paper. @tt8804, 2026-08-11, [#59](https://github.com/hallettmiket/inhibition/issues/59).
Written before tonight's data lands, so the readings are fixed in advance.*

---

## 1. The question, stated as an estimand

A screen ranks **modes** — a molecule's docked poses clustered into binding modes,
each mode a candidate row. Ranking is cheap. The next stage, a 10 ns attack sweep,
costs ~0.5 GPU-hours per mode. With 8,096 ranked modes, sweeping everything costs
~4,000 GPU-hours and is not affordable.

So the screen needs a stopping rule. The configuration for a target currently reads

```yaml
target:     Pin1
domain:     catalytic
pdb:        3IKD
anchor:     Cys113
sweep_depth: ?
```

and the question is what belongs in the last line.

**It is not a number.** Depth as "sweep ranks 1..D" is a quantity with no
transferable meaning: rank is an ordinal within one warhead class of one library,
so D depends on library size, class composition, and the difficulty of that
class's geometric criterion. A D fitted on Pin1/3IKD says nothing about the next
target.

The estimand that *does* transfer is

> **P(the sweep returns a productive mode | a quantity computable at ranking time)**

and the design question is which quantity, and where to cut it.

## 2. Rank does not supply it

Over the 188 swept modes carrying both a rank and a score, on this target:

| predictor (free at ranking time) | ρ vs `frac_attack_ready` | AUC for "productive" |
|---|---:|---:|
| `class_rank` | +0.053 (p = 0.47) | **0.518** |
| `viable_fraction` | +0.198 (p = 0.006) | 0.567 |
| `enrichment` | +0.188 (p = 0.009) | 0.563 |
| `conditional_eb` | +0.133 (p = 0.07) | 0.536 |
| `consensus` | +0.035 (p = 0.63) | 0.552 |
| `anchor_quality_mean` | +0.184 (p = 0.011) | 0.546 |

*productive* = `frac_attack_ready > 0`; base rate 109/192 = 57%. AUC is the
Mann–Whitney statistic, so 0.5 is chance.

Rank is at chance. Note also the **sign**: within the two classes with real sample
size, ρ(class_rank, frac_attack_ready) is **positive** — bdhi_c4 +0.261 (n = 48),
bdhi_c5 +0.239 (n = 40) — i.e. worse-ranked modes swept marginally *better*.
Neither is significant, but nothing here supports sweeping top-down.

## 3. Why that table is not yet evidence that these predictors are useless

**The sample is range-restricted at exactly the end that matters.**

| | n | median `viable_fraction` | p25 | fraction at zero |
|---|---:|---:|---:|---:|
| swept modes | 192 | 0.173 | 0.135 | **2.6%** |
| whole library | 6,579 | 0.116 | 0.078 | **5.6%*** |

\* of all 8,096 ranked modes, **1,869 have `viable_fraction == 0`. Five have ever
been swept** — 0.27% of the stratum.

An AUC computed inside the top of a predictor's range is not an estimate of its
discrimination; restriction of range attenuates every association toward chance.
The honest reading of §2 is **"unmeasured over most of the range"**, not
"useless". This is the single most important caveat in this note, and it is the
reason the experiment in §6 is designed the way it is.

## 4. The proposed parameter

Use **`viable_fraction`** — the fraction of a mode's docked poses that satisfy the
near-attack criterion (distance window plus the mechanism's angular bar).

It has the four properties a target-agnostic threshold needs:

1. **Free.** Computed at ranking time, before any dynamics.
2. **Dimensionless and bounded**, so a threshold is a number between 0 and 1
   rather than a position in a list.
3. **Defined identically for every class and every target.** It does not depend on
   library size, and unlike rank it is not an ordinal.
4. **Mechanistically the right question.** It asks whether the docked ensemble
   contains *any* attack-competent geometry. A mode with `viable_fraction = 0` has
   no pose from which the reaction could proceed; the sweep would have to
   *manufacture* one through dynamics.

Its mechanism-specific null is already computed — `isotropic_null`, 0.0816 for
`sn2_ring_opening`, `michael_addition` and `snar_displacement`, 0.0670 for
`sn2_displacement` — so the comparable form is

> **`enrichment` = `viable_fraction` / `isotropic_null`**

which normalises away the fact that the SN₂ angular criterion (150°) is a far
narrower target than the perpendicular one (30° off-normal), the bias behind
[#47](https://github.com/hallettmiket/inhibition/issues/47). **`enrichment` is the
quantity a threshold should be carried on**; it is already in `rank_v2`.

**Depth becomes derived, not configured.** Sweep every mode with
`enrichment ≥ E*`; the depth D is whatever rank that reaches in a given class, and
differs per class by construction.

## 5. Is depth target-agnostic?

Three separable layers, and they answer differently:

| layer | agnostic? | why |
|---|---|---|
| the **rule** — "sweep while enrichment ≥ E*" | **Yes** | no target-specific term |
| the **threshold** E* | **Probably**, pending evidence | enrichment is already null-normalised per mechanism; if P(productive \| enrichment) has the same shape across classes, one E* serves |
| the **depth** D | **No, and it should not be** | D is the rank at which a given class's enrichment distribution crosses E*. Different libraries cross it in different places, which is the correct behaviour |

So the config becomes

```yaml
target:     Pin1
domain:     catalytic
pdb:        3IKD
anchor:     Cys113
sweep_rule:
  parameter:        enrichment      # viable_fraction / isotropic_null
  floor:            E*              # measured, see §6
  pilot:            stratified      # required before floor is trusted on a new target
  capture_target:   0.95            # fraction of productive modes the floor must retain
```

`sweep_depth` never appears. What appears is a floor and the pilot that justifies
it.

**The open question this note cannot close:** whether E* is one number per target
or one per (target, warhead class). The null is per *mechanism*, which argues that
enrichment is already comparable across classes and one floor should serve. If
tonight's per-stratum curves differ by class beyond their confidence intervals,
that argument fails and the floor becomes per class.

## 6. The experiment, with its readings fixed in advance

A stratified sample over `viable_fraction`, drawn across all nine classes so no
stratum is one class's result in disguise. 72 modes:

| stratum | n | what it is for |
|---|---:|---|
| **v = 0** | **30** | the decisive one. 1,869 modes, 5 ever swept |
| 0 < v ≤ 0.05 | 16 | the transition |
| 0.05 < v ≤ 0.10 | 10 | |
| 0.10 < v ≤ 0.15 | 8 | |
| v > 0.15 | 8 | anchor against the 192 already swept |

Built by `scripts/sweep_gap_worklist.py --v-strata`, interleaved across strata so
a run stopped by the clock still spans the range.

### Power, computed before the run

For a proportion near 0.5, the 95% CI half-width is `1.96·√(p(1−p)/n)`:

| n per stratum | half-width |
|---:|---:|
| 12 | ±28 pp |
| 20 | ±22 pp |
| 30 | ±18 pp |
| 44 | ±15 pp |
| 100 | ±10 pp |

At ~12 per stratum this is **a coarse curve: enough to see a cliff, not enough to
place it precisely**. Resolving a 15 pp difference between adjacent strata needs
n ≈ 44 each, ≈ 264 sweeps ≈ 132 GPU-hours. That is the honest cost of a precise
answer and it is not what tonight buys.

**The zero stratum is the exception, and is why it gets n = 30.** It is a
one-sided question. If 0 of 30 are productive, the rule of three puts the 95%
upper bound on the rate at 3/30 = **10%**, which is decision-grade without a
precise estimate.

### Readings, fixed now

| observation | conclusion |
|---|---|
| **0 of 30 zero-viable modes productive** | Adopt `viable_fraction > 0` as a hard prefilter. Excludes **1,869 modes, 28% of the library**, for free, on every future screen. |
| **1–3 of 30 productive** | Prefilter still worth it as a *deprioritiser*, not a hard gate; quantify the loss and state it |
| **≥ 4 of 30 productive** | The docked ensemble does **not** bound what dynamics can reach. The near-attack criterion is not a filter, and E* cannot be set from docking at all — depth must then be set by cost, not by prediction |
| **monotone rise in productivity across the four non-zero strata** | Fit P(productive \| enrichment), set E* at the 95% capture point |
| **flat across the non-zero strata** | Enrichment discriminates only zero-vs-nonzero. Use it as a prefilter and choose depth by budget |

### What is deliberately not claimed

- *Productive* here means `frac_attack_ready > 0` in a 10 ns sweep. It is a
  geometric readout, not evidence of reaction, and D0075/D0076 record that this
  sweep's own gating has already been wrong once.
- 233 of the 239 previously swept modes are **mode 0**. Non-dominant modes have
  essentially never been swept, so any curve fitted on historical data is a curve
  about dominant modes. Tonight's ladder is the first data that is not.
- The strata are drawn from unswept modes, which are themselves a non-random
  sample of the library — they are what the 2.2.0 selection skipped
  ([#53](https://github.com/hallettmiket/inhibition/issues/53)). Spreading each
  stratum across nine classes mitigates this; it does not remove it.

## 7. What this contributes to the paper

The generalisable claim is not a depth. It is:

> A staged covalent screen should not be configured with a depth. It should be
> configured with a **null-normalised geometric floor** and a **stratified pilot
> that measures where that floor sits for the target at hand**, because the rank
> at which a library stops being worth simulating is a property of the library,
> while the floor is a property of the chemistry.

The falsifiable part is §6's first reading: that a mode whose docked ensemble
contains no near-attack pose does not become attack-ready under 10 ns of
unbiased dynamics. If that holds, it is a cheap and portable prefilter. If it
fails, the docking-derived criterion cannot bound the search and this whole
framing has to be replaced with a cost-based cut.
