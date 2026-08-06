# Results — which ranking metric selects for physical stability?

*Run 2026-08-06 against the pre-registration in
[`elevation_prereg.md`](elevation_prereg.md), which was committed before any
simulation ran. Cohort: `elevation_cohort_1.csv` (29 molecules, 4 groups,
rotatable-bond matched at 4.0 across the three BDHI groups) plus 8
crystallographic Cys113 positives added as the anchor. Receptor: chemist-prepared
**3IKD** (D0059). @tt8804.*

---

## The answer, in one line

**Neither metric predicts stability — and the anchor proves that is a real
result rather than a dead assay.** The three BDHI groups are statistically
indistinguishable from one another on both tiers, while every one of them is
significantly less stable than molecules known to react with Cys113. The
measurement discriminates; enrichment and consensus do not predict what it
discriminates.

That is the prereg's third reading, **A ≈ B ≈ D**, whose fixed conclusion is:

> Neither metric predicts stability. BPMD is measuring something orthogonal, and
> the ranking has no physical support from this experiment.

One amendment is forced by the anchor, and it makes the result *stronger*, not
weaker: **"orthogonal" is too generous.** Tier 1 separates crystallographic
positives from generated candidates at p = 0.007 with Cliff's δ = −0.78. It is
measuring pose survival, which is exactly what it was built to measure. What has
no physical support is the *ranking*, not the *assay*.

---

## What was run

| | |
|---|---|
| molecules | 37 — A 8, B 8, D 8, V 5, REF 8 |
| replicas | 3 per molecule, per tier |
| tier 1 | energy minimisation + 100 ps NVT + 200 ps NPT, **no position restraints** |
| tier 2 | well-tempered BPMD along d(warhead → Cys113 SG), 3 × 3 ns |
| receptor | `3IKD_noligand.pdb`, Cys113 = residue 63 of 115, SG asserted at (13.385, 3.989, −2.040) |
| failures | tier 1: **0 of 111** |

**The anchor's membership rule was fixed before the run**, because
`crystal_positives` returns 15 ligands and the prereg budgets ≤ 8: sorted by
ident, take the first 8. It depends on nothing but the PDB accession codes.
It yields 5 chloroacetamides, 2 naphthoquinones and 1 chloroazine.

---

## Tier 1 — did the docked pose survive plain dynamics?

`gromacs_explicit` applies no position restraints during NVT/NPT, so 300 ps of
unrestrained dynamics runs before any bias. The pre-registered readout is
**|Δd|**, the warhead-to-SG displacement between the docked pose and the frame
production starts from, averaged over 3 replicas. **Smaller is more stable.**

| group | n | median \|Δd\| (nm) | IQR | median replica spread |
|---|---:|---:|---|---:|
| **A** hi-enr / hi-cons | 8 | 0.277 | 0.251 – 0.283 | 0.135 |
| **B** lo-enr / hi-cons | 8 | 0.198 | 0.183 – 0.269 | 0.148 |
| **D** lo-enr / lo-cons | 8 | 0.204 | 0.123 – 0.494 | 0.050 |
| **V** chloroacetamide | 5 | 0.203 | 0.077 – 0.246 | 0.147 |
| **REF** crystallographic | 8 | **0.102** | 0.074 – 0.126 | 0.037 |

### The pre-registered contrasts: nothing

| contrast | median 1 | median 2 | Cliff's δ | p | p (Holm) | verdict |
|---|---:|---:|---:|---:|---:|:--:|
| A vs B | 0.277 | 0.198 | −0.469 | 0.130 | 0.391 | ≈ |
| B vs D | 0.198 | 0.204 | −0.062 | 0.878 | 0.884 | ≈ |
| A vs D | 0.277 | 0.204 | −0.250 | 0.442 | 0.884 | ≈ |

δ is signed so that **positive = group 1 more stable**. All three point
estimates are negative, i.e. the *higher*-enrichment group is if anything the
*less* stable one — but none is significant and **no claim is made from the
direction of a non-significant effect.** The prereg's fourth reading (D ≥ A, B,
"report as a failure") requires a significant contrast in that direction, and
there is none.

### The anchor: everything

| contrast | median 1 | median 2 | Cliff's δ | p | verdict |
|---|---:|---:|---:|---:|:--:|
| A vs REF | 0.277 | 0.102 | −0.781 | **0.0070** | REF more stable |
| B vs REF | 0.198 | 0.102 | −0.750 | **0.0104** | REF more stable |
| D vs REF | 0.204 | 0.102 | −0.594 | **0.0499** | REF more stable |
| V vs REF | 0.203 | 0.102 | −0.300 | — | *descriptive only, n = 5* |

**This is what makes the null interpretable.** Without REF, "A ≈ B ≈ D" is
consistent with an assay that cannot separate anything. With REF, the same assay
separates known reactive molecules from generated candidates in all three
groups — so the null is a statement about the metrics, not about the
measurement.

---

## Post-hoc observations (NOT pre-registered — labelled as such)

These were computed after seeing the data. They do not replace anything above.

**1. The generated poses drift under *dynamics*; the anchor's drift is almost
entirely the energy minimisation.** `min.gro` was already on disk, so the
pre-registered window can be split for free:

| group | drift during minimisation (nm) | drift during the 300 ps |
|---|---:|---:|
| A | 0.055 | 0.226 |
| B | 0.048 | 0.159 |
| D | 0.035 | 0.174 |
| V | 0.018 | 0.173 |
| **REF** | 0.052 | **0.049** |

Every group relaxes the same small amount under minimisation. Only REF then
*stays put* once thermal motion is applied.

**2. The drift is almost perfectly one-directional: away from the sulfur.**
Signed Δd was positive in **110 of 111** replicas, spanning +0.036 to +0.762 nm.
The single exception (`t4_85b653525430` rep 3, −0.013 nm) is smaller than that
molecule's own replica spread. **All 37 molecules have a positive mean.** No
docked pose in this cohort tightens its near-attack geometry under dynamics —
the docking places the warhead closer to Cys113 than the force field will hold
it, in every group including REF.

**3. On "is it still a near-attack conformation", V matches REF and the BDHI
groups do not.** Fraction of replicas whose warhead is still inside the
0.28–0.42 nm window at the start of production:

| group | A | B | D | V | REF |
|---|---:|---:|---:|---:|---:|
| fraction still in window | 0.08 | 0.21 | 0.08 | **0.53** | **0.54** |

This is the one place the prereg's *V ≈ REF* statement finds support, and it is
**on a readout the prereg did not name.** On the readout it *did* name (|Δd|), V
sits at 0.203 against REF's 0.102 and is not equivalent. Both are reported;
neither is allowed to stand in for the other.

**4. REF's single worst member is its one SNAr chloroazine** (`xtal:9INP:A1D9X`,
|Δd| = 0.302). Its 5 chloroacetamides sit at a median of 0.090 and its 2
naphthoquinones at 0.093. The anchor is not uniform, and a 5/2/1 class mix at
n = 8 cannot resolve that.

---

## Tier 2 — BPMD

*(pending — see below)*

---

## What this does not settle

Carried forward from the prereg, unchanged, plus what the run added:

- **BDHI-only for the A/B/D contrasts.** BDHI has zero crystallographic
  positives, so a result here does not transfer to other classes.
- **Stability is not reactivity.** A stable near-attack pose is necessary for
  the reaction, not sufficient. Nothing here measures whether a molecule reacts.
- **n = 5 for group V**, descriptive by construction.
- **n = 8 supports only large effects.** The A-vs-B contrast at δ = −0.469 is a
  moderate effect the design cannot resolve either way; "≈" here means "not
  distinguished", not "shown to be equal".
- **The BDHI groups all failed the same way.** With every group's poses leaving
  the window, the contrasts are being drawn between degrees of failure. A
  cohort where some group survived would test the metrics harder than this one
  could.
