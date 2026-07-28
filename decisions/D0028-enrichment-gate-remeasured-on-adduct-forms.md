---
id: D0028
title: Enrichment gate re-measured on adduct forms — D0015's decisive evidence does not survive
date: 2026-07-28
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/run_enrichment_gate.py
  - decisions/D0015-covalent-ranking-uses-affinity-not-cnnaffinity.md
  - data/reference/decoys_covalent_3.csv
  - config/gates.yaml
evidence:
  - 'affinity_kcal on adduct forms: ROC-AUC 0.718, CI [0.483, 0.944] — the CI now INCLUDES 0.5'
  - 'D0015 measured the same metric at ROC-AUC 0.815, CI [0.667, 0.931], which EXCLUDED 0.5'
  - 'cnn_affinity on adduct forms: ROC-AUC 0.392, CI [0.181, 0.645], EF1% 0.0 (was 0.707)'
  - 'affinity EF1% rose 16.7 -> 19.0; BEDROC 0.333 vs cnn 0.146'
  - 'both verdicts remain UNDERPOWERED: 4 independent chemotypes against a floor of 6'
  - 'decoy warhead classes were assigned by the NARROW reactive-atom SMARTS; [CH2][Cl] matched nitrogen mustards (cyclophosphamide) and nitrosoureas (lomustine) as chloroacetamides'
  - 'only 112 of 294 decoys carry a whole warhead group AND survive the adduct transform'
  - 'verified decoys are 104 acrylamide, 5 naphthoquinone_c2, 3 chloroacetamide'
  - 'two of six actives (sulfamate_acetamide) and one (snar_chloroazine) have NO same-class decoy'
runbook: null
---

## Context

D0022 changed what gets docked: the adduct form rather than the pre-reaction
ligand. That record noted the enrichment gate had been measured through the old
protocol and said it "should be re-measured rather than assumed", since actives
and decoys were treated identically and ROC-AUC was unlikely to move much.

That prediction was wrong in the way that matters.

## What the re-measurement found

**The gate is weaker on adduct forms, and D0015's decisive evidence no longer
holds.**

| metric | D0015 (pre-reaction) | now (adduct) |
|---|---|---|
| `affinity_kcal` ROC-AUC | 0.815 | **0.718** |
| its 95% CI | [0.667, 0.931] — excludes 0.5 | **[0.483, 0.944] — includes 0.5** |
| `affinity_kcal` EF1% | 16.7 | 19.0 |
| `cnn_affinity` ROC-AUC | 0.707 | 0.392 |

D0015 chose `affinity_kcal` over `CNNaffinity` on the explicit grounds that
affinity's interval excluded 0.5 while CNNaffinity's did not. **On adduct forms
neither interval excludes 0.5.** The comparison that decided the rank metric is
gone.

The *direction* survives and is if anything sharper: affinity beats CNNaffinity
on every statistic (AUC 0.718 vs 0.392, EF1% 19.0 vs 0.0, BEDROC 0.333 vs
0.146), and CNNaffinity is now worse than random at the point estimate. gnina's
own warning that CNN scoring is uncalibrated for covalent docking stands.

## A second, larger problem: the decoy set

Re-running exposed a defect in the decoys themselves, not in the docking.

Decoy warhead classes were assigned by the **narrow reactive-atom SMARTS**. That
pattern is written to tell gnina which atom to bond, not to identify a
chemotype, and it is far too permissive for the latter: `[CH2][Cl]` matches any
primary alkyl chloride, so **cyclophosphamide (a nitrogen mustard) and lomustine
(a nitrosourea) were both labelled `chloroacetamide`**. Their 2-chloroethyl
groups are genuinely electrophilic but are a different warhead class entirely.

Requiring the *whole* warhead group and a successful adduct transform leaves
**112 of 294** decoys, distributed 104 acrylamide / 5 naphthoquinone_c2 /
3 chloroacetamide. Against six actives spanning four classes, that means:

- `sulfamate_acetamide` — 2 actives, **0** decoys
- `snar_chloroazine` — 1 active, **0** decoys
- `acrylamide` — **0** actives, 104 decoys

**D0020 established that `affinity_kcal` is not comparable across warhead
classes.** A gate that scores actives of one class against decoys of another is
therefore measuring exactly the quantity D0020 says is meaningless. Some unknown
share of the separation — in the old measurement and the new one — is chemotype
separation rather than binding discrimination.

## Decision

**D0015's conclusion is retained; its stated justification is not.**
`affinity_kcal` remains T_3 and T_4's rank metric, because it beats the
alternative on every statistic and because CNNaffinity is uncalibrated for
covalent docking by the tool's own admission. But the record must not keep
citing an interval that excludes 0.5, because it no longer does.

**The covalent gate's verdict stays UNDERPOWERED**, which is what it has always
said. The graded-verdict floor (6 independent chemotypes, gates.yaml) refused to
claim more than UNDERPOWERED even when the point estimates looked good, and that
refusal is now vindicated: the flattering point estimate did not survive a
change in the ligand form.

`decoys_covalent_3.csv` records the verified class assignment for future runs.

## Consequences

**T_4's ranking is weakly supported and should be presented that way.** It was
already displayed with its gate verdict; the verdict has not changed, only the
confidence one should place in the numbers behind it.

**The decoy set needs regenerating before the gate can say anything stronger.**
Property matching alone is not enough — decoys must be matched on warhead class
too, or the gate cannot separate binding discrimination from chemotype
discrimination. That is D0014's territory and is not attempted here.

**Assigning a chemotype with a reactive-atom SMARTS is a category error** and it
appeared twice in this project: here, and in the alert attribution, where
excusing only the reactive atoms let `alpha_halo_carbonyl` straddle the boundary
(D0025). The reactive-atom pattern says where a bond forms. The whole-warhead
fragment says what the chemistry is. They are not interchangeable and the
library carries both for that reason.
