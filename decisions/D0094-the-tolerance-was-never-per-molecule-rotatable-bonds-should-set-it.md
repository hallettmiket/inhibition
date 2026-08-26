---
id: D0094
title: The splitting tolerance was never per-molecule — the ensemble ranks atoms, not molecules, and rotatable-bond count should set it
date: 2026-08-26
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/pose_contacts.py
  - exp/15_rmsf_predictor/run_all.py
  - exp/16_contact_clustering/run_all.py
  - exp/18_rmsf_calibration/run_all.py
  - exp/18_rmsf_calibration/tolerance_model.py
  - decisions/D0092-contact-space-is-fixed-the-group-count-climbs-because-6000-poses-undersample-it.md
evidence:
  - '@tt8804: "if we are callibrating to one experiental constant you better make sure we ran a large enough and robust experiment"'
  - 'THE CONSTANT: 2.207 with a cluster-bootstrap 95% CI of [1.95, 2.51]; per-molecule ratios span 0.90-6.80 and only 35% of molecules sit within +-25% of it'
  - 'THE ROWS WERE NOT INDEPENDENT: 147 swept modes come from 119 molecules; 28 rows share an ident, same prediction, different trajectory. A row-wise bootstrap would have understated the interval'
  - 'THE TWO CORRELATIONS ARE DIFFERENT NUMBERS: within a molecule across atoms rho = +0.657 (validates the WEIGHTS); across molecules in absolute scale rho = +0.112, 95% CI [-0.06, +0.27], CROSSING ZERO (what the TOLERANCE uses)'
  - 'the prediction barely moves between molecules: CV 0.15 against the measurement CV 0.45 -- the truth varies 3.0x more than the predictor does'
  - 'IT DOES NOT BEAT ONE NUMBER: median relative error 31.2% for ensemble/2.21 against 33.8% for a flat 0.61 A, Wilcoxon signed-rank p = 0.515'
  - 'ROTATABLE BONDS DO: rho = +0.523 against the ensemble +0.124, over the same 119 molecules, from a descriptor costing nothing'
  - 'OUT OF SAMPLE (20x5-fold CV grouped by ident, 119 molecules): rotatable bonds + ensemble 24.9%, rotatable bonds + heavy 25.9%, rotatable bonds 26.2%, shipped 32.2%, flat constant 33.1%'
  - 'THE CEILING: the same molecule re-measured on a different trajectory moves at CV 0.24, so 24.9% is AT the measurement own reproducibility -- no model can be scored below it on this data'
  - 'THE PREDICTOR IS NOT BROKEN, it is converged and stable: 50 conformers is within 1.4-2.1% of 100 and 200, seed-to-seed 3.5%. It measures the wrong quantity for this purpose, precisely and repeatably'
runbook: python exp/18_rmsf_calibration/run_all.py; python exp/18_rmsf_calibration/tolerance_model.py
---

# D0094 — the tolerance was never per-molecule

## What was claimed

`pose_contacts.py` said, and exp/16 and D0092 repeated: *"THE SCALE IS THE
MOLECULE'S OWN … a floppy molecule gets a looser cut than a rigid one, and nobody
picks a number."* The tolerance is `median(ensemble RMSF) / 2.21`, and `2.21` is
the one experimental constant the whole grouping rests on.

Nobody picks a number. But the number that emerges carries almost no
molecule-specific information.

## What the constant actually was

One line in exp/15's summary block — `d.pred_med.median() / d.meas_med.median()`
— printed to two decimals, with no interval, no stratification, and no check that
the quantity it calibrates is the quantity it was validated on. Audited:

| | |
|---|---|
| ratio of medians (shipped) | 2.207, 95% CI **[1.95, 2.51]** |
| median of ratios | 2.122, 95% CI [1.91, 2.44] |
| per-molecule ratio | range **0.90 – 6.80**, IQR [1.57, 3.13] |
| molecules within ±25% of it | **35%** |

The interval alone would not condemn it. The next section does.

## The two correlations are not the same number, and only one was measured

exp/15 reports **ρ = 0.657**. That is *within* a molecule, *across* its atoms:
"does the ensemble rank this molecule's floppy atoms correctly?" It validates the
per-atom **weights**, which is what it was built for, and the weights are fine.

The **tolerance** uses a different quantity — the absolute scale of one molecule
against another. Measured for the first time here:

* across molecules, **ρ = +0.112, 95% CI [−0.06, +0.27]** — crossing zero
* predicted RMSF varies at **CV 0.15**; measured varies at **CV 0.45**

The predictor is nearly constant across molecules while the truth is not. So
`median(rmsf)/2.21` is approximately one number plus noise, and the noise is
uncorrelated with the thing it is meant to track.

**Both quantities are "the RMSF predictor works". Both are populated and
plausible. One was measured and the other was assumed.** That is disguise #1 in
its purest form — a value selected because its name described what was wanted.

## It does not beat writing one number down

| | median relative error |
|---|---|
| `ensemble / 2.21` | 31.2% |
| a flat 0.61 Å for every molecule | 33.8% |

Wilcoxon signed-rank **p = 0.515**. Two and a half points, indistinguishable from
nothing. Whatever the derivation looks like, the shipped tolerance is a constant.

## Rotatable-bond count does beat it

Out of sample, 20×5-fold cross-validation grouped by ident, 119 molecules:

| model | out-of-sample median error |
|---|---|
| **rotatable bonds + ensemble** | **24.9%** |
| rotatable bonds + heavy atoms | 25.9% |
| rotatable bonds | 26.2% |
| ensemble / k (shipped) | 32.2% |
| flat constant (the floor) | 33.1% |

And the ceiling matters as much as the floor: the same molecule re-measured on a
different trajectory moves at **CV 0.24**. The combined model at 24.9% is *at the
measurement's own reproducibility* — as good as this data can support. The
ensemble is not useless; it adds 1.3 points on top of rotatable bonds. It is just
not sufficient alone, and it was carrying the tolerance alone.

## What is NOT wrong

* **The atom weights.** ρ = 0.657 within a molecule is the right validation for
  them, and it holds.
* **The predictor's stability.** 50 conformers sits within 1.4–2.1% of 100 and
  200; seed-to-seed movement is 3.5%. It is converged. It measures the wrong
  quantity for this purpose, precisely and repeatably — which is why more
  conformers would not have helped and why nothing looked broken.
* **D0092's conclusions.** The saturation result swept tolerances from 0.73 Å to
  3.5 Å and found no plateau at any of them, so it does not depend on which one is
  chosen. The persistence and fixed-region results likewise hold for any fixed
  absolute cut.

## What this decides

1. **Correct the claim in the code and the documents.** Done: `pose_contacts.py`
   now states which correlation validates what.
2. **`tolerance_from_descriptors` is built and NOT adopted**, alongside
   `pose_cluster.py`'s precedent. Adopting it re-groups every cloud and would make
   every number measured under the current tolerance stale; it belongs with the
   re-screen (#79).
3. **The coefficients are deliberately not pinned in source.** A fitted constant
   written into a module is the stale-pin family this project has hit five times;
   the function raises unless the caller passes a fit.
4. **Report the tolerance with its provenance.** Any artefact quoting a tolerance
   should say whether it came from the ensemble (a near-constant) or from a fitted
   descriptor model.

## Why it looked right

Because the predictor works — at the job it was validated for. exp/15 is careful
work: it defeats a bond-perception trap, a coordinate-matching trap, a
non-bijective assignment and a kekulization failure, and it reports an honest
ρ = 0.657. Nothing in it is wrong. The error was downstream, in reading "the RMSF
predictor is validated" as licence for a second use it was never tested on — and
the second use produced tolerances of 0.48–0.82 Å that look exactly like
per-molecule numbers, because they are per-molecule numbers. They are just not
*informative* per-molecule numbers.

## Guard

`exp/18_rmsf_calibration/run_all.py` prints the two correlations side by side,
labelled with what each validates, and scores the predictor against a flat
constant. A predictor that cannot beat one number is not a predictor, and that
comparison is the check that would have failed on day one.
