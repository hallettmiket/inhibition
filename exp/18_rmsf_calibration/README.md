# exp/18 — is the 2.21 RMSF calibration robust, and is the tolerance per-molecule?

**Verdict: the constant is wide and the tolerance is not per-molecule.** The
conformer ensemble ranks *atoms within* a molecule (ρ = +0.657 — this validates
the weights, and they are fine). It does not rank *molecules against each other*
(ρ = +0.112, CI [−0.06, +0.27]), which is what the tolerance uses. Dividing it by
2.21 does not beat writing one number down (p = 0.515). Rotatable-bond count does.
Record: [D0094](../../decisions/D0094-the-tolerance-was-never-per-molecule-rotatable-bonds-should-set-it.md).

## Run

Needs `pyarrow` for the candidate table, so **`dwi_cheminf`**, not `dwi_admet` —
under `dwi_admet` the descriptor screen is skipped and the section is silently
absent from the report.

```bash
E=/data/lab_vm/envs/dwi_cheminf/bin/python
$E run_all.py            # the constant, both correlations, noise floor, descriptors, convergence
$E tolerance_model.py    # 20x5-fold CV: what should set the tolerance instead
```

Outputs: `00_outputs/blacksmith/rmsf_calibration/`.

## Reading the numbers

* **The two correlations are different questions.** Quote which one you mean. A
  predictor validated on one is not validated on the other.
* **The floor is a flat constant.** Any model that cannot beat one number is not
  a model, and that comparison is what would have caught this on day one.
* **The ceiling is CV 0.24** — the same molecule re-measured on another
  trajectory. Nothing can be scored meaningfully below it.
* The rows are **not independent**: 147 modes from 119 molecules. Bootstraps here
  resample idents.

## Caveats

* 119 molecules, all T_3/T_4 covalent candidates — the descriptor result may not
  transfer to another chemistry.
* Measured RMSF comes from the triage-sweep trajectories, so the calibration is
  tied to that length; RMSF grows with trajectory time.
* `tolerance_from_descriptors` is **proposed, not adopted** — adopting it
  re-groups every cloud (#79).
