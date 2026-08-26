#!/usr/bin/env python3
"""
Purpose: what should set the splitting tolerance, if the conformer ensemble does not?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: exp/15's predicted-vs-measured table + the candidate SMILES
Output: 00_outputs/blacksmith/rmsf_calibration/

run_all.py established that `median(ensemble RMSF)/2.21` does NOT beat writing one
number down for every molecule (Wilcoxon p = 0.515), because the ensemble barely
varies between molecules (CV 0.15) while the truth varies three times as much.
It also found that ROTATABLE BOND COUNT ranks molecules at rho = +0.523 where the
ensemble manages +0.124 -- four times the signal, from a descriptor that costs
nothing.

A CORRELATION IS NOT A MODEL, so this fits the candidates and scores them
OUT OF SAMPLE. rho = 0.523 measured on the same 119 molecules a model is fitted
to would be the leakage this project has already ruled Boltz-2 out for; grouped
5-fold cross-validation by IDENT is the honest version.

THE FLOOR AND THE CEILING ARE BOTH REPORTED. The floor is a flat constant -- any
model that cannot beat one number is not a model. The ceiling is the measurement's
own reproducibility: the same molecule re-measured on a different trajectory moves
by CV 0.24, so no predictor can be scored below roughly that and the comparison is
meaningless without it on screen.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                  # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("tolerance-model")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--repeats", type=int, default=20)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from scipy.stats import spearmanr
    import nac_rank as nr
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    fs = sorted(glob.glob(str(rp.BLACKSMITH / "rmsf_predictor/rmsf_predictor_*.csv")),
                key=os.path.getmtime)
    d = max((pd.read_csv(f) for f in fs), key=len)
    d = d[(d.pred_med > 0) & (d.meas_med > 0)]
    # ONE ROW PER MOLECULE. 147 modes come from 119 molecules and the prediction is
    # identical within a molecule, so scoring per row would weight some molecules
    # more than others for no reason.
    g = d.groupby("ident").agg(meas=("meas_med", "median"),
                               ens=("pred_med", "median")).reset_index()

    smi = {c.ident: c.smiles for c in nr.load_candidates()}
    rows = []
    for r in g.itertuples():
        s = smi.get(r.ident)
        m = Chem.MolFromSmiles(s) if s else None
        if m is None:
            continue
        rows.append(dict(ident=r.ident, meas=r.meas, ens=r.ens,
                         rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
                         heavy=m.GetNumHeavyAtoms()))
    D = pd.DataFrame(rows)
    log.info("%d molecules with SMILES and a measurement", len(D))

    MODELS = {
        "flat constant (the floor)":        [],
        "ensemble RMSF / k (shipped)":      ["ens"],
        "rotatable bonds":                  ["rotb"],
        "rotatable bonds + heavy atoms":    ["rotb", "heavy"],
        "rotatable bonds + ensemble":       ["rotb", "ens"],
    }

    def fit_predict(tr, te, cols):
        """Least squares on log(meas) -- RMSF is positive and right-skewed, and a
        linear fit on the raw scale would let one floppy molecule set the line."""
        y = np.log(tr.meas.values)
        if not cols:
            return np.full(len(te), float(np.exp(y.mean())))
        X = np.column_stack([np.ones(len(tr))] + [tr[c].values.astype(float) for c in cols])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        Xt = np.column_stack([np.ones(len(te))] + [te[c].values.astype(float) for c in cols])
        return np.exp(Xt @ beta)

    rng = np.random.default_rng(a.seed)
    scores = {k: [] for k in MODELS}
    for rep in range(a.repeats):
        order = rng.permutation(len(D))
        folds = np.array_split(order, a.folds)
        for f in folds:
            te = D.iloc[f]
            tr = D.iloc[np.setdiff1d(order, f)]
            for name, cols in MODELS.items():
                p = fit_predict(tr, te, cols)
                scores[name].append(np.median(np.abs(p - te.meas.values) / te.meas.values))

    res = pd.DataFrame([dict(model=k, err=float(np.mean(v) * 100),
                             sd=float(np.std(v) * 100)) for k, v in scores.items()])
    res = res.sort_values("err")
    t = sout.Topic("blacksmith", "rmsf_calibration")
    res.to_csv(t.write("tolerance_model_cv", ".csv"), index=False)
    D.to_csv(t.write("tolerance_model_data", ".csv"), index=False)

    noise = 0.24   # measured in run_all.py section 3
    P = print
    P("\n" + "=" * 78)
    P("  WHAT SHOULD SET THE TOLERANCE?  "
      f"{a.repeats}x{a.folds}-fold CV, {len(D)} molecules")
    P("=" * 78)
    P(f"\n  out-of-sample median relative error in predicting a molecule's RMSF\n")
    P(f"    {'model':<34}{'error':>9}{'  vs flat':>11}")
    flat = float(res[res.model.str.startswith("flat")].err.iloc[0])
    for _, r in res.iterrows():
        mark = "  <- shipped" if "shipped" in r.model else ""
        P(f"    {r.model:<34}{r.err:8.1f}% {(r.err - flat):+9.1f}pp{mark}")
    P(f"\n    measurement's own reproducibility (the ceiling): ~{noise * 100:.0f}%")
    P("    a model cannot be scored meaningfully below that.")
    best = res.iloc[0]
    P(f"\n  BEST: {best.model} at {best.err:.1f}%")
    P(f"  Spearman(rotatable bonds, measured RMSF) = "
      f"{spearmanr(D.rotb, D.meas)[0]:+.3f}   "
      f"ensemble = {spearmanr(D.ens, D.meas)[0]:+.3f}")
    P("\n" + "=" * 78)
    P(f"  written to {t.dir}")
    P("=" * 78 + "\n")


if __name__ == "__main__":
    main()
