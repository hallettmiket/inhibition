#!/usr/bin/env python3
"""
Purpose: is the 2.21 RMSF calibration robust, and is the tolerance really per-molecule?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: exp/15's per-molecule predicted vs measured RMSF table; the candidate table
Output: 00_outputs/blacksmith/rmsf_calibration/

`pose_contacts.RMSF_CALIBRATION = 2.21` is the ONE experimental constant the whole
grouping rests on, and @tt8804 asked whether the experiment behind it was large
and robust enough. It was a single point estimate -- `pred_med.median() /
meas_med.median()` in exp/15's summary block, printed to two decimals, with no
interval, no stratification, and no check that the quantity it calibrates is the
quantity it was validated on.

THE TWO CORRELATIONS ARE NOT THE SAME NUMBER, AND ONLY ONE WAS MEASURED. exp/15
reports rho = 0.657: WITHIN a molecule, ACROSS its atoms -- "does the ensemble
rank this molecule's floppy atoms correctly". That validates the per-atom
WEIGHTS, which is what it was built for. The TOLERANCE uses a different quantity:
the ABSOLUTE scale of one molecule against another. Nothing measured that. Both
are "the RMSF predictor works", both are populated and plausible, and the second
is the one `median(rmsf)/2.21` actually depends on.

THE ROWS ARE NOT INDEPENDENT. 147 swept modes come from 119 molecules, so several
rows share an ident -- same molecule, different pose, therefore IDENTICAL
prediction and a different measurement. A naive bootstrap over rows would
understate the interval; the bootstrap here resamples IDENTS. The repeats are
also a free control: they measure how much `meas_med` moves for one molecule
between trajectories, which is the noise floor no predictor can beat.

WHAT WOULD MAKE THIS PASS WHEN IT SHOULD FAIL: comparing the predictor only
against itself. So it is compared against a FLAT CONSTANT -- if dividing a
near-constant prediction by 2.21 is no better than writing down one number for
every molecule, then "the tolerance is the molecule's own" is not true, however
principled the derivation looks.
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
from shared import pose_contacts as pc              # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("rmsf-calibration")


def exp15_table() -> pd.DataFrame:
    fs = sorted(glob.glob(str(rp.BLACKSMITH / "rmsf_predictor/rmsf_predictor_*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit("run exp/15_rmsf_predictor first")
    best = max((pd.read_csv(f) for f in fs), key=len)
    return best[(best.pred_med > 0) & (best.meas_med > 0)].copy()


def cluster_bootstrap(d: pd.DataFrame, stat, n: int = 4000, seed: int = 7):
    """Resample IDENTS, not rows -- rows sharing an ident are one observation."""
    rng = np.random.default_rng(seed)
    groups = [g for _, g in d.groupby("ident")]
    out = []
    for _ in range(n):
        pick = rng.integers(0, len(groups), len(groups))
        out.append(stat(pd.concat([groups[i] for i in pick], ignore_index=True)))
    return np.array(out)


def descriptors(smiles: str) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {}
    return dict(
        heavy=m.GetNumHeavyAtoms(),
        rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
        mw=Descriptors.MolWt(m),
        tpsa=rdMolDescriptors.CalcTPSA(m),
        rings=rdMolDescriptors.CalcNumRings(m),
        fsp3=rdMolDescriptors.CalcFractionCSP3(m),
        logp=Descriptors.MolLogP(m),
        # rotatable bonds per heavy atom: flexibility DENSITY, which is closer to
        # what a per-atom fluctuation should scale with than a raw count
        rot_density=rdMolDescriptors.CalcNumRotatableBonds(m) / max(m.GetNumHeavyAtoms(), 1),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--boot", type=int, default=4000)
    ap.add_argument("--conv-molecules", type=int, default=18)
    ap.add_argument("--conv-levels", default="10,25,50,100,200")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from scipy.stats import spearmanr, pearsonr

    d = exp15_table()
    d["ratio"] = d.pred_med / d.meas_med
    log.info("%d rows from %d molecules", len(d), d.ident.nunique())
    t = sout.Topic("blacksmith", "rmsf_calibration")

    # ---------- 1. the constant itself ---------------------------------- #
    rom = lambda x: float(x.pred_med.median() / x.meas_med.median())   # noqa: E731
    mor = lambda x: float((x.pred_med / x.meas_med).median())          # noqa: E731
    bs_rom = cluster_bootstrap(d, rom, a.boot, a.seed)
    bs_mor = cluster_bootstrap(d, mor, a.boot, a.seed)

    # ---------- 2. within vs across molecule ---------------------------- #
    rho_across = spearmanr(d.pred_med, d.meas_med)[0]
    bs_rho = cluster_bootstrap(
        d, lambda x: float(spearmanr(x.pred_med, x.meas_med)[0]), 1500, a.seed)
    rho_within = float(d.spearman.median())

    # ---------- 3. the noise floor, from repeated idents ----------------- #
    rep = d.groupby("ident").filter(lambda g: len(g) > 1)
    if len(rep):
        spread = rep.groupby("ident").meas_med.agg(
            lambda v: (v.max() - v.min()) / v.median())
        within_mol_cv = rep.groupby("ident").meas_med.agg(
            lambda v: v.std(ddof=1) / v.mean()).median()
    else:
        spread, within_mol_cv = pd.Series(dtype=float), float("nan")

    # ---------- 4. does it beat a flat constant? ------------------------ #
    k = rom(d)
    flat = float(d.meas_med.median())
    err_pred = np.abs(d.pred_med / k - d.meas_med) / d.meas_med
    err_flat = np.abs(flat - d.meas_med) / d.meas_med
    from scipy.stats import wilcoxon
    try:
        p_beat = float(wilcoxon(err_pred, err_flat).pvalue)
    except ValueError:
        p_beat = float("nan")

    # ---------- 5. does any cheap descriptor do better? ----------------- #
    desc_rows = []
    try:
        import nac_rank as nr
        smi = {c.ident: c.smiles for c in nr.load_candidates()}
        got = [(i, smi[i]) for i in d.ident.unique() if i in smi]
        log.info("descriptors for %d of %d molecules", len(got), d.ident.nunique())
        D = pd.DataFrame([dict(ident=i, **descriptors(s)) for i, s in got]).dropna()
        md = d.groupby("ident").agg(meas=("meas_med", "median"),
                                    pred=("pred_med", "median")).reset_index()
        D = D.merge(md, on="ident")
        for c in [c for c in D.columns if c not in ("ident", "meas", "pred")]:
            desc_rows.append(dict(descriptor=c, n=len(D),
                                  rho=float(spearmanr(D[c], D.meas)[0])))
        desc_rows.append(dict(descriptor="ensemble RMSF (shipped)", n=len(D),
                              rho=float(spearmanr(D.pred, D.meas)[0])))
    except Exception as exc:                                    # noqa: BLE001
        log.warning("descriptor screen skipped: %s", str(exc)[:90])
        D = pd.DataFrame()
    desc = pd.DataFrame(desc_rows).sort_values(
        "rho", key=lambda s: s.abs(), ascending=False) if desc_rows else pd.DataFrame()

    # ---------- 6. has the predictor converged at 50 conformers? -------- #
    conv = pd.DataFrame()
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "e15", REPO / "exp" / "15_rmsf_predictor" / "run_all.py")
        e15 = importlib.util.module_from_spec(spec)
        sys.modules["e15"] = e15
        spec.loader.exec_module(e15)
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        levels = [int(x) for x in a.conv_levels.split(",")]
        rng = np.random.default_rng(a.seed)
        cands = [i for i in d.ident.unique()
                 if (rp.allposes_dir() / f"{i}.sdf").is_file()]
        pick = rng.choice(cands, size=min(a.conv_molecules, len(cands)), replace=False)
        rows = []
        for n_, ident in enumerate(pick, 1):
            ms = Chem.SDMolSupplier(str(rp.allposes_dir() / f"{ident}.sdf"),
                                    removeHs=False, sanitize=True)
            tmpl = next((m for m in ms if m is not None), None)
            if tmpl is None:
                continue
            hv = [x.GetIdx() for x in tmpl.GetAtoms() if x.GetAtomicNum() > 1]
            for L in levels:
                for s in (a.seed, a.seed + 101):
                    try:
                        r = e15.predict_rmsf(tmpl, hv, L, s)
                    except Exception:                            # noqa: BLE001
                        continue
                    rows.append(dict(ident=ident, n_conf=L, seed=s,
                                     pred_med=float(np.median(r))))
            log.info("  convergence %2d/%d %s", n_, len(pick), ident)
        conv = pd.DataFrame(rows)
    except Exception as exc:                                    # noqa: BLE001
        log.warning("convergence sweep skipped: %s", str(exc)[:90])

    # ---------- write ---------------------------------------------------- #
    d.to_csv(t.write("per_molecule", ".csv"), index=False)
    if len(desc):
        desc.to_csv(t.write("descriptor_screen", ".csv"), index=False)
    if len(conv):
        conv.to_csv(t.write("conformer_convergence", ".csv"), index=False)
    pd.DataFrame(dict(ratio_of_medians=bs_rom, median_of_ratios=bs_mor)).to_csv(
        t.write("bootstrap", ".csv"), index=False)

    # ---------- report ---------------------------------------------------- #
    P = print
    P("\n" + "=" * 80)
    P("  IS 2.21 ROBUST, AND IS THE TOLERANCE REALLY PER-MOLECULE?")
    P("=" * 80)
    P(f"\n  {len(d)} swept modes from {d.ident.nunique()} molecules "
      f"({len(d) - d.ident.nunique()} repeated idents — not independent rows)")

    P("\n  1. THE CONSTANT")
    P(f"     ratio of medians (shipped)  {rom(d):.3f}   "
      f"95% CI [{np.percentile(bs_rom, 2.5):.2f}, {np.percentile(bs_rom, 97.5):.2f}]")
    P(f"     median of ratios            {mor(d):.3f}   "
      f"95% CI [{np.percentile(bs_mor, 2.5):.2f}, {np.percentile(bs_mor, 97.5):.2f}]")
    P(f"     mean of ratios              {d.ratio.mean():.3f}")
    P(f"     per-molecule ratio: IQR [{d.ratio.quantile(.25):.2f}, "
      f"{d.ratio.quantile(.75):.2f}], full range [{d.ratio.min():.2f}, {d.ratio.max():.2f}]")
    P(f"     molecules within +-25% of the constant: "
      f"{((d.ratio > .75 * k) & (d.ratio < 1.25 * k)).mean() * 100:.0f}%")

    P("\n  2. THE TWO CORRELATIONS — only the first was ever measured")
    P(f"     WITHIN a molecule, across atoms (validates the WEIGHTS)   "
      f"rho = {rho_within:+.3f}")
    P(f"     ACROSS molecules, absolute scale (what the TOLERANCE uses) "
      f"rho = {rho_across:+.3f}"
      f"   95% CI [{np.percentile(bs_rho, 2.5):+.2f}, {np.percentile(bs_rho, 97.5):+.2f}]")
    P(f"     spread across molecules: predicted CV {d.pred_med.std() / d.pred_med.mean():.2f}, "
      f"measured CV {d.meas_med.std() / d.meas_med.mean():.2f}")
    P("     -> the prediction barely varies between molecules; the truth varies "
      f"{d.meas_med.std() / d.meas_med.mean() / (d.pred_med.std() / d.pred_med.mean()):.1f}x more")

    if len(spread):
        P(f"\n  3. NOISE FLOOR — {len(spread)} molecules measured more than once")
        P(f"     same molecule, different trajectory: median CV of measured RMSF "
          f"{within_mol_cv:.2f}")
        P(f"     median spread (max-min)/median: {spread.median():.2f}")

    P("\n  4. DOES IT BEAT WRITING DOWN ONE NUMBER?")
    P(f"     predicted / {k:.2f}          median relative error {np.median(err_pred) * 100:5.1f}%")
    P(f"     a flat {flat:.2f} A for every molecule   median relative error "
      f"{np.median(err_flat) * 100:5.1f}%")
    P(f"     Wilcoxon signed-rank p = {p_beat:.3f}   -> "
      f"{'the predictor beats a constant' if p_beat < 0.05 and np.median(err_pred) < np.median(err_flat) else 'NO significant improvement over a flat constant'}")

    if len(desc):
        P("\n  5. DOES ANY CHEAP DESCRIPTOR PREDICT MEASURED RMSF BETTER?")
        for _, r in desc.iterrows():
            P(f"     {r.descriptor:<26} rho = {r.rho:+.3f}  (n={int(r.n)})")

    if len(conv):
        P("\n  6. HAS THE PREDICTOR CONVERGED AT 50 CONFORMERS?")
        g = conv.groupby("n_conf").pred_med.agg(["median", "std"])
        base = conv[conv.n_conf == 50].groupby("ident").pred_med.median()
        for L, r in g.iterrows():
            at = conv[conv.n_conf == L].groupby("ident").pred_med.median()
            j = at.index.intersection(base.index)
            drift = float(np.median(np.abs(at[j] - base[j]) / base[j]) * 100) if len(j) else float("nan")
            P(f"     n_conf={L:4d}: median pred RMSF {r['median']:.3f} A   "
              f"median |change| vs 50 conformers {drift:5.1f}%")
        s1 = conv[conv.seed == a.seed].groupby(["ident", "n_conf"]).pred_med.median()
        s2 = conv[conv.seed != a.seed].groupby(["ident", "n_conf"]).pred_med.median()
        j = s1.index.intersection(s2.index)
        if len(j):
            P(f"     seed-to-seed: median |change| "
              f"{float(np.median(np.abs(s1[j] - s2[j]) / s1[j]) * 100):.1f}%")

    P("\n" + "=" * 80)
    P(f"  written to {t.dir}")
    P("=" * 80 + "\n")


if __name__ == "__main__":
    main()
