#!/usr/bin/env python3
"""
Purpose: is the flat energy landscape the pocket's, or our reactive receptor's?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: --candidates (ids), docked into BOTH the reactive and the plain 3IKD
Output: 00_outputs/blacksmith/flat_landscape_cause/

THE CHECK D0090 FLAGGED AS DECISIVE. D0090 measured that a molecule's whole pose
cloud spans a median 3.96 kcal/mol -- inside the scoring function's own 2-3
kcal/mol error -- and argued the resulting flat landscape is why the cloud never
saturates: a search with no basins to fall into ends somewhere new every run.

But the screen docks into a REACTIVE receptor whose van der Waals parameters were
softened deliberately (`R_EQ_12 = 3.2`, `EPS_12 = 1.0`, `nac_screen`) so the
warhead can approach Cys113. Softening the potential flattens the landscape near
the warhead BY CONSTRUCTION. So the finding may be a property of our setup rather
than of the pocket, and D0090 says so and stays `proposed` until this runs.

WHAT IS COMPARED, AND WHAT IS NOT. Both arms dock the same molecule at the same
depth and are scored by the same AutoDock function; only the receptor
parameterisation differs -- softened + flexible Cys113 versus stock maps. That is
the variable of interest.

It is NOT a clean single-variable experiment and should not be reported as one:
the reactive arm also carries reactive atom typing and a flexible sidechain,
because those are not separable from the reactive setup as the screen builds it.
The claim this can support is "the reactive SETUP flattens the landscape", not
"the softened vdW term alone does".

READOUTS. Energy spread within the cloud (the D0090 quantity), and the covering
number's saturation exponent at the pipeline's own 3.5 A tolerance. If the plain
receptor gives a wider spread AND a smaller exponent, the flatness is ours.
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                  # noqa: E402
from shared import run_paths as rp                  # noqa: E402
from shared import target_config as tc              # noqa: E402
import nac_screen as ns                             # noqa: E402
import nac_rank as nr                               # noqa: E402

log = logging.getLogger("flat-cause")
LADDER = (100, 250, 500, 1000, 2000)


def cover(coords: np.ndarray, r: float) -> int:
    n = len(coords)
    if n == 0:
        return 0
    dmin = np.full(n, np.inf)
    chosen, nxt = 0, 0
    while True:
        d = np.sqrt(((coords - coords[nxt]) ** 2).sum(axis=2).mean(axis=1))
        dmin = np.minimum(dmin, d)
        chosen += 1
        far = int(np.argmax(dmin))
        if dmin[far] <= r or chosen >= n:
            return chosen
        nxt = far


def _rebuild(dlg: Path):
    """(heavy-atom coords per pose, energies) from a .dlg, no SMARTS needed."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    pm = PDBQTMolecule.from_file(str(dlg), is_dlg=True, skip_typing=True)
    mols = [m for m in RDKitMolCreate.from_pdbqt_mol(pm) if m is not None]
    if not mols:
        raise ValueError(f"nothing rebuilt from {dlg}")
    mol = mols[0]
    idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    xyz = np.array([mol.GetConformer(c).GetPositions()[idx]
                    for c in range(mol.GetNumConformers())])
    en = np.array([e for e in ns.pose_energies(dlg)], dtype=float)
    return xyz, en


def dock_reactive(cand, nrun: int, gpu: str, seed: int):
    rec_dir = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    w = Path(tempfile.mkdtemp(prefix="fc_rx_"))
    ligs = list(ns.prepare_ligand(cand, w / "lig.pdbqt"))
    if not ligs:
        raise RuntimeError("ligand preparation produced nothing")
    dlg = ns.dock(ligs[0], rec_dir, w / "c0", nrun, gpu, seed=seed)
    return _rebuild(dlg)


def dock_plain(cand, nrun: int, gpu: str, seed: int):
    """Stock maps, stock ligand -- the same call `enrichment_gate_3ikd` makes."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    RDLogger.DisableLog("rdApp.*")
    plain = rp.receptor_plain()
    w = Path(tempfile.mkdtemp(prefix="fc_pl_"))
    mol = ns.largest_fragment(cand.smiles)
    if mol is None:
        raise ValueError("unparseable smiles")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
        raise ValueError("embed failed")
    AllChem.MMFFOptimizeMolecule(mol)
    txt, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation()(mol)[0])
    if not ok:
        raise ValueError(err)
    lig = w / "lig.pdbqt"
    lig.write_text(txt)
    cmd = [str(ns.AUTODOCK), "-M", "rec.maps.fld", "-L", str(lig),
           "--nrun", str(nrun), "--resnam", str((w / "out").resolve())]
    if seed is not None:
        cmd += ["--seed", str(int(seed))]
    r = subprocess.run(cmd, cwd=plain, capture_output=True, text=True,
                       env=dict(os.environ, CUDA_VISIBLE_DEVICES=gpu))
    dlg = w / "out.dlg"
    out = (r.stdout or "") + (r.stderr or "")
    if r.returncode != 0 or not dlg.is_file() or "not successful" in out:
        raise RuntimeError(f"plain dock failed (exit {r.returncode}): "
                           + "\n".join(out.strip().splitlines()[-3:]))
    return _rebuild(dlg)


def fit_b(n, y):
    ok = (np.asarray(y) > 0) & (np.asarray(n) > 0)
    if ok.sum() < 3:
        return float("nan")
    ln, ly = np.log(np.asarray(n, float)[ok]), np.log(np.asarray(y, float)[ok])
    A = np.vstack([np.ones_like(ln), ln]).T
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    return float(coef[1])


def summarise(xyz, en, tol, seed):
    en = en[~np.isnan(en)]
    best = float(en.min()) if len(en) else float("nan")
    rng = np.random.default_rng(seed)
    ladder = sorted(set([k for k in LADDER if k <= len(xyz)] + [len(xyz)]))
    cov = [cover(xyz[rng.choice(len(xyz), size=k, replace=False)], tol)
           for k in ladder]
    return dict(
        n_poses=len(xyz), best=round(best, 2),
        median=round(float(np.median(en)), 2),
        span=round(float(en.max() - en.min()), 2),
        within_1=round(float((en <= best + 1.0).mean()), 3),
        within_2=round(float((en <= best + 2.0).mean()), 3),
        b_tol=round(fit_b(ladder, cov), 3),
        centres_at_500=(cov[ladder.index(500)] if 500 in ladder else None))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--nrun", type=int, default=2000)
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--seed", type=int, default=13)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    tol = float(tc.get("md.sweep_survivor_rmsd_nm")) * 10.0
    cands = {c.ident: c for c in nr.load_candidates()}
    rows = []
    for ident in a.candidates:
        if ident not in cands:
            log.warning("%s not in the candidate table", ident)
            continue
        for arm, fn in (("reactive", dock_reactive), ("plain", dock_plain)):
            try:
                xyz, en = fn(cands[ident], a.nrun, a.gpu, a.seed)
            except Exception as exc:                               # noqa: BLE001
                log.warning("  %s / %s failed: %s", ident, arm, exc)
                continue
            s = summarise(xyz, en, tol, a.seed)
            rows.append(dict(ident=ident, arm=arm, **s))
            log.info("  %s %-8s n=%d span=%.2f within2=%.0f%% b@%.1fA=%.3f",
                     ident, arm, s["n_poses"], s["span"],
                     s["within_2"] * 100, tol, s["b_tol"])

    d = pd.DataFrame(rows)
    t = sout.Topic("blacksmith", "flat_landscape_cause")
    d.to_csv(t.write("flat_landscape_cause", ".csv"), index=False)
    print("\n" + "=" * 76)
    print("  IS THE FLAT LANDSCAPE THE POCKET'S, OR OUR REACTIVE RECEPTOR'S?")
    print("=" * 76 + "\n")
    print(d.to_string(index=False))
    if {"reactive", "plain"} <= set(d.arm):
        p = d.pivot_table(index="ident", columns="arm",
                          values=["span", "within_2", "b_tol"])
        print("\n  reactive - plain, per molecule:")
        for c in ("span", "within_2", "b_tol"):
            if (c, "reactive") in p and (c, "plain") in p:
                dd = (p[(c, "reactive")] - p[(c, "plain")]).dropna()
                print(f"    {c:10s} mean {dd.mean():+.3f}   "
                      f"({'reactive flatter' if (c=='span' and dd.mean()<0) or (c!='span' and dd.mean()>0) else 'plain flatter'})")
    print()


if __name__ == "__main__":
    main()
