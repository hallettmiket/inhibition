#!/usr/bin/env python3
"""
Purpose: does a molecule's chemistry predict whether its pose cloud saturates?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: --molecules (a stratified list), docked fresh at --nrun
Output: 00_outputs/blacksmith/saturation_vs_chemistry/

@tt8804: "setup a genuine test with rigidity and other things considered to see
what we might need to move up the pipeline as a filter maybe? I dont want to
become completely biased to ridgid or small or large molecules though just find
some kind of reasonable cutoff."

WHAT IS MEASURED. For each molecule: dock, keep the PoseBusters-valid poses, and
fit the covering number against depth as a power law, `centres = a * n^b`. The
exponent b IS the saturation measure -- b = 0 means the cloud is fully explored
and deeper docking finds nothing new; b = 1 means every new pose is a new place.
D0090 measured b = 0.42 at 3.5 A on one molecule, so the question is whether that
varies with chemistry and whether any of it is predictable cheaply.

THE RESOLUTION IS NOT INVENTED HERE. 3.5 A is the pipeline's OWN tolerance:
`md.sweep_survivor_rmsd_nm = 0.35` nm is the ligand RMSD at which the sweep still
calls a pose "held". Two docked poses closer than that are within the tolerance
of the very next stage and would produce the same verdict, so resolving below it
is resolving distinctions the pipeline then discards (@tt8804). 2.0 A is carried
alongside as the Astex convention, and 5.0 A as a loose bound.

A COVERING SET IS A PARTITION OF THE VOLUME, which is what @tt8804 asked for:
"define the volume of poses and partion this volume into pose modes". The greedy
centres ARE the representatives -- unlike HDBSCAN this has an explicit length
scale, so its count is bounded for a bounded volume rather than growing with
density (D0090).

THE BIAS THIS EXPERIMENT COULD CAUSE, STATED UP FRONT. If b turns out to track
rigidity, the tempting move is to filter the library on rotatable bonds. That
would buy a saturating search by discarding flexible chemistry, which is a
selection on convenience rather than on merit -- and this library already cannot
separate the two, since rotatable bonds and heavy-atom count correlate at 0.51
across the in-scope set. Any cutoff proposed from this must be justified by what
it costs in coverage of real chemistry, not by how much tidier the clouds get.
"""

from __future__ import annotations

import argparse
import logging
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

log = logging.getLogger("sat-vs-chem")

RADII = (1.0, 2.0, 3.5, 5.0)
LADDER = (100, 250, 500, 1000, 2000)


def cover(coords: np.ndarray, r: float) -> int:
    """Greedy farthest-point covering number at radius r (see D0090)."""
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


def descriptors(smiles: str) -> dict:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors as rdd, Descriptors
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {}
    h = m.GetNumHeavyAtoms()
    return dict(rotb=rdd.CalcNumRotatableBonds(m), heavy=h,
                rot_frac=rdd.CalcNumRotatableBonds(m) / max(1, h),
                rings=rdd.CalcNumRings(m), fsp3=rdd.CalcFractionCSP3(m),
                tpsa=rdd.CalcTPSA(m), logp=Descriptors.MolLogP(m))


def dock_and_filter(cand, nrun: int, gpu: str, seed: int):
    """(heavy-atom coords of the PoseBusters-valid poses, pass rate)."""
    from rdkit import Chem
    from posebusters import PoseBusters
    rec_dir = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    work = Path(tempfile.mkdtemp(prefix="satchem_"))
    ligs = list(ns.prepare_ligand(cand, work / "lig.pdbqt"))
    if not ligs:
        raise RuntimeError(f"{cand.ident}: ligand preparation produced nothing")
    dlg = ns.dock(ligs[0], rec_dir, work / "c0", nrun, gpu, seed=seed)
    mol, _match = ns.rebuild_and_match(dlg, cand)
    heavy_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    xyz = np.array([mol.GetConformer(c).GetPositions()[heavy_idx]
                    for c in range(mol.GetNumConformers())])
    f = work / "poses.sdf"
    w = Chem.SDWriter(str(f))
    for c in range(mol.GetNumConformers()):
        w.write(mol, confId=c)
    w.close()
    df = PoseBusters(config="dock").bust([f], None, rp.receptor_prep())
    cols = [c for c in df.columns if df[c].dtype == bool]
    keep = df[cols].all(axis=1).to_numpy()
    return xyz[keep], float(keep.mean())


def fit_b(n: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Exponent and R^2 of y = a * n^b, fitted in log-log."""
    ok = (y > 0) & (n > 0)
    if ok.sum() < 3:
        return float("nan"), float("nan")
    ln, ly = np.log(n[ok]), np.log(y[ok])
    A = np.vstack([np.ones_like(ln), ln]).T
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ coef
    r2 = 1 - ((ly - pred) ** 2).sum() / max(1e-12, ((ly - ly.mean()) ** 2).sum())
    return float(coef[1]), float(r2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--molecules", required=True,
                    help="file of candidate ids, one per line")
    ap.add_argument("--nrun", type=int, default=2000)
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    want = [l.strip() for l in Path(a.molecules).read_text().splitlines() if l.strip()]
    cands = {c.ident: c for c in nr.load_candidates()}
    missing = [w for w in want if w not in cands]
    if missing:
        raise SystemExit(f"not in the candidate table: {missing}")

    tol = float(tc.get("md.sweep_survivor_rmsd_nm")) * 10.0
    log.info("pipeline tolerance = %.1f A (md.sweep_survivor_rmsd_nm)", tol)

    rows = []
    for i, ident in enumerate(want, 1):
        cand = cands[ident]
        log.info("[%d/%d] %s", i, len(want), ident)
        try:
            xyz, rate = dock_and_filter(cand, a.nrun, a.gpu, a.seed)
        except Exception as exc:                                   # noqa: BLE001
            log.warning("  %s failed: %s", ident, exc)
            continue
        log.info("  %d/%d PoseBusters-valid (%.1f%%)",
                 len(xyz), a.nrun, rate * 100)
        rng = np.random.default_rng(a.seed)
        ladder = [k for k in LADDER if k <= len(xyz)] + [len(xyz)]
        ladder = sorted(set(ladder))
        curves = {r: [] for r in RADII}
        for k in ladder:
            idx = rng.choice(len(xyz), size=k, replace=False)
            for r in RADII:
                curves[r].append(cover(xyz[idx], r))
        row = dict(ident=ident, pb_rate=round(rate, 4), n_valid=len(xyz),
                   **descriptors(cand.smiles))
        n = np.array(ladder, float)
        for r in RADII:
            b, r2 = fit_b(n, np.array(curves[r], float))
            row[f"b_{r}a"] = round(b, 4)
            row[f"r2_{r}a"] = round(r2, 4)
            row[f"centres500_{r}a"] = (curves[r][ladder.index(500)]
                                       if 500 in ladder else None)
        rows.append(row)
        log.info("  b @ %.1f A = %.3f   centres at 500 poses = %s",
                 tol, row.get(f"b_{tol}a"), row.get(f"centres500_{tol}a"))

    d = pd.DataFrame(rows)
    t = sout.Topic("blacksmith", "saturation_vs_chemistry")
    d.to_csv(t.write("saturation_vs_chemistry", ".csv"), index=False)

    print("\n" + "=" * 76)
    print("  DOES CHEMISTRY PREDICT WHETHER THE CLOUD SATURATES?")
    print("=" * 76)
    cols = ["ident", "rotb", "heavy", "rings", "fsp3", "pb_rate",
            f"b_{tol}a", f"centres500_{tol}a"]
    print("\n" + d[[c for c in cols if c in d.columns]].to_string(index=False))
    from scipy.stats import spearmanr
    print(f"\n  saturation exponent b at {tol:.1f} A against cheap descriptors "
          f"(n = {len(d)}):")
    print("  (a LOWER b means the cloud saturates sooner)")
    for c in ("rotb", "rot_frac", "heavy", "rings", "fsp3", "tpsa", "logp"):
        if c not in d.columns or d[c].notna().sum() < 4:
            continue
        rho, p = spearmanr(d[c], d[f"b_{tol}a"], nan_policy="omit")
        print(f"    {c:10s} rho = {rho:+.3f}   p = {p:.3f}"
              f"{'   *' if p < 0.05 else ''}")
    print(f"\n  b across the set: min {d[f'b_{tol}a'].min():.3f}  "
          f"median {d[f'b_{tol}a'].median():.3f}  max {d[f'b_{tol}a'].max():.3f}")
    print("\n  NOTE: n is small and this library spans heavy atoms 19-29 with "
          "rotb\n  correlated to size at 0.51 -- it cannot separate flexibility "
          "from bulk.\n")


if __name__ == "__main__":
    main()
