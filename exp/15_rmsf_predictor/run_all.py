#!/usr/bin/env python3
"""
Purpose: can a cheap conformer ensemble predict the per-atom flexibility MD measures?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: the 147 swept modes' trajectories under the sweep work root
Output: 00_outputs/blacksmith/rmsf_predictor/

WHY A PREDICTOR AND NOT THE MEASUREMENT (@tt8804). The splitting tolerance is
needed BEFORE the sweep, to decide what gets swept. Measured RMSF is only
available afterwards, and only for molecules that were already selected -- the
set that no longer needs splitting. So the tolerance must be predicted, and the
147 measured trajectories are the VALIDATION SET, not the input. That is the
only arrangement that is not circular.

THE PREDICTOR. Embed ~50 conformers with ETKDG, MMFF-minimise, align on heavy
atoms, take the per-atom spread. It needs no protein, no trajectory and no
docking -- seconds per molecule -- and it measures the thing the tolerance is
for: how far an atom can move without the molecule becoming a different molecule.

BOND ORDERS COME FROM THE POSE, NOT FROM THE COORDINATES. A first attempt built
the molecule from the MD's heavy-atom coordinates alone and let RDKit perceive
bonds. It failed loudly enough to notice -- "Final molecular charge (0) does not
match input (3)" -- and produced [S+2], [C+] and [N-] where the real molecule has
none, so the conformers came from a molecule that is not ours. Correlation was
still +0.51, which is exactly the trap: a broken input that yields a plausible
number. The template route (`AssignBondOrdersFromTemplate`, as
`crystal_pose_audit` already uses) takes chemistry from the pose SDF, which RDKit
reads correctly.

ATOM ORDER IS MATCHED BY INTERNAL GEOMETRY, NOT BY NAME OR BY POSITION. The MD
topology's ligand atoms come from an antechamber mol2 whose atom TYPES (`ca`,
`nb`, `os`) parse as ELEMENTS -- RDKit refuses the file outright and openbabel
silently returns calcium, niobium and osmium, 14 heavy atoms where there are 25.
So no name or type is trusted.

Raw coordinates do not work either: GROMACS centres the system in its box, so the
ligand is translated by tens of angstroms (measured: 16.8, 27.4, 37.1 A on the
first molecule) and a nearest-neighbour match on absolute positions is not even a
bijection. What IS invariant is the molecule's own shape, so each atom is matched
on its SORTED VECTOR OF DISTANCES TO EVERY OTHER LIGAND ATOM -- unchanged by any
translation or rotation. Measured on the same molecule: a perfect bijection at
0.030 A signature error and 0.005 A centred RMSD.

THE ASSIGNMENT IS HUNGARIAN, NOT NEAREST-NEIGHBOUR. A first version took each
atom's best partner independently, which is not a bijection whenever two atoms
share a signature -- and 28 of 147 molecules did. None is symmetric: they are
combinatorial products. What they have is LOCAL topological equivalence -- a
para-substituted phenyl's two ortho carbons, a sulfone's two oxygens, a
gem-dimethyl -- verified: all 28 carry topologically equivalent atoms, median 6,
and the colliding elements are ring carbons and paired oxygens. Equivalent atoms
are interchangeable, so either assignment is correct; an optimal global
assignment simply picks one and is consistent about it.

NO HEAVY-ONLY MOLECULE IS EVER BUILT. Deleting hydrogens to get a heavy-atom
molecule makes an aromatic N-H unkekulizable -- the failure `redock_04_rmsd`
already documents, and it cost 9 more molecules here. Only an INDEX MAPPING is
carried; conformers are embedded on the intact molecule and read at the mapped
positions.

The match is asserted to be a bijection and to reproduce the conformation; it
raises rather than proceeding on a partial correspondence.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("rmsf-predictor")

GMX = Path("/data/lab_vm/envs/dwi_gromacs_cuda/bin.AVX2_256/gmx")
#: how close a topology atom must be to its SDF partner for the match to be
#: accepted. The system is built FROM the pose, so this should be ~0.
MATCH_TOL_A = 0.35


def lig_indices(rep: Path) -> list[int]:
    out = []
    for l in (rep / "lig.ndx").read_text().splitlines()[1:]:
        out += [int(x) for x in l.split() if x.isdigit()]
    return out


def gro_ligand(gro: Path, nums: list[int]):
    """(coords in A, atom names) for the given 1-based atom numbers, in order."""
    lines = gro.read_text().splitlines()
    n = int(lines[1])
    want = set(nums)
    sel = {}
    for l in lines[2:2 + n]:
        an = int(l[15:20])
        if an in want:
            sel[an] = (float(l[20:28]) * 10, float(l[28:36]) * 10,
                       float(l[36:44]) * 10, l[10:15].strip())
    if len(sel) != len(nums):
        raise ValueError(f"{len(sel)} of {len(nums)} ligand atoms found in {gro.name}")
    return (np.array([sel[i][:3] for i in nums]),
            [sel[i][3] for i in nums])


def measured_rmsf(rep: Path, work: Path) -> np.ndarray:
    """Per-atom RMSF in nm from the production trajectory, topology order."""
    out = work / "rmsf.xvg"
    r = subprocess.run(
        [str(GMX), "rmsf", "-f", str(rep / "whole.xtc"), "-s", str(rep / "prod.tpr"),
         "-n", str(rep / "lig.ndx"), "-o", str(out), "-res", "no"],
        input="0\n", capture_output=True, text=True, timeout=900)
    if not out.is_file():
        raise RuntimeError("gmx rmsf produced nothing: "
                           + "\n".join((r.stderr or "").strip().splitlines()[-3:]))
    return np.array([float(l.split()[1]) for l in out.read_text().splitlines()
                     if l and not l.startswith(("@", "#"))])


def pose_mol(ident: str, pose_rank: int):
    """The simulated pose, with correct chemistry, from the representatives SDF."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    f = rp.poses_dir() / f"{ident}.sdf"
    if not f.is_file():
        return None
    for m in Chem.SDMolSupplier(str(f), removeHs=False, sanitize=True):
        if m is None or not m.HasProp("pose_rank"):
            continue
        if int(m.GetProp("pose_rank")) == int(pose_rank):
            return m
    return None


def topology_order_map(md_xyz: np.ndarray, tmpl):
    """(indices into `tmpl` in MD topology order, centred RMSD of the match).

    An index map, not a rebuilt molecule -- see the module docstring on why a
    heavy-only copy cannot be made safely.
    """
    heavy = [a.GetIdx() for a in tmpl.GetAtoms() if a.GetAtomicNum() > 1]
    txyz = np.array(tmpl.GetConformer().GetPositions())[heavy]
    if len(heavy) != len(md_xyz):
        raise ValueError(f"{len(heavy)} SDF heavy atoms vs {len(md_xyz)} in the topology")

    def signature(x):
        """Each atom's sorted distances to every other atom -- pose-invariant."""
        d = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(-1))
        return np.sort(d, axis=1)

    from scipy.optimize import linear_sum_assignment
    C = np.sqrt(((signature(md_xyz)[:, None, :]
                  - signature(txyz)[None, :, :]) ** 2).sum(-1))
    # OPTIMAL GLOBAL ASSIGNMENT. argmin per row is not a bijection whenever two
    # atoms share a signature, which local topological equivalence guarantees.
    rows, order = linear_sum_assignment(C)
    if not np.array_equal(rows, np.arange(len(md_xyz))) or \
            len(set(order.tolist())) != len(order):
        raise ValueError("assignment is not a bijection")
    mc = md_xyz - md_xyz.mean(0)
    tc = txyz - txyz.mean(0)
    worst = float(np.sqrt(((mc - tc[order]) ** 2).sum(1).mean()))
    if worst > MATCH_TOL_A:
        raise ValueError(f"matched atoms give {worst:.2f} A centred RMSD; the MD "
                         f"system should hold this exact pose")
    return [heavy[i] for i in order], worst


def predict_rmsf(tmpl, order: list[int], n_conf: int, seed: int) -> np.ndarray:
    """Per-atom spread over an aligned conformer ensemble, in MD topology order.

    Embedded on the INTACT molecule and read at `order`; no heavy-only copy is
    made, so an aromatic N-H never loses the hydrogen that makes it kekulizable.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    mh = Chem.Mol(tmpl)
    cids = AllChem.EmbedMultipleConfs(mh, numConfs=n_conf, randomSeed=seed)
    if len(cids) == 0:
        cids = AllChem.EmbedMultipleConfs(mh, numConfs=n_conf, randomSeed=seed,
                                          useRandomCoords=True)
    if len(cids) < 5:
        raise ValueError(f"only {len(cids)} conformers embedded")
    AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=300)
    amap = [(i, i) for i in order]
    for c in list(cids)[1:]:
        rdMolAlign.AlignMol(mh, mh, prbCid=c, refCid=cids[0], atomMap=amap)
    X = np.array([mh.GetConformer(c).GetPositions()[order] for c in cids])
    return np.sqrt(((X - X.mean(0)) ** 2).sum(-1).mean(0))


def sweep_rows() -> pd.DataFrame:
    fs = sorted(glob.glob(str(rp.sweep_dir() / "attack_sweep_*.csv")),
                key=os.path.getmtime)
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d.drop_duplicates("ident", keep="last")
    return d[d.status.astype(str).str.startswith("ok")]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--conformers", type=int, default=50)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    sw = sweep_rows()
    log.info("swept modes with an ok result: %d", len(sw))
    root = rp.sweep_work()
    rows, skipped = [], {}
    from scipy.stats import spearmanr, pearsonr

    items = list(sw.itertuples())
    if a.limit:
        items = items[:a.limit]
    for n, r in enumerate(items, 1):
        ident = str(r.parent_ident)
        # THE DIRECTORY IS CHOSEN BY RANK, NOT BY SORT ORDER. `attack_sweep`
        # writes each run to `rank{pose_rank}_{ps}ps/`, so a molecule swept at
        # several ranks has several directories. Taking `sorted(...)[0]` read
        # `rank1_` for a row recording rank 5 -- and the ligand then failed the
        # conformation check at 4-6 A, which looked exactly like a pose-provenance
        # defect in the pipeline. It was selection by position, here.
        hits = sorted(root.glob(f"rank{int(r.pose_rank)}_*ps/{ident}/md"))
        if not hits:
            skipped["no md dir for this rank"] = \
                skipped.get("no md dir for this rank", 0) + 1
            continue
        md = hits[0]
        rep = md / "rep1"
        if not (rep / "whole.xtc").is_file():
            skipped["no trajectory"] = skipped.get("no trajectory", 0) + 1
            continue
        try:
            nums = lig_indices(rep)
            md_xyz, _names = gro_ligand(md / "sys.gro", nums)
            tmpl = pose_mol(ident, r.pose_rank)
            if tmpl is None:
                skipped["no pose"] = skipped.get("no pose", 0) + 1
                continue
            order, worst = topology_order_map(md_xyz, tmpl)
            pred = predict_rmsf(tmpl, order, a.conformers, a.seed)
            with tempfile.TemporaryDirectory() as td:
                meas = measured_rmsf(rep, Path(td)) * 10.0     # nm -> A
            if len(meas) != len(pred):
                skipped["length mismatch"] = skipped.get("length mismatch", 0) + 1
                continue
            rho = spearmanr(pred, meas)[0]
            rows.append(dict(ident=ident, n_atoms=len(pred), match_a=round(worst, 3),
                             spearman=rho, pearson=pearsonr(pred, meas)[0],
                             pred_med=float(np.median(pred)),
                             meas_med=float(np.median(meas)),
                             pred_min=float(pred.min()), pred_max=float(pred.max()),
                             meas_min=float(meas.min()), meas_max=float(meas.max())))
        except Exception as exc:                               # noqa: BLE001
            k = str(exc)[:60]
            skipped[k] = skipped.get(k, 0) + 1
            continue
        if n % 20 == 0:
            log.info("  %d/%d, %d usable", n, len(items), len(rows))

    d = pd.DataFrame(rows)
    if d.empty:
        print("\nnothing usable. reasons:", skipped)
        raise SystemExit(1)
    t = sout.Topic("blacksmith", "rmsf_predictor")
    d.to_csv(t.write("rmsf_predictor", ".csv"), index=False)

    print("\n" + "=" * 76)
    print("  CAN A CONFORMER ENSEMBLE PREDICT THE FLEXIBILITY MD MEASURES?")
    print("=" * 76)
    print(f"\n  molecules validated: {len(d)}   atoms: {int(d.n_atoms.sum()):,}")
    if skipped:
        print("  skipped:")
        for k, v in sorted(skipped.items(), key=lambda x: -x[1])[:6]:
            print(f"    {v:>4}  {k}")
    print(f"\n  PER-MOLECULE rank correlation (predicted vs measured, across atoms):")
    print(f"    median {d.spearman.median():+.3f}   "
          f"IQR {d.spearman.quantile(.25):+.3f} to {d.spearman.quantile(.75):+.3f}")
    print(f"    fraction positive: {(d.spearman > 0).mean()*100:.0f}%   "
          f"fraction above +0.5: {(d.spearman > 0.5).mean()*100:.0f}%")
    print(f"\n  MAGNITUDE (does it get the scale right, or only the order?):")
    print(f"    predicted median per-atom RMSF {d.pred_med.median():.2f} A")
    print(f"    measured  median per-atom RMSF {d.meas_med.median():.2f} A")
    print(f"    ratio pred/meas: {d.pred_med.median()/d.meas_med.median():.2f}")
    print(f"\n  coordinate match quality: worst {d.match_a.max():.3f} A "
          f"(tolerance {MATCH_TOL_A})")
    print()


if __name__ == "__main__":
    main()
