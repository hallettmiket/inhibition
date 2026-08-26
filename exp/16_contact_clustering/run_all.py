#!/usr/bin/env python3
"""
Purpose: does contact-space grouping produce modes whose poses actually resemble each other?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: --molecules (default: 20 spanning the library), their persisted clouds
Output: 00_outputs/blacksmith/contact_clustering/

THE TEST THAT KILLED THE LAST TWO PROPOSALS, APPLIED FIRST THIS TIME. A 3 A
spatial partition put 1,758 of 6,000 poses into one cell whose members differed
by 9.13 A (D0091); the shipped DBSCAN rule put 137 poses in a 9.3 A bag (D0088).
Both were reported before anyone asked whether a group's members resemble each
other. So the headline here is WITHIN-GROUP CARTESIAN RMSD, and it is computed
for every group of every molecule.

Complete linkage guarantees the within-group CONTACT distance is under the cut --
that is structural and is asserted, not hoped. It does NOT guarantee anything
about Cartesian RMSD, and whether the two agree is exactly the open question.

RESIDUES COME FROM exp/14, NOT FROM A LITERATURE LIST OR A CUTOFF. The greedy
non-redundant ranking, whose order is identical across five independent dockings
(Spearman 1.000) and whose top 15 span 90% of the contact matrix's variance.
WATERS ARE DROPPED: `A:40:HOH` ranked 13th, and a landmark that is modelled
inconsistently between structures is not a landmark.

THE TOLERANCE IS THE MOLECULE'S OWN. Predicted per-atom RMSF from a conformer
ensemble (exp/15: rho = 0.657 over 147 molecules, 100% positive), calibrated by
the measured 2.21x overestimate. A floppy molecule gets a looser cut than a rigid
one, and nobody picks a number.
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
sys.path.insert(0, str(REPO / "exp" / "15_rmsf_predictor"))

from shared import outputs as sout                  # noqa: E402
from shared import pose_contacts as pc              # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("contact-clustering")


def key_residues(n: int) -> list:
    """Top `n` from exp/14's greedy non-redundant pick, waters excluded."""
    fs = sorted(glob.glob(str(rp.BLACKSMITH / "residue_selection_*/greedy_pick_*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit("run exp/14_residue_selection first")
    d = pd.read_csv(fs[-1]).sort_values("order")
    keep = [r for r in d.residue if not r.endswith(":HOH")]
    if len(keep) < n:
        raise SystemExit(f"only {len(keep)} non-water residues in {fs[-1]}")
    return keep[:n]


def receptor_coords(names: list) -> list:
    """Heavy-atom coordinates for each named residue, in the order given."""
    want = {}
    for nm in names:
        c, i, rn = nm.split(":")
        want[(c, i, rn)] = []
    for ln in rp.receptor_prep().read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        k = (ln[21:22].strip() or "A", ln[22:26].strip(), ln[17:20].strip())
        if k in want:
            want[k].append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    out = []
    for nm in names:
        c, i, rn = nm.split(":")
        v = want[(c, i, rn)]
        if not v:
            raise SystemExit(f"{nm} has no heavy atoms in the receptor")
        out.append(np.array(v))
    return out


def load_cloud(ident: str):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    f = rp.allposes_dir() / f"{ident}.sdf"
    if not f.is_file():
        return None, None
    ms = [m for m in Chem.SDMolSupplier(str(f), removeHs=False, sanitize=True)
          if m is not None]
    if len(ms) < 30:
        return None, None
    hv = [a.GetIdx() for a in ms[0].GetAtoms() if a.GetAtomicNum() > 1]
    xyz = np.array([m.GetConformer().GetPositions()[hv] for m in ms])
    return xyz, (ms[0], hv)


def cartesian_rmsd(xyz: np.ndarray, idx: np.ndarray) -> tuple:
    """(median, max) unaligned heavy-atom RMSD inside one group."""
    if len(idx) < 2:
        return 0.0, 0.0
    s = xyz[idx] if len(idx) <= 200 else xyz[np.random.default_rng(0).choice(idx, 200, replace=False)]
    D = np.array([np.sqrt(((s - s[i]) ** 2).sum(-1).mean(-1)) for i in range(len(s))])
    iu = np.triu_indices(len(s), 1)
    return float(np.median(D[iu])), float(D[iu].max())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--n-molecules", type=int, default=20)
    ap.add_argument("--residues", type=int, default=15)
    ap.add_argument("--conformers", type=int, default=50)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    from run_all import predict_rmsf                                # exp/15
    names = key_residues(a.residues)
    res = receptor_coords(names)
    log.info("landmarks (%d): %s", len(names), ", ".join(names))

    fs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_v2/rank_v2_T4_{rp.topic()}_conditional_eb_*.csv")),
                key=os.path.getmtime)
    rk = pd.read_csv(fs[-1])
    mols = [m for m in rk.parent_ident.dropna().unique()
            if (rp.allposes_dir() / f"{m}.sdf").is_file()]
    rng = np.random.default_rng(a.seed)
    mols = list(rng.choice(mols, size=min(a.n_molecules, len(mols)), replace=False))

    rows, gr = [], []
    for n, ident in enumerate(mols, 1):
        xyz, meta = load_cloud(ident)
        if xyz is None:
            continue
        tmpl, hv = meta
        try:
            rmsf = predict_rmsf(tmpl, hv, a.conformers, a.seed)
        except Exception as exc:                                    # noqa: BLE001
            log.warning("  %s: no RMSF (%s)", ident, str(exc)[:50])
            continue
        w = pc.atom_weights(rmsf)
        tol = float(np.median(rmsf) / pc.RMSF_CALIBRATION)
        T = pc.contact_tensor(xyz, res)
        D = pc.pose_distances(T, w)
        lab = pc.group(D, tol)
        worst_contact = pc.within_group_max(D, lab)
        assert worst_contact <= tol + 1e-6, (
            f"complete linkage broke its own guarantee: {worst_contact:.3f} > {tol:.3f}")
        sizes = np.bincount(lab)
        med, mx = [], []
        for k in range(len(sizes)):
            idx = np.flatnonzero(lab == k)
            m1, m2 = cartesian_rmsd(xyz, idx)
            med.append(m1); mx.append(m2)
            gr.append(dict(ident=ident, group=k, size=len(idx),
                           cart_median=m1, cart_max=m2))
        rows.append(dict(
            ident=ident, poses=len(xyz), atoms=xyz.shape[1], tol_a=tol,
            groups=len(sizes), largest=int(sizes.max()),
            singletons=int((sizes == 1).sum()),
            worst_contact_a=worst_contact,
            cart_median_of_groups=float(np.median([m for m, s in zip(med, sizes) if s > 1] or [0])),
            cart_max_any_group=float(max(mx) if mx else 0.0)))
        log.info("  %2d/%d %s: %d poses -> %d groups, largest %d, tol %.2f A, "
                 "worst within-group Cartesian %.2f A",
                 n, len(mols), ident, len(xyz), len(sizes), sizes.max(), tol, max(mx))

    d = pd.DataFrame(rows); g = pd.DataFrame(gr)
    t = sout.Topic("blacksmith", "contact_clustering")
    d.to_csv(t.write("per_molecule", ".csv"), index=False)
    g.to_csv(t.write("per_group", ".csv"), index=False)

    print("\n" + "=" * 78)
    print("  CONTACT-SPACE GROUPING — do a group's poses resemble each other?")
    print("=" * 78)
    print(f"\n  molecules: {len(d)}   poses: {int(d.poses.sum()):,}   "
          f"landmarks: {len(names)}")
    print(f"\n  GROUPS PER MOLECULE: median {d.groups.median():.0f}  "
          f"range {d.groups.min()}-{d.groups.max()}")
    print(f"  largest group:        median {d.largest.median():.0f}  "
          f"max {d.largest.max()}  ({d.largest.max()/d.poses.max()*100:.0f}% of a cloud)")
    print(f"  singleton groups:     median {d.singletons.median():.0f} "
          f"({(g['size']==1).mean()*100:.0f}% of all groups)")
    print(f"  tolerance used:       median {d.tol_a.median():.2f} A  "
          f"range {d.tol_a.min():.2f}-{d.tol_a.max():.2f}")
    print(f"\n  THE GUARANTEE (contact distance within a group <= tol):")
    print(f"    largest violation across every group of every molecule: "
          f"{(d.worst_contact_a - d.tol_a).max():+.2e} A   -- holds")
    print(f"\n  THE OPEN QUESTION (Cartesian RMSD inside a group):")
    print(f"    median within-group RMSD, over multi-pose groups: "
          f"{g[g['size']>1].cart_median.median():.2f} A")
    print(f"    90th percentile:                                  "
          f"{g[g['size']>1].cart_median.quantile(.9):.2f} A")
    print(f"    WORST any group anywhere:                         "
          f"{g.cart_max.max():.2f} A")
    print(f"\n  for comparison, the two rules this replaces:")
    print(f"    shipped DBSCAN  : largest 137 poses spanning 9.3 A  (D0088)")
    print(f"    3 A partition   : largest 1,758 poses spanning 9.13 A (D0091)")
    print()


if __name__ == "__main__":
    main()
