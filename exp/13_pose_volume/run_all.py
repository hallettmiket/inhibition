#!/usr/bin/env python3
"""
Purpose: does the VOLUME the poses occupy grow with sampling, or is it bounded?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: a persisted deep cloud, or --candidate to dock one
Output: 00_outputs/blacksmith/pose_volume_<candidate>/

@tt8804: "my intuition is that the use of pose busters and the nature of docking
sets a bounded pose space that is contained within the receptor volume and we can
partition this space by 3A and assign poses per partition. We will need to test
this first, whether pose space volume grows."

WHY THIS IS A DIFFERENT QUESTION FROM D0090's, AND WHY THAT MATTERS. D0090
measured the COVERING NUMBER -- how many balls of radius r are needed to cover
the poses -- and found it grows as n^0.42 with no asymptote, so pose space never
saturates. But a covering number can grow while the region being covered stays
fixed: finer and finer resolution WITHIN a bounded volume, rather than new
territory. Those are opposite conclusions and the covering number cannot tell
them apart.

The volume can. If the occupied region stops expanding, docking is exploring a
bounded space and the covering growth is granularity; and a bounded space can be
partitioned once, at a stated resolution, with a count that does not move.

TWO OCCUPANCY MODELS, BECAUSE THE NAIVE ONE OVERSTATES GROWTH. Counting a voxel
as occupied when an atom CENTRE falls in it is boundary-sensitive: an added pose
tips cells that were nearly occupied, so the count creeps for a reason that is
about the lattice rather than about the molecule reaching anywhere new. Counting
an atom as its van der Waals SPHERE is both more physical and less
boundary-sensitive, and the gap between the two curves is the size of the
artefact. Both are reported.

POSEBUSTERS FIRST, because the premise is about the space PoseBusters leaves. It
removes ~10% of poses (D0089), disproportionately ones clashing into protein --
which sit at the EDGE of the occupied region, so filtering should shrink the
volume and bound it further. Running unfiltered would overstate the extent.

NOT PARTITIONED ON THE REACTIVE ATOM. Its position IS the distance term the
criterion scores on, so grouping by it and then scoring would be D0088's
circularity in a new coordinate system. All-atom occupancy and the centroid are
reported instead; the reactive atom is carried ONLY as a diagnostic and is
labelled as such.
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

log = logging.getLogger("pose-volume")

LADDER = (100, 250, 500, 1000, 2000, 3500)
VOXELS = (1.0, 1.5, 3.0)
#: van der Waals radii, A. Only the elements this library contains.
VDW = {6: 1.70, 7: 1.55, 8: 1.52, 16: 1.80, 9: 1.47, 17: 1.75, 35: 1.85, 53: 1.98}


def point_cells(pts: np.ndarray, s: float) -> set:
    """Voxels containing an atom CENTRE. Boundary-sensitive by construction."""
    return set(map(tuple, np.floor(np.asarray(pts).reshape(-1, 3) / s).astype(int)))


def sphere_cells(pts: np.ndarray, radii: np.ndarray, s: float) -> set:
    """Voxels overlapped by an atom's van der Waals SPHERE.

    The physical occupancy: a carbon is 1.7 A across, not a point, so a voxel it
    touches is occupied whether or not the centre landed inside. Enumerated over
    the sphere's bounding box and tested against the voxel centre, which is the
    standard approximation and is exact to within half a voxel diagonal.
    """
    out = set()
    pts = np.asarray(pts).reshape(-1, 3)
    radii = np.asarray(radii).reshape(-1)
    for p, r in zip(pts, radii):
        lo = np.floor((p - r) / s).astype(int)
        hi = np.floor((p + r) / s).astype(int)
        for i in range(lo[0], hi[0] + 1):
            for j in range(lo[1], hi[1] + 1):
                for k in range(lo[2], hi[2] + 1):
                    c = (np.array([i, j, k]) + 0.5) * s
                    if ((c - p) ** 2).sum() <= r * r:
                        out.add((i, j, k))
    return out


def cover(coords: np.ndarray, r: float) -> int:
    """Greedy farthest-point covering number, from the medoid (D0090, exp/8)."""
    n = len(coords)
    tot = np.zeros(n)
    for i in range(n):
        tot += np.sqrt(((coords - coords[i]) ** 2).sum(axis=2).mean(axis=1))
    nxt = int(np.argmin(tot))
    dmin = np.full(n, np.inf)
    chosen = 0
    while True:
        d = np.sqrt(((coords - coords[nxt]) ** 2).sum(axis=2).mean(axis=1))
        dmin = np.minimum(dmin, d)
        chosen += 1
        far = int(np.argmax(dmin))
        if dmin[far] <= r or chosen >= n:
            return chosen
        nxt = far


def load_cloud(path: Path):
    """(coords (n,a,3), vdW radii per heavy atom, the mols)."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    ms = [m for m in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
          if m is not None]
    if not ms:
        raise SystemExit(f"no poses in {path}")
    heavy = [i for i, a in enumerate(ms[0].GetAtoms()) if a.GetAtomicNum() > 1]
    radii = np.array([VDW.get(ms[0].GetAtomWithIdx(i).GetAtomicNum(), 1.7)
                      for i in heavy])
    xyz = np.array([m.GetConformer().GetPositions()[heavy] for m in ms])
    return xyz, radii, ms


def pb_mask(mols) -> np.ndarray:
    from rdkit import Chem
    from posebusters import PoseBusters
    td = Path(tempfile.mkdtemp(prefix="pv_pb_"))
    f = td / "p.sdf"
    w = Chem.SDWriter(str(f))
    for m in mols:
        w.write(m)
    w.close()
    df = PoseBusters(config="dock").bust([f], None, rp.receptor_prep())
    cols = [c for c in df.columns if df[c].dtype == bool]
    return df[cols].all(axis=1).to_numpy()


def fit_b(n, y) -> float:
    n, y = np.asarray(n, float), np.asarray(y, float)
    ok = (n > 0) & (y > 0)
    if ok.sum() < 3:
        return float("nan")
    A = np.vstack([np.ones(ok.sum()), np.log(n[ok])]).T
    coef, *_ = np.linalg.lstsq(A, np.log(y[ok]), rcond=None)
    return float(coef[1])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--cloud", default="", help="path to a persisted cloud SDF")
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--no-posebusters", action="store_true",
                    help="skip the filter; the premise is about the filtered "
                         "space, so this is a diagnostic only")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    path = (Path(a.cloud) if a.cloud else
            rp.BLACKSMITH / f"deep_cloud_{a.candidate}" / "cloud_1.sdf")
    xyz, radii, mols = load_cloud(path)
    log.info("cloud: %d poses x %d heavy atoms from %s", len(xyz), xyz.shape[1],
             path.name)

    if not a.no_posebusters:
        keep = pb_mask(mols)
        log.info("PoseBusters: %d of %d valid (%.1f%%)",
                 int(keep.sum()), len(keep), 100 * keep.mean())
        xyz = xyz[keep]

    tol = float(tc.get("md.sweep_survivor_rmsd_nm")) * 10.0
    rng = np.random.default_rng(a.seed)
    ladder = sorted(set([k for k in LADDER if k <= len(xyz)] + [len(xyz)]))
    rows = []
    for k in ladder:
        idx = rng.choice(len(xyz), size=k, replace=False)
        sub = xyz[idx]
        row = {"poses": k}
        for s in VOXELS:
            row[f"pt_cells_{s}"] = len(point_cells(sub, s))
            row[f"pt_vol_{s}"] = len(point_cells(sub, s)) * s ** 3
        # the physical volume: one voxel size is enough, and the fine one is
        # where the boundary artefact would be worst
        sph = sphere_cells(sub.reshape(-1, 3),
                           np.tile(radii, len(sub)), 1.0)
        row["sphere_vol_1.0"] = len(sph) * 1.0
        row["centroid_cells_3.0"] = len(point_cells(sub.mean(axis=1), 3.0))
        row[f"cover_{tol}"] = cover(sub, tol)
        rows.append(row)
        log.info("  %5d poses: volume %.0f A^3 (spheres), %d cells @3A, "
                 "%d cover @%.1fA", k, row["sphere_vol_1.0"],
                 row["pt_cells_3.0"], row[f"cover_{tol}"], tol)

    d = pd.DataFrame(rows)
    t = sout.Topic("blacksmith", f"pose_volume_{a.candidate}")
    d.to_csv(t.write("pose_volume", ".csv"), index=False)

    n = d.poses.values
    print("\n" + "=" * 76)
    print(f"  DOES THE OCCUPIED VOLUME GROW? — {a.candidate}"
          f"{'' if not a.no_posebusters else '  (UNFILTERED)'}")
    print("=" * 76 + "\n")
    print(d.to_string(index=False))
    print("\n  power-law exponent b in a*n^b   (1 = linear, 0 = bounded)")
    for c in [c for c in d.columns if c != "poses"]:
        print(f"    {c:22s} b = {fit_b(n, d[c].values):+.3f}")
    print("\n  the boundary artefact is the gap between point and sphere "
          "occupancy:")
    bp, bs = fit_b(n, d["pt_vol_1.0"].values), fit_b(n, d["sphere_vol_1.0"].values)
    print(f"    point-occupancy b {bp:+.3f}   sphere-occupancy b {bs:+.3f}   "
          f"difference {bp - bs:+.3f}")
    print()


if __name__ == "__main__":
    main()
