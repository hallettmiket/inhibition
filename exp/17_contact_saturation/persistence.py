#!/usr/bin/env python3
"""
Purpose: do the groups found at shallow depth survive at deep depth, or does depth rewrite them?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: the 6,000-pose RAW deep cloud
Output: 00_outputs/blacksmith/contact_saturation/

THE COUNT DOES NOT SATURATE (tolerance_sweep.py: b = +0.69 at 0.73 A, +0.57 at
1.0 A, no finite plateau at any tolerance that keeps groups tight). That is only
fatal if it means the ANSWER changes with depth. It does not have to: a fixed
absolute tolerance carves fixed regions of contact space, so a deeper cloud can
only add regions -- it cannot move the ones already found. If every shallow group
still exists at 6,000 poses and the growth is all new, sparse, low-population
regions, then non-saturation is a statement about the tail and the shortlist is
safe. If shallow groups DISSOLVE, the count is the least of the problems.

MATCHED BY MEDOID IN CONTACT SPACE, at the same tolerance used to build them --
the same rule that decides membership decides identity, so a "match" cannot be a
looser test than a "group".
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import pose_contacts as pc              # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("persistence")
_spec = importlib.util.spec_from_file_location(
    "x17", Path(__file__).with_name("run_all.py"))
_M = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_M)


def flat(T, w):
    """The contact vectors themselves, in the space the metric is Euclidean in."""
    W = np.sqrt(w / w.sum())[None, :, None]
    return (T * W).reshape(len(T), -1) / np.sqrt(T.shape[2])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--residues", type=int, default=15)
    ap.add_argument("--shallow", type=int, default=500)
    ap.add_argument("--draws", type=int, default=5)
    ap.add_argument("--tol", type=float, default=None, help="default: the RMSF tolerance")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    res = _M.receptor_coords(_M.key_residues(a.residues))
    xyz, meta = _M.load_sdf(rp.BLACKSMITH / f"deep_cloud_{a.candidate}" / "cloud_1.sdf")
    rmsf = _M.predict_rmsf(meta[0], meta[1], 50, a.seed)
    w = pc.atom_weights(rmsf)
    tol = a.tol if a.tol else float(np.median(rmsf) / pc.RMSF_CALIBRATION)

    T_all = pc.contact_tensor(xyz, res)
    V = flat(T_all, w)                                  # one vector per pose
    D_all = pc.pose_distances(T_all, w)
    lab_all = pc.group(D_all, tol)
    n_deep = lab_all.max() + 1
    deep_med = np.array([V[np.flatnonzero(lab_all == k)].mean(0)
                         for k in range(n_deep)])
    deep_sz = np.bincount(lab_all)
    log.info("deep: %d poses -> %d groups at tol %.2f A", len(xyz), n_deep, tol)

    rng = np.random.default_rng(a.seed)
    rows = []
    for d in range(a.draws):
        idx = rng.choice(len(xyz), size=a.shallow, replace=False)
        D = pc.pose_distances(T_all[idx], w)
        lab = pc.group(D, tol)
        for k in range(lab.max() + 1):
            mem = idx[np.flatnonzero(lab == k)]
            cen = V[mem].mean(0)
            dist = np.linalg.norm(deep_med - cen, axis=1)
            j = int(dist.argmin())
            rows.append(dict(draw=d, shallow_group=k, size=len(mem),
                             nearest_deep=j, dist=float(dist[j]),
                             matched=bool(dist[j] <= tol),
                             deep_size=int(deep_sz[j])))
    df = pd.DataFrame(rows)
    t = sout.Topic("blacksmith", "contact_saturation")
    df.to_csv(t.write("persistence", ".csv"), index=False)

    print("\n" + "=" * 80)
    print("  DO SHALLOW GROUPS SURVIVE AT DEPTH?  "
          f"{a.shallow} poses vs {len(xyz)}, tol {tol:.2f} A")
    print("=" * 80)
    print(f"\n  shallow groups per draw: {df.groupby('draw').size().mean():.0f}"
          f"   deep groups: {n_deep}")
    print(f"  shallow groups with a deep counterpart within {tol:.2f} A: "
          f"{df.matched.mean() * 100:.1f}%")
    print(f"  median centre displacement: {df.dist.median():.3f} A "
          f"(tolerance {tol:.2f} A)")
    big = df[df["size"] > 1]
    print(f"\n  restricted to NON-SINGLETON shallow groups (n={len(big)}): "
          f"{big.matched.mean() * 100:.1f}% matched, "
          f"median displacement {big.dist.median():.3f} A")
    top = df[df["size"] >= 5]
    if len(top):
        print(f"  restricted to shallow groups of >= 5 poses (n={len(top)}): "
              f"{top.matched.mean() * 100:.1f}% matched")
    print(f"\n  where the growth is: deep groups of size 1 = "
          f"{(deep_sz == 1).sum()} of {n_deep} ({(deep_sz == 1).mean() * 100:.0f}%), "
          f"size <= 3 = {(deep_sz <= 3).sum()} ({(deep_sz <= 3).mean() * 100:.0f}%)")
    print(f"  poses living in groups of >= 5: "
          f"{deep_sz[deep_sz >= 5].sum() / len(xyz) * 100:.0f}% of the cloud, "
          f"in {(deep_sz >= 5).sum()} groups")
    print("\n" + "=" * 80 + "\n")


if __name__ == "__main__":
    main()
