#!/usr/bin/env python3
"""
Purpose: is a COVERING SET reproducible across independent dockings, as modes are?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: the 5 independent replicate clouds under election_<candidate>_r{1..5}_allposes
Output: printed table (and --csv to persist)

COMPANION TO run_all.py, WHICH ASKED THE SAME QUESTION OF HDBSCAN MODES. If a
covering set is to replace mode clustering as the way representatives are chosen
(@tt8804: "define the volume of poses and partition this volume into pose
modes"), it has to survive the test the incumbent already passed: a
representative that does not come back in an independent draw of the cloud is a
partition of one sample, not a place the molecule sits.

RESOLUTION IS THE VARIABLE, and that is the point. HDBSCAN has no length scale
(D0090), so its reproducibility is whatever density happens to give. A covering
set states the scale, which means the reproducibility/count trade-off becomes
something you choose rather than something you inherit.

THE START IS THE MEDOID, NOT INDEX 0. Farthest-point traversal is deterministic
after its first centre, so the first choice decides the whole set. Taking index 0
made the construction depend on the order AutoDock happened to write its poses --
a positional choice, which is this project's defining defect shape. The medoid
(the pose minimising total distance to the rest) is a property of the CLOUD, so
two independent dockings of the same molecule start from the same place in the
pocket rather than from whichever pose was written first.

`--start first` keeps the old behaviour, because the difference between them is
itself a measurement: it says how much of the reproducibility was an artefact of
write order.
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from shared import run_paths as rp                  # noqa: E402
from shared import target_config as tc              # noqa: E402


def load(topic: str, cand: str) -> np.ndarray:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    f = rp.BLACKSMITH / f"{topic}_allposes" / f"{cand}.sdf"
    ms = [m for m in Chem.SDMolSupplier(str(f), removeHs=False, sanitize=False)
          if m is not None]
    if not ms:
        raise SystemExit(f"no poses in {f}")
    h = [i for i, a in enumerate(ms[0].GetAtoms()) if a.GetAtomicNum() > 1]
    return np.array([m.GetConformer().GetPositions()[h] for m in ms])


def medoid(c: np.ndarray) -> int:
    """The pose with the smallest total distance to every other pose."""
    n = len(c)
    tot = np.zeros(n)
    for i in range(n):
        tot += np.sqrt(((c - c[i]) ** 2).sum(axis=2).mean(axis=1))
    return int(np.argmin(tot))


def cover_set(c: np.ndarray, r: float, start: str = "medoid") -> np.ndarray:
    """Greedy farthest-point centres at radius r (the representatives)."""
    n = len(c)
    dmin = np.full(n, np.inf)
    sel, nxt = [], (medoid(c) if start == "medoid" else 0)
    while True:
        d = np.sqrt(((c - c[nxt]) ** 2).sum(axis=2).mean(axis=1))
        dmin = np.minimum(dmin, d)
        sel.append(nxt)
        far = int(np.argmax(dmin))
        if dmin[far] <= r or len(sel) >= n:
            return c[sel]
        nxt = far


def _rms(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--radii", default="")
    ap.add_argument("--start", choices=("medoid", "first"), default="medoid",
                    help="first centre: the cloud's medoid (order-independent) "
                         "or pose 0 (whatever AutoDock wrote first)")
    ap.add_argument("--csv", default="")
    a = ap.parse_args()

    tol = float(tc.get("md.sweep_survivor_rmsd_nm")) * 10.0
    radii = ([float(x) for x in a.radii.split(",")] if a.radii
             else [2.0, tol, 5.0])
    rows = []
    for r in radii:
        sets = {i: cover_set(load(f"election_{a.candidate}_r{i}", a.candidate),
                             r, a.start)
                for i in range(1, a.replicates + 1)}
        fr = []
        for i, j in itertools.combinations(range(1, a.replicates + 1), 2):
            A, B = sets[i], sets[j]
            D = np.array([[_rms(x, y) for y in B] for x in A])
            fr += [(D.min(axis=1) <= r).mean(), (D.min(axis=0) <= r).mean()]
        core = sum(1 for x in sets[1]
                   if all(any(_rms(x, y) <= r for y in sets[k])
                          for k in range(2, a.replicates + 1)))
        rows.append(dict(radius_a=r, centres=[len(v) for v in sets.values()],
                         pairwise_mean=round(float(np.mean(fr)), 4),
                         pairwise_min=round(float(np.min(fr)), 4),
                         core=core, of=len(sets[1]),
                         core_frac=round(core / len(sets[1]), 4)))
        print(f"cover @ {r:>4.1f} A: centres per replicate {rows[-1]['centres']}")
        print(f"   pairwise recovery mean {rows[-1]['pairwise_mean']*100:.1f}%  "
              f"min {rows[-1]['pairwise_min']*100:.1f}%")
        print(f"   centres of r1 found in ALL {a.replicates}: "
              f"{core}/{len(sets[1])} ({core/len(sets[1])*100:.0f}%)\n")
    print("for comparison, HDBSCAN modes (run_all.py): 88.6% pairwise, "
          "41/59 (69%) in all 5, 54-63 modes")
    if a.csv:
        pd.DataFrame(rows).to_csv(a.csv, index=False)


if __name__ == "__main__":
    main()
