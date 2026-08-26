#!/usr/bin/env python3
"""
Purpose: at what tolerance does the contact-space group count stop climbing?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: the 6,000-pose RAW deep cloud
Output: 00_outputs/blacksmith/contact_saturation/

run_all.py measured ONE tolerance -- the molecule's predicted RMSF, 0.73 A -- and
found the count still climbing as n^0.70 at 6,000 poses. That is a fact about
0.73 A, not about the method. The tolerance is the only free parameter and it
sets the covering number directly, so the question "does the count taper" has no
answer until it is asked across the range.

THE RAW CLOUD, DELIBERATELY. `<topic>_allposes` is not all poses: nac_screen_v2
writes only poses whose DBSCAN label is in `mode_ids`, so ~21% of every
production cloud -- the scattered poses that failed to join a mode -- is already
gone. Measuring a replacement for DBSCAN on clouds DBSCAN has already cleaned
answers an easier question than the one asked. exp/14, exp/15 and exp/16 all
read that file and inherit the same caveat.

WHAT IS ENTAILED AND THEREFORE NOT A FINDING. A fixed absolute tolerance bounds
the count by the covering number at that scale, so SOME tolerance must flatten
the curve; D0091 is the record of reporting an apparatus-entailed bound as a
result. The finding here is the RATE -- which tolerance flattens it, and whether
that tolerance is one the downstream stages can tell apart.
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
from shared import target_config as tc              # noqa: E402

log = logging.getLogger("tolerance-sweep")

_spec = importlib.util.spec_from_file_location(
    "x17", Path(__file__).with_name("run_all.py"))
_M = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_M)

LADDER = [100, 200, 350, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--residues", type=int, default=15)
    ap.add_argument("--draws", type=int, default=3)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    res = _M.receptor_coords(_M.key_residues(a.residues))
    xyz, meta = _M.load_sdf(rp.BLACKSMITH / f"deep_cloud_{a.candidate}" / "cloud_1.sdf")
    if xyz is None:
        raise SystemExit("no deep cloud")
    rmsf = _M.predict_rmsf(meta[0], meta[1], 50, a.seed)
    w = pc.atom_weights(rmsf)
    rmsf_tol = float(np.median(rmsf) / pc.RMSF_CALIBRATION)
    sweep_bar = float(tc.md_survivor_rmsd_nm()) * 10.0

    tols = sorted({round(rmsf_tol, 3), 1.0, 1.5, 2.0, 2.5, 3.0, round(sweep_bar, 3)})
    log.info("cloud %d poses; tolerances %s (RMSF %.2f, sweep bar %.2f)",
             len(xyz), tols, rmsf_tol, sweep_bar)

    rng = np.random.default_rng(a.seed)
    rows = []
    for k in [x for x in LADDER if x <= len(xyz)]:
        nd = 1 if k == len(xyz) else a.draws
        for d in range(nd):
            idx = (np.arange(len(xyz)) if k == len(xyz)
                   else rng.choice(len(xyz), size=k, replace=False))
            sub = xyz[idx]
            # ONE distance matrix, every tolerance. The metric does not depend on
            # the cut, so recomputing it per tolerance would be the same numbers
            # seven times.
            D = pc.pose_distances(pc.contact_tensor(sub, res), w)
            for t_ in tols:
                lab = pc.group(D, t_)
                sizes = np.bincount(lab)
                rows.append(dict(tol=t_, poses=k, draw=d, groups=len(sizes),
                                 largest=int(sizes.max()),
                                 singletons=int((sizes == 1).sum())))
        log.info("  n=%5d done", k)

    df = pd.DataFrame(rows)
    t = sout.Topic("blacksmith", "contact_saturation")
    df.to_csv(t.write("tolerance_sweep", ".csv"), index=False)

    print("\n" + "=" * 84)
    print("  WHICH TOLERANCE MAKES THE COUNT STOP CLIMBING?  raw 6,000-pose cloud")
    print("=" * 84)
    print(f"\n  {'tol (A)':>8} {'n=500':>8} {'n=6000':>8} {'x growth':>9} "
          f"{'exponent b':>11} {'plateau':>9} {'n95':>10} {'largest':>8} {'singl.':>7}")
    for t_ in tols:
        s = df[df.tol == t_]
        agg = s.groupby("poses")["groups"].mean()
        b, _ = _M.fit_power(agg.index.values, agg.values)
        G, K, n95, _ = _M.fit_saturating(agg.index.values, agg.values)
        g500, g6k = agg.get(500, np.nan), agg.iloc[-1]
        top = s[s.poses == agg.index[-1]]
        plateau = "none" if not np.isfinite(G) or G < g6k else f"{G:,.0f}"
        n95s = "-" if plateau == "none" else f"{n95:,.0f}"
        mark = "  <- RMSF" if abs(t_ - rmsf_tol) < 1e-6 else (
            "  <- sweep bar" if abs(t_ - sweep_bar) < 1e-6 else "")
        print(f"  {t_:8.2f} {g500:8.0f} {g6k:8.0f} {g6k / g500:8.2f}x "
              f"{b:+11.3f} {plateau:>9} {n95s:>10} {top.largest.mean():8.0f} "
              f"{top.singletons.mean() / g6k * 100:6.0f}%{mark}")
    print("\n  b = 1.0 every new pose is a new group; b = 0.0 the count is flat.")
    print("  plateau/n95 from the species-accumulation fit; 'none' means the fit "
          "implies\n  no finite ceiling within the measured range.")
    print("\n" + "=" * 84)
    print(f"  written to {t.dir}")
    print("=" * 84 + "\n")


if __name__ == "__main__":
    main()
