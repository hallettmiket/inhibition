#!/usr/bin/env python3
"""
Purpose: does the occupied region of residue-contact space grow with pose count?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: the 6,000-pose RAW deep cloud
Output: 00_outputs/blacksmith/contact_saturation/

WHAT IS ENTAILED, STATED FIRST. Every coordinate of contact space is a distance
capped at pose_contacts.CAP_A (10 A) and floored at van der Waals contact, so the
region is bounded BY CONSTRUCTION before a single pose is docked. D0091 is the
record of docking a cloud into a 26 A box, measuring that it stayed inside 26 A,
and reporting the bound as a finding. Reporting "contact space is bounded" would
be the same error in a different coordinate system.

So the finding is never the bound. It is:

  EXTENT      -- how fast the occupied diameter approaches its ceiling.
  DIMENSION   -- the covering number N(eps) ~ eps^-d gives the effective
                 dimension the poses actually occupy, out of the 420 the
                 coordinates offer. A rigid ligand has 6 rigid-body degrees of
                 freedom plus torsions, so d near 6 would say the metric is
                 tracking pose and not noise; d near 420 would say it is
                 tracking noise.
  FILL        -- n against N(eps). The group count climbs (b = +0.69) for one of
                 two reasons: the region is still growing, or the region is fixed
                 and 6,000 poses is far too few to cover it. These are different
                 problems with different fixes and this separates them.
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

log = logging.getLogger("space-growth")
_spec = importlib.util.spec_from_file_location(
    "x17", Path(__file__).with_name("run_all.py"))
_M = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_M)

LADDER = [100, 200, 350, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000]
EPS = [0.73, 1.0, 1.5, 2.0, 2.5, 3.0]


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
    rmsf = _M.predict_rmsf(meta[0], meta[1], 50, a.seed)
    w = pc.atom_weights(rmsf)
    T = pc.contact_tensor(xyz, res)
    D = pc.pose_distances(T, w)
    log.info("cloud %d poses, contact tensor %s, cap %.1f A", len(xyz), T.shape, pc.CAP_A)

    rng = np.random.default_rng(a.seed)
    rows = []
    for k in [x for x in LADDER if x <= len(xyz)]:
        for d in range(1 if k == len(xyz) else a.draws):
            idx = (np.arange(len(xyz)) if k == len(xyz)
                   else rng.choice(len(xyz), size=k, replace=False))
            sub = D[np.ix_(idx, idx)]
            iu = np.triu_indices(len(idx), 1)
            r = dict(poses=k, draw=d, diameter=float(sub[iu].max()),
                     mean_dist=float(sub[iu].mean()),
                     p99=float(np.percentile(sub[iu], 99)))
            for e in EPS:
                r[f"N_{e}"] = int(pc.group(sub, e).max() + 1)
            rows.append(r)
        log.info("  n=%5d done", k)
    df = pd.DataFrame(rows)
    t = sout.Topic("blacksmith", "contact_saturation")
    df.to_csv(t.write("space_growth", ".csv"), index=False)

    agg = df.groupby("poses").mean(numeric_only=True)

    # effective dimension from the covering-number slope, at each depth
    def dim_at(row) -> float:
        n = np.array([row[f"N_{e}"] for e in EPS], float)
        ok = n > 2
        if ok.sum() < 3:
            return float("nan")
        A = np.vstack([np.ones(ok.sum()), np.log(np.array(EPS)[ok])]).T
        c, *_ = np.linalg.lstsq(A, np.log(n[ok]), rcond=None)
        return float(-c[1])

    print("\n" + "=" * 80)
    print("  DOES CONTACT SPACE ITSELF GROW?   (bounded by the 10 A cap "
          "BY CONSTRUCTION)")
    print("=" * 80)
    print(f"\n  EXTENT — is the region still expanding?\n")
    print(f"    {'poses':>7} {'diameter':>10} {'99th pct':>10} {'mean':>8} "
          f"{'diam vs n=100':>14}")
    d0 = agg.diameter.iloc[0]
    for k, r in agg.iterrows():
        print(f"    {k:7,} {r.diameter:9.2f}A {r.p99:9.2f}A {r.mean_dist:7.2f}A "
              f"{r.diameter / d0:13.2f}x")
    b_d, _ = _M.fit_power(agg.index.values, agg.diameter.values)
    b_m, _ = _M.fit_power(agg.index.values, agg.mean_dist.values)
    print(f"\n    diameter exponent {b_d:+.3f}, mean-distance exponent {b_m:+.3f} "
          f"(0.0 = fixed region)")
    print(f"    60x more poses widened the diameter by "
          f"{agg.diameter.iloc[-1] / d0:.2f}x and moved the MEAN pose separation "
          f"by {agg.mean_dist.iloc[-1] / agg.mean_dist.iloc[0]:.2f}x")

    print(f"\n  DIMENSION — how many degrees of freedom the poses actually use\n")
    print(f"    {'poses':>7} {'effective d':>12}   covering numbers N(eps)")
    for k, r in agg.iterrows():
        ns = "  ".join(f"{e}:{r[f'N_{e}']:.0f}" for e in EPS)
        print(f"    {k:7,} {dim_at(r):12.2f}   {ns}")
    print(f"\n    coordinates available: {T.shape[1] * T.shape[2]} "
          f"({T.shape[1]} atoms x {T.shape[2]} residues)")
    print(f"    a rigid ligand has 6 rigid-body degrees of freedom; "
          f"{meta[0].GetNumAtoms()} atoms with torsions add more")

    print(f"\n  FILL — is the region growing, or just not yet covered?\n")
    top = agg.iloc[-1]
    dd = dim_at(top)
    for e in EPS:
        N = top[f"N_{e}"]
        # poses needed to cover a d-dimensional region at this scale, from the
        # coupon-collector-like scaling N * ln(N)
        need = N * np.log(max(N, 2))
        print(f"    eps={e:4.2f} A: {N:6.0f} cells occupied by 6,000 poses "
              f"-> ~{need:,.0f} poses to cover what is already visible")
    print(f"\n    the region is {'FIXED' if abs(b_d) < 0.05 else 'still growing'} "
          f"(diameter exponent {b_d:+.3f}); the climbing group count is therefore "
          f"{'undersampling of a fixed region' if abs(b_d) < 0.05 else 'genuine expansion'}")
    print("\n" + "=" * 80 + "\n")


if __name__ == "__main__":
    main()
