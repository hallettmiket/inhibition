"""
Purpose: recompute a size-stable anchoring statistic for every existing mode, without re-docking.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: 00_outputs/blacksmith/nac_v3/poses_*.csv (every pose's distance + angle)
Output: 00_outputs/blacksmith/nac_v3/anchor_p90_<N>.csv — one row per mode

WHY THIS CAN BE DONE AT ALL. `anchor_quality` is a pure function of a pose's
distance, angle and mechanism (`shared/nac_criterion.py`), and the 2.2.0 screen
persisted all three for **every** pose. So the statistic can be recomputed from
what is on disk — no docking, no GPU, and the numbers correspond to exactly the
poses the current scores were computed from. That last property is what makes
this legitimate where a re-dock would not be: a re-docked pose cloud is a
different sample and cannot sit beside the existing scores.

THE DEFECT BEING REPAIRED (#43). `anchor_quality_max` is the maximum over a
mode's poses, so it grows with how many poses the mode has. Measured on T4:

    rho(anchor_quality_max,  mode size) = +0.740
    rho(anchor_quality_mean, mode size) = +0.445
    rho(enrichment_conditional, size)   = -0.021

Median max was 0.828 for mode 0 against 0.019 for mode 1 — a 40x gap that is
mostly the number of draws. That matters beyond display: the max decides which
extra modes get swept (`weekend_worklist.py`'s MODE_GAIN rule) and which pose
represents a mode, so the bias reaches selection, not just the table.

A QUANTILE ESTIMATES THE SAME PROPERTY AND CONVERGES. "How good is this mode's
well-anchored tail" is a question about the distribution, and the 90th percentile
answers it with an estimate whose expectation does not drift with n. The max
answers a question about the sample instead of the population.

Reported alongside, never overwriting: `anchor_quality_max` stays on disk so the
current ranking remains reproducible and the two can be compared.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import nac_criterion as nac              # noqa: E402
from shared import outputs as sout                   # noqa: E402

log = logging.getLogger("anchor-p90")
DATA = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default="nac_v3")
    ap.add_argument("--quantile", type=float, default=90.0)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from scipy.stats import spearmanr

    fs = sorted(glob.glob(str(DATA / args.topic / "poses_*.csv")))
    if not fs:
        raise SystemExit(f"no per-pose files under {DATA / args.topic}")
    log.info("reading %d per-pose shard files", len(fs))
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    log.info("%d poses over %d molecules", len(d), d.ident.nunique())

    # Poses labelled -1 belong to no mode; they are noise by construction and
    # must not contribute to any mode's statistic.
    d = d[d["mode"] >= 0].copy()

    # Vectorised over the frame rather than per row: `anchor_quality` is scalar,
    # and 2.9M Python-level calls would take longer than the docking did.
    aq = np.array([nac.anchor_quality(dd, aa, mm)
                   for dd, aa, mm in zip(d.distance.values, d.angle.values,
                                         d.mechanism.values)])
    d["anchor"] = aq

    g = d.groupby(["ident", "mode"])
    out = g.agg(n_poses_mode=("anchor", "size"),
                anchor_quality_max=("anchor", "max"),
                anchor_quality_mean=("anchor", "mean")).reset_index()
    out[f"anchor_quality_p{int(args.quantile)}"] = (
        g["anchor"].quantile(args.quantile / 100.0).values)
    out["parent_ident"] = out.ident
    out["ident"] = out.ident + "_m" + out["mode"].astype(int).astype(str)

    col = f"anchor_quality_p{int(args.quantile)}"
    print("\n" + "=" * 72)
    print(f"  size dependence — rho(statistic, mode size), 0 is unbiased")
    print("=" * 72)
    for c in ("anchor_quality_max", "anchor_quality_mean", col):
        r, _ = spearmanr(out[c], out.n_poses_mode, nan_policy="omit")
        print(f"    {c:26s} {r:+.3f}")
    print("\n  median by mode:")
    out["_m"] = out["mode"].clip(upper=3)
    print(out.groupby("_m")[["n_poses_mode", "anchor_quality_max", col]]
             .median().round(3).to_string())

    dest = sout.Topic("blacksmith", args.topic).write("anchor_p90", ".csv")
    out.drop(columns=["_m"]).to_csv(dest, index=False)
    print(f"\n  {len(out):,} modes -> {dest}")
    print("  `anchor_quality_max` is preserved so the current ranking stays "
          "reproducible.\n")


if __name__ == "__main__":
    main()
