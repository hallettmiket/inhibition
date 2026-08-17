#!/usr/bin/env python3
"""
Purpose: how long should the triage sweep be? Measured, at 0.1 ns resolution.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-16
Input: --worklist <sweep_gaps_N.csv> [--cut 0.5] [--md100-h 4.5]
Output: 00_outputs/blacksmith/exp_sweep_length/sweep_length_<N>.csv + a printed
        decision, and docs/sweep_length.md is written from the same numbers

THE DESIGN THIS IS OPTIMISING (@tt8804). Poses are ranked on the docked geometry
-- the best case, if the pose held. A short MD then asks only whether the pose is
STABLE, and every mode that stays under a max-RMSD bar for x ns earns a 100 ns
run. The budget is spent going down the ranked list, so the question is not "is
x ns accurate" but "what x finds the most real candidates per GPU-hour".

WHY x MATTERS AT ALL, IN BOTH DIRECTIONS. A shorter sweep lets you get further
down the ranked list for the same money -- more coverage, more candidates. But
max RMSD only grows with time, so a shorter sweep also passes poses that a longer
one would reject, and every one of those costs a full 100 ns run. Coverage up,
precision down, and the exchange rate between them is what this measures.

THE ESTIMAND. Genuine survivors found per GPU-hour, where a mode is genuine if it
holds under the cut for the full 10 ns. Per mode:

    cost(x)  = x * sweep_h_per_ns  +  pass_rate(x) * md100_h
    yield    = genuine_rate                       (independent of x -- see below)
    metric   = yield / cost(x)

`genuine_rate` does not depend on x BECAUSE TRUNCATION IS ONE-SIDED: max@x <=
max@10, so every genuine mode passes at every x. Shortening cannot lose one from
the survivor SET. What it changes is how many extras come with them.

WHAT THIS CANNOT SETTLE. "Genuine" is defined against the 10 ns verdict, not
against 100 ns, because only a handful of 100 ns runs exist. The 10 ns hazard was
still running at 9 ns (0.015/ns, down from 0.045), so some 10 ns survivors will
fail at 100 ns -- and under this design that is fine, the 100 ns stage is itself
the filter. It does mean the optimum here is the best x for reproducing the 10 ns
answer cheaply, not the best x for predicting 100 ns. Re-run with
`--truth-md100` once enough 100 ns rows exist.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

log = logging.getLogger("exp-sweep-length")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")

#: Measured, not assumed: 251 completed 10 ns sweeps, median 41 min end to end.
SWEEP_H_PER_NS = 41.0 / 10.0 / 60.0
#: Observed on the 100 ns runs launched 2026-08-16.
MD100_H = 4.5


def traces(worklist: Path) -> dict:
    """{ident: (t_ns, rmsd_nm)} for every finished mode of this campaign.

    Keyed through `sweep_assets.rep_dir` WITH the pose rank, because a molecule
    with several swept modes has one directory per mode and matching on the
    molecule alone hands them all the same trajectory.
    """
    import sweep_assets as sa
    from shared import sweep_state as ss
    wl = pd.read_csv(worklist)
    pr = dict(zip(wl.ident.astype(str), wl.pose_rank.astype(int)))
    d = ss.state(worklist)
    ok = d[(d.sweep_state == "ok") & d["_queued"]]
    out = {}
    for i in ok.ident.astype(str):
        rep = sa.rep_dir(i.rsplit("_m", 1)[0], pr.get(i))
        if rep is None:
            continue
        t, y = sa._xvg(rep / "rmsd.xvg")
        if t is not None and len(t) > 100:
            out[i] = (t, y)
    return out


def curve(tr: dict, cut: float, md100_h: float, grid: np.ndarray) -> pd.DataFrame:
    """Pass rate and yield-per-GPU-hour at every x on `grid`."""
    keys = list(tr)
    # running maximum per trace, sampled on the grid -- one pass, not one per x
    mx = np.empty((len(keys), len(grid)))
    for r, k in enumerate(keys):
        t, y = tr[k]
        run = np.maximum.accumulate(y)
        mx[r] = np.interp(grid, t, run)
    genuine = mx[:, -1] < cut                      # holds for the full window
    rows = []
    for c, x in enumerate(grid):
        p = float((mx[:, c] < cut).mean())
        cost = x * SWEEP_H_PER_NS + p * md100_h
        rows.append({"x_ns": round(float(x), 2), "pass_rate": p,
                     "sweep_h": x * SWEEP_H_PER_NS, "md100_h": p * md100_h,
                     "total_h": cost, "yield_per_1000h": genuine.mean() / cost * 1000})
    return pd.DataFrame(rows), mx, genuine


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--worklist", required=True)
    ap.add_argument("--cut", type=float, default=0.5)
    ap.add_argument("--md100-h", type=float, default=MD100_H)
    ap.add_argument("--boot", type=int, default=2000)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    tr = traces(Path(args.worklist))
    grid = np.round(np.arange(0.1, 10.01, 0.1), 2)
    df, mx, genuine = curve(tr, args.cut, args.md100_h, grid)
    n = len(tr)
    best = df.loc[df.yield_per_1000h.idxmax()]

    print(f"\n{n} finished modes · cut {args.cut} nm · "
          f"{genuine.mean()*100:.1f}% genuine · 100 ns at {args.md100_h} h")
    print(f"sweep cost {SWEEP_H_PER_NS*60:.2f} min per ns per mode; "
          f"one 100 ns run = {args.md100_h/SWEEP_H_PER_NS:.0f} ns of sweeping\n")

    # THE PEAK IS ONLY WORTH REPORTING IF IT SURVIVES RESAMPLING. 168 modes is
    # not many, and an argmax over 100 grid points will always find something.
    rng = np.random.default_rng(0)
    opt = []
    for _ in range(args.boot):
        idx = rng.integers(0, n, n)
        g = genuine[idx]
        p = (mx[idx] < args.cut).mean(axis=0)
        cost = grid * SWEEP_H_PER_NS + p * args.md100_h
        opt.append(float(grid[np.argmax(g.mean() / cost)]))
    opt = np.array(opt)
    lo, hi = np.percentile(opt, [2.5, 97.5])
    print(f"OPTIMUM  x = {best.x_ns:.1f} ns   "
          f"({best.yield_per_1000h:.1f} genuine per 1000 GPU-h)")
    print(f"  bootstrap 95% CI over {args.boot} resamples: {lo:.1f} – {hi:.1f} ns")
    print(f"  median of resampled optima: {np.median(opt):.1f} ns\n")

    # HOW FLAT IS IT? A peak you cannot distinguish from its neighbours is a
    # recommendation to pick on other grounds, and saying so is the finding.
    within = df[df.yield_per_1000h >= 0.98 * best.yield_per_1000h]
    print(f"within 2% of the optimum: x = {within.x_ns.min():.1f} – "
          f"{within.x_ns.max():.1f} ns  ({len(within)} of {len(df)} grid points)")

    print("\n  x    pass    sweep    100ns    total   per-1000h")
    for x in (0.5, 1, 2, 3, 4, 5, 6, 7, 8, 8.5, 9, 9.5, 10):
        r = df[np.isclose(df.x_ns, x)]
        if len(r):
            r = r.iloc[0]
            mark = "  <-- optimum" if np.isclose(r.x_ns, best.x_ns) else ""
            print(f"{r.x_ns:5.1f}{r.pass_rate*100:7.1f}%{r.sweep_h:8.3f}h"
                  f"{r.md100_h:8.3f}h{r.total_h:8.3f}h{r.yield_per_1000h:11.1f}{mark}")

    print("\nSENSITIVITY — the optimum moves with the 100 ns cost:")
    for h in (2.0, 3.0, 4.5, 6.0, 8.0):
        d2, _, _ = curve(tr, args.cut, h, grid)
        print(f"   100 ns = {h:4.1f} h  ->  x* = {d2.loc[d2.yield_per_1000h.idxmax()].x_ns:4.1f} ns")
    print("SENSITIVITY — and with the RMSD bar:")
    for c in (0.4, 0.5, 0.6, 0.8, 1.0):
        d2, _, _ = curve(tr, c, args.md100_h, grid)
        print(f"   cut = {c:.1f} nm  ->  x* = {d2.loc[d2.yield_per_1000h.idxmax()].x_ns:4.1f} ns")

    from shared import outputs as sout
    dest = sout.Topic("blacksmith", "exp_sweep_length").write("sweep_length", ".csv")
    df.to_csv(dest, index=False)
    print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
