#!/usr/bin/env python3
"""
Purpose: plot the exp/17 saturation, tolerance, extent and persistence results
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: 00_outputs/blacksmith/contact_saturation/*.csv (written by the four runners)
Output: 00_outputs/blacksmith/contact_saturation/figures_*.png

THE PANEL THAT CARRIES THE ARGUMENT IS THE THIRD ONE. Panels 1 and 2 show a count
that climbs and never flattens, which read alone says the method is unbounded.
Panel 3 puts the region's extent on the same x-axis and it is a horizontal line.
Shown apart, the two invite opposite conclusions; the contrast IS the finding, so
they are drawn together and the extent panel shares the pose axis.

AXES ARE PINNED TO WHAT THEY MEAN. The extent panel starts at zero, because a
y-axis auto-scaled to a flat series turns +0.019 into a dramatic slope -- the same
figure would then argue the opposite of what it measures.
"""

from __future__ import annotations

import glob
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from shared import outputs as sout                  # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("plots")
DIR = rp.BLACKSMITH / "contact_saturation"


def newest(stem: str) -> pd.DataFrame:
    fs = sorted(glob.glob(str(DIR / f"{stem}_*.csv")), key=os.path.getmtime)
    if not fs:
        raise SystemExit(f"no {stem}_*.csv in {DIR} -- run the exp/17 scripts first")
    return pd.read_csv(fs[-1])


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    deep, shallow = newest("deep_ladder"), newest("shallow_ladders")
    tol, space, pers = newest("tolerance_sweep"), newest("space_growth"), newest("persistence")

    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 3, hspace=0.33, wspace=0.38)
    C = plt.cm.viridis

    # ---- 1. the ladder ---------------------------------------------------- #
    ax = fig.add_subplot(gs[0, 0])
    for ident, s in shallow.groupby("ident"):
        m = s.groupby("poses")["groups"].mean()
        ax.plot(m.index, m.values, "-", color="0.75", lw=0.9, zorder=1)
    m = deep.groupby("poses")["groups"].mean()
    ax.plot(m.index, m.values, "o-", color="#1f4e79", lw=2, ms=5, zorder=3,
            label="deep cloud (raw, 6,000)")
    ax.plot([], [], "-", color="0.75", lw=0.9, label="12 production molecules")
    n = np.array([100, 6000])
    ax.plot(n, m.iloc[0] * (n / 100) ** 1.0, "--", color="#c0392b", lw=1.2,
            label="b = 1.0 (every pose a new group)")
    ax.plot(n, [m.iloc[0]] * 2, ":", color="#27ae60", lw=1.4, label="b = 0.0 (flat)")
    ax.set(xscale="log", yscale="log", xlabel="poses docked", ylabel="groups")
    ax.set_title("1. The count never flattens\nb = +0.69 deep, +0.67 median across molecules",
                 fontsize=10.5, loc="left")
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
    ax.grid(alpha=0.25, which="both")

    # ---- 2. tolerance sweep ------------------------------------------------ #
    ax = fig.add_subplot(gs[0, 1])
    tols = sorted(tol.tol.unique())
    for i, t_ in enumerate(tols):
        m = tol[tol.tol == t_].groupby("poses")["groups"].mean()
        ax.plot(m.index, m.values, "o-", ms=3.5, lw=1.5,
                color=C(i / max(len(tols) - 1, 1)), label=f"{t_:.2f} Å")
    ax.set(xscale="log", yscale="log", xlabel="poses docked", ylabel="groups")
    ax.set_title("2. No tolerance flattens it\nlooser only lowers the intercept",
                 fontsize=10.5, loc="left")
    ax.legend(fontsize=7, title="tolerance", title_fontsize=7, ncol=2)
    ax.grid(alpha=0.25, which="both")

    # ---- 3. THE CONTRAST --------------------------------------------------- #
    ax = fig.add_subplot(gs[0, 2])
    sp = space.groupby("poses").mean(numeric_only=True)
    ax.plot(sp.index, sp.diameter, "o-", color="#c0392b", lw=2, ms=4.5, label="diameter")
    ax.plot(sp.index, sp.p99, "s-", color="#e67e22", lw=1.6, ms=4, label="99th pct")
    ax.plot(sp.index, sp.mean_dist, "^-", color="#8e44ad", lw=1.6, ms=4, label="mean")
    ax.set(xscale="log", xlabel="poses docked",
           ylabel="contact-space distance (Å)", ylim=(0, sp.diameter.max() * 1.35))
    ax.set_title("3. …but the REGION is fixed\ndiameter exponent +0.019, mean +0.001",
                 fontsize=10.5, loc="left")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.25, which="both")
    ax2 = ax.twinx()
    ax2.plot(m.index if False else sp.index,
             [tol[(tol.tol == tols[0]) & (tol.poses == k)].groups.mean() for k in sp.index],
             "--", color="0.45", lw=1.3)
    ax2.set_ylabel("groups (dashed, right axis)", color="0.45", fontsize=8.5)
    ax2.tick_params(axis="y", colors="0.45", labelsize=8)

    # ---- 4. the trade-off -------------------------------------------------- #
    ax = fig.add_subplot(gs[1, 0])
    top = tol[tol.poses == tol.poses.max()].groupby("tol").mean(numeric_only=True)
    ax.plot(top.index, top.groups, "o-", color="#1f4e79", lw=2, ms=6, label="groups")
    ax.set(xlabel="tolerance (Å)", ylabel="groups at n = 6,000", yscale="log")
    ax.grid(alpha=0.25)
    ax3 = ax.twinx()
    ax3.plot(top.index, top.largest, "s--", color="#c0392b", lw=2, ms=6)
    ax3.set_ylabel("largest group (poses)", color="#c0392b")
    ax3.tick_params(axis="y", colors="#c0392b")
    ax3.axhline(137, color="#c0392b", ls=":", lw=1.2)
    ax3.text(top.index.max(), 150, "D0088's condemned bag (137)", fontsize=7,
             color="#c0392b", ha="right")
    ax.set_title("4. The trade-off has no sweet spot\nfewer groups only by rebuilding the bag",
                 fontsize=10.5, loc="left")

    # ---- 5. persistence ---------------------------------------------------- #
    ax = fig.add_subplot(gs[1, 1])
    t_used = float(pers.dist.max() * 0 + 0.729)
    ax.hist(pers.dist, bins=45, color="#1f4e79", alpha=0.85)
    ax.axvline(t_used, color="#c0392b", lw=2,
               label=f"tolerance {t_used:.2f} Å")
    ax.axvline(pers.dist.median(), color="#27ae60", lw=2, ls="--",
               label=f"median {pers.dist.median():.3f} Å")
    ax.set(xlabel="displacement of a shallow group's centre at depth (Å)",
           ylabel="shallow groups")
    ax.set_title(f"5. Groups never move\n100% of {len(pers):,} n=500 groups persist at n=6,000",
                 fontsize=10.5, loc="left")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    # ---- 6. where the growth is -------------------------------------------- #
    ax = fig.add_subplot(gs[1, 2])
    for i, k in enumerate([500, 1000, 3000, 6000]):
        s = tol[(tol.tol == tols[0]) & (tol.poses == k)]
        if s.empty:
            continue
        frac_single = s.singletons.mean() / s.groups.mean() * 100
        ax.bar(i, frac_single, color=C(i / 3), width=0.62)
        ax.text(i, frac_single + 1.4, f"{frac_single:.0f}%", ha="center", fontsize=9)
    ax.set(xticks=range(4), xticklabels=["500", "1,000", "3,000", "6,000"],
           xlabel="poses docked", ylabel="% of groups holding one pose",
           ylim=(0, 70))
    ax.set_title("6. Growth is all tail\nsingletons FALL with depth — the cloud consolidates",
                 fontsize=10.5, loc="left")
    ax.grid(alpha=0.25, axis="y")

    fig.suptitle("exp/17 — contact-space grouping: the count climbs, the region does not "
                 "(t4_716800c125a7, raw cloud)", fontsize=13, y=0.975)
    t = sout.Topic("blacksmith", "contact_saturation")
    out = t.write("figures", ".png")
    fig.savefig(out, dpi=155, bbox_inches="tight")
    log.info("wrote %s", out)
    print(out)


if __name__ == "__main__":
    main()
