#!/usr/bin/env python3
"""
Purpose: how warhead-engagement decays with rank position, and where the sweep cuts.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-27
Input: the run's rank_v2 T_4 engagement table (topic from config)
Output: 00_outputs/blacksmith/engagement_curve_<topic>/engagement_curve_<N>.png

Companion to `rank_curve_plot.py`, which plots the OLD ordering. Kept separate
rather than parameterised because the two answer different questions and a reader
must not mistake one for the other: `conditional_eb` is a frequency shrunk toward
a prior, `engagement` is a geometry score on 0-1, and D0098 measured the first at
rho = -0.015 against the MD outcome and the second at +0.652.

THE CURVE IS MONOTONE BY CONSTRUCTION and that is not a finding. Rank IS the sort
on engagement, so "score falls with rank" is arithmetic. What the plot is for is
the SHAPE -- where the cliff is, how many groups sit above it, and whether the
cap or the floor is what actually selects.

RANKED WITHIN WARHEAD CLASS, so the curves are drawn per class. An SN2 backside
attack and a perpendicular approach to an sp2 carbon do not span the same angles,
and `isotropic_null` differs between them for exactly that reason; one pooled
curve would compare geometries that are not comparable.
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

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                 # noqa: E402
from shared import run_paths as rp                 # noqa: E402
from shared import target_config as tc             # noqa: E402

log = logging.getLogger("engagement-curve")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=None)
    ap.add_argument("--zoom", type=int, default=300,
                    help="right panel x-limit: the region the sweep can reach")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    topic = a.topic or rp.topic()
    fs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_v2/rank_v2_T4_{topic}_engagement_*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit(f"no engagement ranking for topic {topic!r}")
    d = pd.read_csv(fs[-1])
    fam = tc.family_of()
    d["family"] = d.warhead_class.map(fam)
    floor, cap = tc.sweep_budget_floor(), tc.sweep_max_depth()
    log.info("%s: %d groups, %d warhead classes, %d in a swept family",
             topic, len(d), d.warhead_class.nunique(), int(d.family.notna().sum()))

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.2))
    cls = sorted(d.warhead_class.dropna().unique())
    cmap = plt.cm.tab10
    for ax, xmax, title in (
            (axes[0], None, f"All {len(d):,} groups — the whole ordering"),
            (axes[1], a.zoom, f"Top {a.zoom} per class — what the sweep can reach")):
        for i, c in enumerate(cls):
            g = d[d.warhead_class == c].sort_values("engagement", ascending=False)
            y = g.engagement.to_numpy()
            x = np.arange(1, len(y) + 1)
            if xmax:
                x, y = x[:xmax], y[:xmax]
            swept = c in fam
            ax.plot(x, y, lw=2.0 if swept else 1.0, color=cmap(i % 10),
                    alpha=1.0 if swept else 0.35,
                    label=f"{c}{'  (swept)' if swept else ''}")
        ax.axhline(floor, color="#c0392b", ls="--", lw=1.4)
        ax.text(ax.get_xlim()[1], floor, f" sweep floor {floor}", color="#c0392b",
                va="bottom", ha="right", fontsize=9)
        if xmax:
            ax.axvline(cap, color="#1a6b50", ls=":", lw=1.8)
            ax.text(cap, 1.0, f" cap {cap}/family", color="#1a6b50",
                    va="top", fontsize=9)
        ax.set(xlabel="rank within warhead class", ylabel="engagement (0–1)",
               title=title, ylim=(0, 1.02))
        if xmax is None:
            ax.set_xscale("log")
            ax.set_xlabel("rank within warhead class (log)")
        ax.grid(alpha=0.25)
    axes[1].legend(fontsize=8, loc="upper right")

    # what actually selects: the cap or the floor
    inscope = d[d.family.notna()]
    above = int((inscope.engagement >= floor).sum())
    sel = (inscope[inscope.engagement >= floor]
           .sort_values("engagement", ascending=False)
           .groupby("family").head(cap))
    fig.suptitle(
        f"{topic} — warhead engagement vs rank   ·   "
        f"{above:,} of {len(inscope):,} in-scope groups clear the {floor} floor, "
        f"but the {cap}/family cap takes {len(sel)} "
        f"(engagement {sel.engagement.min():.2f}–{sel.engagement.max():.2f})",
        fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    t = sout.Topic("blacksmith", f"engagement_curve_{topic}")
    out = t.write("engagement_curve", ".png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(out)
    print(f"\n  THE CAP SELECTS, NOT THE FLOOR: {above:,} groups clear {floor}; "
          f"the cap admits {len(sel)}.")
    print(f"  The floor is doing no work here -- the lowest SELECTED engagement is "
          f"{sel.engagement.min():.3f}, {sel.engagement.min()/floor:.0f}x the floor.")


if __name__ == "__main__":
    main()
