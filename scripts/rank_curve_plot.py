#!/usr/bin/env python3
"""
Purpose: how the ranking score decays with rank position, and where the sweep floor cuts.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: the run's rank_v2 T_4 table (topic from config)
Output: 00_outputs/blacksmith/rank_curve_<topic>/rank_curve_<N>.png

@tt8804: "generate a graph to show how ranking score is related to ranking place?
I want to see where this curve drops off."

TWO PANELS, NOT TWO Y-AXES. The ordering is by `conditional_eb`; the sweep floor
cuts on `enrichment`. They are different quantities on different scales, so they
get one panel each sharing the x -- a dual-axis chart would invite reading a
crossing that does not exist.

THE ENRICHMENT PANEL IS POINTS, NOT A LINE, and that is the finding rather than
a style choice. The list is ordered by `conditional_eb`, so `enrichment` is NOT
monotonic along it -- a line would draw a smooth decay over data that zigzags
between 0.5 and 9 at adjacent ranks, and would imply the floor cuts a PREFIX of
the ranked list. It does not: it cuts a scattered subset. A rolling median is
drawn over the points so the trend is still readable.

LINEAR X BY DEFAULT (@tt8804). `--xscale log` is kept because the decay is in the
first ~50 ranks and the tail runs to 774, so a linear axis gives the part being
asked about ~6% of the width -- but linear is what the eye reads as "how many
modes", and that is the question the floor is answering.

Colours are the first three categorical slots of the reference palette, in fixed
order by class name -- never cycled, and never reassigned when a filter changes
the series count.
"""

from __future__ import annotations

import glob
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import run_paths as rp                  # noqa: E402
from shared import target_config as tc              # noqa: E402

log = logging.getLogger("rank-curve")

#: Reference-palette categorical slots 1-3, assigned in fixed order.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, INK_2, MUTED = "#0b0b0b", "#3a3a3a", "#6b6b6b"
SURFACE, GRID, RULE = "#fcfcfb", "#e8e8e6", "#cfcfcc"


def load() -> pd.DataFrame:
    topic = tc.topic() if hasattr(tc, "topic") else tc.get("run.topic")
    fs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_v2/rank_v2_T4_{topic}_conditional_eb_*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit(f"no ranked table for topic {topic}")
    d = pd.read_csv(fs[-1])
    return d[d.class_rank.notna()].copy(), topic, os.path.basename(fs[-1])


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--xscale", choices=("linear", "log"), default="linear")
    ap.add_argument("--xmax", type=int, default=0,
                    help="clip the rank axis; 0 = the full class")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d, topic, src = load()
    floor = float(tc.sweep_budget_floor())
    classes = sorted(d.warhead_class.unique())
    colour = {c: SERIES[i] for i, c in enumerate(classes)}

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 8.4), sharex=True,
                             facecolor=SURFACE,
                             gridspec_kw={"hspace": 0.16})
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(True, which="major", color=GRID, linewidth=0.8, zorder=0)
        ax.grid(True, which="minor", color=GRID, linewidth=0.4, alpha=0.6, zorder=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(RULE)
        ax.tick_params(colors=INK_2, labelsize=9)

    panels = [("conditional_eb", "conditional_eb  —  what the ranking orders on", axes[0]),
              ("enrichment", "enrichment  —  what the sweep floor cuts on", axes[1])]

    cross = {}
    for col, title, ax in panels:
        for c in classes:
            g = d[d.warhead_class == c].sort_values("class_rank")
            if col == "enrichment":
                ax.plot(g.class_rank, g[col], ".", color=colour[c], markersize=3.4,
                        alpha=0.55, zorder=3, label=c, markeredgewidth=0)
                med = g[col].rolling(25, center=True, min_periods=8).median()
                ax.plot(g.class_rank, med, color=colour[c], linewidth=2.0,
                        solid_capstyle="round", zorder=4)
            else:
                ax.plot(g.class_rank, g[col], color=colour[c], linewidth=2.0,
                        solid_capstyle="round", zorder=3, label=c)
                # DIRECT LABEL WHERE THE LINES ARE SEPARATED, not at the end --
                # all three converge toward 0 and end-labels overprinted.
                # STAGGERED along x, because two of the three run within
                # 0.4 units of each other and labels at one x overprinted.
                at = ([7, 20, 55] if args.xscale == "log"
                      else [90, 200, 330])[classes.index(c)]
                at = min(at, len(g))
                x0, y0 = g.class_rank.iloc[at - 1], g[col].iloc[at - 1]
                ax.plot([x0], [y0], "o", color=colour[c], markersize=5, zorder=5)
                ax.annotate(f" {c}", (x0, y0), color=INK_2, fontsize=9,
                            va="bottom", ha="left", xytext=(7, 4),
                            textcoords="offset points", zorder=5)
            if col == "enrichment":
                cross[c] = int((g[col] >= floor).sum())
        ax.set_title(title, color=INK, fontsize=11, loc="left", pad=8)
        ax.set_xscale(args.xscale)
        hi = args.xmax if args.xmax else int(d.class_rank.max())
        ax.set_xlim(1, hi * (2.3 if args.xscale == "log" else 1.10))

    ax = axes[1]
    ax.axhline(floor, color=MUTED, linewidth=1.4, linestyle=(0, (5, 3)), zorder=2)
    # Right-hand end: the left edge is where every series is still high and the
    # label sat on top of the bdhi_c5 line.
    ax.annotate(f"sweep budget_floor = {floor:g}",
                (ax.get_xlim()[1], floor), color=MUTED, fontsize=9,
                va="bottom", ha="right", xytext=(-4, 3),
                textcoords="offset points")
    ax.set_xlabel("rank within warhead class"
                  + ("  (log scale)" if args.xscale == "log" else ""),
                  color=INK_2, fontsize=10)

    axes[0].legend(frameon=False, fontsize=9, labelcolor=INK_2, loc="upper right")

    sub = "  ·  ".join(f"{c}: {cross[c]} of "
                       f"{int((d.warhead_class == c).sum())} modes clear the floor"
                       for c in classes)
    fig.suptitle(f"Ranking score against rank position — {topic}",
                 color=INK, fontsize=13.5, x=0.045, ha="left", y=0.975)
    fig.text(0.045, 0.938, sub, color=MUTED, fontsize=9.5, ha="left")
    fig.text(0.045, 0.030,
             f"{len(d):,} ranked modes from {src}. Ranking is WITHIN class; the three are "
             f"not comparable across (#47).",
             color=MUTED, fontsize=8.5, ha="left")
    fig.text(0.045, 0.012,
             "Lower panel: each point is one mode, the line a 25-mode rolling median. "
             "Enrichment is NOT monotonic in rank.",
             color=MUTED, fontsize=8.5, ha="left")
    fig.subplots_adjust(top=0.895, bottom=0.105, left=0.075, right=0.965)

    t = sout.Topic("blacksmith", f"rank_curve_{topic}")
    dest = t.write(f"rank_curve_{args.xscale}", ".png")
    fig.savefig(dest, dpi=170, facecolor=SURFACE)
    print(f"\n  wrote {dest}")
    for c in classes:
        g = d[d.warhead_class == c].sort_values("class_rank")
        print(f"  {c:11s} n={len(g):4d}  top eb={g.conditional_eb.max():6.2f}  "
              f"{cross.get(c)} clear enrichment >= {floor:g}  "
              f"(deepest at rank {int(g[g.enrichment >= floor].class_rank.max()) if cross.get(c) else 0})")


if __name__ == "__main__":
    main()
