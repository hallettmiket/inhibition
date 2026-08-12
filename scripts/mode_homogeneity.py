#!/usr/bin/env python3
"""
Purpose: measure whether a first-stage binding mode is homogeneous in the ONE
         quantity the criterion scores it on.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-12
Input: the run's aggregates and per-pose rows, via `rank_v2.load_v2`
Output: 00_outputs/blacksmith/mode_homogeneity/mode_homogeneity_<N>.csv
        (one row per first-stage cluster that was sub-split) + a summary line

THE QUESTION (#65). @tt8804, reading the ranking view: *"how are the sub modes so
diff from eachother, these should be individual modes not submodes."*

The first stage clusters on the reactive atom's POSITION and the DIRECTION the
warhead faces (`pose_modes.features`). The criterion scores a pose on the
reactive atom's DISTANCE to Cys113 SG and its approach ANGLE. Those are the same
two quantities. So poses grouped into one first-stage mode should score alike,
and sub-modes of one mode should differ only away from the reactive end -- which
is exactly what the second stage was introduced to separate (#61).

They do not. This measures by how much.

WHAT MAKES IT POSSIBLE. The first stage is DBSCAN at `eps = 3.0` on
`positional_separation + 2.0 * angular_difference_in_radians`. Two poses in the
same spot are neighbours up to 1.5 rad -- 86 degrees -- apart, and DBSCAN is a
CHAINING method: A joins B, B joins C, and A and C share a mode having never
been within eps of each other. A "mode" is a connected component, not a ball.

THE UNIT IS THE FIRST-STAGE CLUSTER, NOT THE MOLECULE. A molecule with two
genuinely distinct modes is not the failure being looked for; a single cluster
holding modes on opposite sides of the criterion window is. Only clusters that
produced more than one sub-mode can show it, so only those are emitted.

MEDIANS, NOT EXTREMES. A cluster's spread is taken between the sub-mode MEDIAN
distances, not between its furthest poses. One stray pose 8 A out would make
every cluster look heterogeneous; a sub-mode whose typical pose is 6.8 A out
while a sibling's is 3.6 A is the thing worth reporting, because those are the
values the score is actually computed from.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from shared import nac_criterion as nac              # noqa: E402
from shared import outputs as sout                   # noqa: E402
from shared import target_config as tc               # noqa: E402

log = logging.getLogger("mode-homogeneity")

B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")


def measure() -> tuple[pd.DataFrame, dict]:
    """One row per sub-split first-stage cluster, with its criterion spread."""
    import rank_v2 as rv                                          # noqa: PLC0415

    topic = str(tc.get("run.topic"))
    # Bound from `run.topic`, never left at rank_v2's own nac_v2 default --
    # D0080 is what that default costs when a script imports rather than runs it.
    rv.SRC = B / topic
    if not rv.SRC.is_dir():
        raise FileNotFoundError(f"no such topic directory: {rv.SRC}")
    agg, poses = rv.load_v2()
    if agg.empty or poses.empty:
        raise SystemExit(f"no aggregates or per-pose rows under {rv.SRC}")

    # The typical pose of each mode, in the criterion's own coordinates. Joined
    # on (molecule, mode) -- `rank_v2.load_v2` has already applied the
    # molecule-level supersession rule, so a re-screened molecule contributes
    # one run rather than two.
    per = (poses.groupby(["ident", "mode"])
                .agg(dist_med=("distance", "median"),
                     angle_med=("angle", "median"),
                     n=("mode", "size"))
                .reset_index()
                .rename(columns={"ident": "parent_ident"}))
    a = agg[agg["parent_mode"].notna()].merge(per, on=["parent_ident", "mode"])

    g = (a.groupby(["parent_ident", "parent_mode"])
          .agg(n_sub=("sub_index", "size"),
               warhead_class=("warhead_class", "first"),
               n_poses=("n", "sum"),
               dist_lo=("dist_med", "min"), dist_hi=("dist_med", "max"),
               angle_lo=("angle_med", "min"), angle_hi=("angle_med", "max"),
               enrich_lo=("enrichment", "min"), enrich_hi=("enrichment", "max"))
          .reset_index())
    g = g[g["n_sub"] > 1].copy()
    if g.empty:
        raise SystemExit("no first-stage cluster was sub-split; nothing to measure")

    lo, hi = nac.NAC_DIST_MIN, nac.NAC_DIST_MAX
    width = hi - lo
    g["dist_span"] = g["dist_hi"] - g["dist_lo"]
    g["angle_span"] = g["angle_hi"] - g["angle_lo"]
    g["enrich_span"] = g["enrich_hi"] - g["enrich_lo"]
    # WIDER THAN THE WINDOW ITSELF: the cluster's sub-modes disagree about the
    # scored quantity by more than the whole range the criterion accepts.
    g["spans_more_than_window"] = g["dist_span"] > width
    # STRADDLES: some sub-modes typically inside the window, some typically out.
    # This is the sharper statement, because it does not depend on the spread
    # being large -- only on the cluster spanning the decision boundary.
    g["straddles_window"] = (g["dist_lo"] <= hi) & (g["dist_hi"] > hi)

    summary = {
        "run_topic": topic,
        "window_a": [lo, hi],
        "n_clusters_subsplit": int(len(g)),
        "median_dist_span_a": round(float(g["dist_span"].median()), 3),
        "p90_dist_span_a": round(float(g["dist_span"].quantile(0.90)), 3),
        "max_dist_span_a": round(float(g["dist_span"].max()), 3),
        "n_spans_more_than_window": int(g["spans_more_than_window"].sum()),
        "frac_spans_more_than_window": round(
            float(g["spans_more_than_window"].mean()), 4),
        "n_straddles_window": int(g["straddles_window"].sum()),
        "frac_straddles_window": round(float(g["straddles_window"].mean()), 4),
    }
    return g.sort_values("dist_span", ascending=False), summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top", type=int, default=10,
                    help="worst clusters to print")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    g, s = measure()
    dest = sout.Topic("blacksmith", "mode_homogeneity").write(
        "mode_homogeneity", ".csv")
    g.to_csv(dest, index=False)

    lo, hi = s["window_a"]
    print("\n" + "=" * 76)
    print(f"  FIRST-STAGE MODE HOMOGENEITY — {s['run_topic']}")
    print("=" * 76)
    print(f"\n  the criterion accepts {lo}-{hi} A, a window {hi - lo:.1f} A wide")
    print(f"  first-stage clusters that were sub-split: "
          f"{s['n_clusters_subsplit']:,}\n")
    print("  spread of sub-mode MEDIAN distance-to-sulfur within one cluster:")
    print(f"    median {s['median_dist_span_a']} A   "
          f"90th pct {s['p90_dist_span_a']} A   max {s['max_dist_span_a']} A")
    print(f"\n  clusters spanning MORE than the whole window: "
          f"{s['n_spans_more_than_window']:,} "
          f"({100 * s['frac_spans_more_than_window']:.1f}%)")
    print(f"  clusters with sub-modes both INSIDE and OUTSIDE the window: "
          f"{s['n_straddles_window']:,} "
          f"({100 * s['frac_straddles_window']:.1f}%)")
    print(f"\n  worst {a.top}:")
    cols = ["parent_ident", "parent_mode", "n_sub", "warhead_class",
            "dist_lo", "dist_hi", "dist_span", "enrich_lo", "enrich_hi"]
    print(g.head(a.top)[cols].to_string(index=False))
    print(f"\n  written to {dest}\n")


if __name__ == "__main__":
    main()
