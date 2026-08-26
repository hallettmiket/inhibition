#!/usr/bin/env python3
"""
Purpose: the ranking curve with the 12-pose gate, and without it, on bdhi.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: every bdhi molecule's persisted cloud + its per-pose table (this run)
Output: 00_outputs/blacksmith/gate_vs_all_<topic>/  (csv + png)

@tt8804: "show me side by side the same ranking score to ranking place graph
comparing on bdhi only if we only consider consensus cutoff modes or all poses
with modes as a decluttering function. how many more modes are we looking at."

THE ONE THING THAT CHANGES BETWEEN PANELS IS THE GATE. Both sides use the SAME
HDBSCAN grouping, the same scores, the same code path. The left keeps only
groups of >= 12 poses (`ranking.mode_gate`, D0084); the right keeps every group
INCLUDING SINGLETONS, which is the de-duplication reading -- a pose HDBSCAN
calls noise is a group of one, not a pose that failed to exist.

Isolating the gate matters because the alternative comparison -- today's shipped
pipeline against HDBSCAN-ungated -- moves two things at once and could not say
which one produced the difference.

SCORES ARE RECOMPUTED, NOT REUSED. `conditional_eb` for these groups cannot come
off the rank table: that table's modes are the SHIPPED clustering's. Every
ingredient is rebuilt per group from the per-pose verdicts and pushed through
`rank_v2.derived_scores`, so the y axis means the same thing it does in
production.
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

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import nac_criterion as nac              # noqa: E402
from shared import outputs as sout                   # noqa: E402
from shared import pose_cluster as pc                # noqa: E402
from shared import run_paths as rp                   # noqa: E402
from shared import target_config as tc               # noqa: E402
import rank_v2 as rv                                 # noqa: E402

log = logging.getLogger("gate-vs-all")

SERIES = ["#eb6834", "#1baf7a"]          # palette slots 2,3 — bdhi_c4, bdhi_c5
INK, INK_2, MUTED = "#0b0b0b", "#3a3a3a", "#6b6b6b"
SURFACE, GRID, RULE = "#fcfcfb", "#e8e8e6", "#cfcfcc"


def groups_for(ident: str, poses: pd.DataFrame) -> pd.DataFrame | None:
    """One row per HDBSCAN group of this molecule's cloud, singletons included."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    sdf = rp.allposes_dir() / f"{ident}.sdf"
    ms = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=False)
          if m is not None]
    if not ms:
        return None
    t = poses[poses.ident == ident]
    t = t[t["mode"] != -1].sort_values(["mode", "pose_idx"]).reset_index(drop=True)
    if len(t) != len(ms):
        return None                      # refuse to pair; do not guess (§1.6a)
    heavy = [i for i, a in enumerate(ms[0].GetAtoms()) if a.GetAtomicNum() > 1]
    xyz = np.array([m.GetConformer().GetPositions()[heavy] for m in ms])
    lab = pc.cluster(xyz)
    # NOISE IS A SINGLETON: give each noise pose its own group id.
    lab = lab.copy()
    nxt = (lab.max() + 1) if (lab >= 0).any() else 0
    for i in np.flatnonzero(lab == -1):
        lab[i] = nxt
        nxt += 1
    v = t["viable"].to_numpy(bool)
    ir = t["in_range"].to_numpy(bool)
    null = float(nac.isotropic_null(str(t["mechanism"].iloc[0])))
    out = []
    for k in sorted(set(lab)):
        s = lab == k
        n = int(s.sum())
        out.append(dict(
            parent_ident=ident, group=int(k), warhead_class=None,
            n_poses=len(lab), n_poses_mode=n, consensus=n / len(lab),
            n_in_range=int((ir & s).sum()),
            n_viable=int((v & s).sum()),
            n_viable_given_in_range=int((v & ir & s).sum()),
            viable_fraction=float(v[s].mean()),
            isotropic_null=null,
            is_singleton=(n == 1)))
    return pd.DataFrame(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    topic = tc.get("run.topic")
    rv.SRC = rp.BLACKSMITH / topic
    agg, poses = rv.load_v2()
    cls = agg.drop_duplicates("parent_ident").set_index("parent_ident")["warhead_class"]
    bdhi = sorted(i for i, c in cls.items() if str(c).startswith("bdhi"))
    if a.limit:
        bdhi = bdhi[:a.limit]
    log.info("bdhi molecules: %d", len(bdhi))

    frames, skipped = [], 0
    for n, ident in enumerate(bdhi, 1):
        g = groups_for(ident, poses)
        if g is None:
            skipped += 1
            continue
        g["warhead_class"] = str(cls[ident])
        frames.append(g)
        if n % 50 == 0:
            log.info("  %d/%d  (%d skipped)", n, len(bdhi), skipped)
    d = pd.concat(frames, ignore_index=True)
    log.info("groups: %d over %d molecules; %d molecules skipped "
             "(cloud/table length mismatch)", len(d), d.parent_ident.nunique(), skipped)

    d["enrichment"] = d.viable_fraction / d.isotropic_null
    d["enrichment_conditional"] = d.apply(rv.conditional_enrichment, axis=1)
    d = rv.derived_scores(d)

    gate = int(tc.rank_min_mode_poses())
    views = {f"gated  (>= {gate} poses)": d[d.n_poses_mode >= gate].copy(),
             "all groups (singletons kept)": d.copy()}
    for name, v in views.items():
        v["class_rank"] = v.groupby("warhead_class")["conditional_eb"] \
                           .rank(ascending=False, method="min")

    t = sout.Topic("blacksmith", f"gate_vs_all_{topic}")
    d.to_csv(t.write("groups", ".csv"), index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    classes = sorted(d.warhead_class.unique())
    colour = {c: SERIES[i] for i, c in enumerate(classes)}
    floor = float(tc.sweep_budget_floor())

    # TWO ROWS. The top is the whole curve; the bottom is the SAME data clipped
    # to y >= floor, which is the only part the sweep ever sees (@tt8804). The
    # bottom row's x is clipped too -- on the full rank axis the survivors are a
    # sliver at the left edge and the shape that matters is unreadable.
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 9.6), facecolor=SURFACE)
    for row in (0, 1):
        for col, (name, v) in enumerate(views.items()):
            ax = axes[row][col]
            ax.set_facecolor(SURFACE)
            ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"):
                ax.spines[sp].set_color(RULE)
            ax.tick_params(colors=INK_2, labelsize=9)
            xmax = 0
            for c in classes:
                g = v[(v.warhead_class == c) & v.conditional_eb.notna()] \
                    .sort_values("class_rank")
                ax.plot(g.class_rank, g.conditional_eb, color=colour[c],
                        linewidth=2.0, solid_capstyle="round", zorder=3, label=c)
                over = g[g.conditional_eb >= floor]
                if len(over):
                    xmax = max(xmax, float(over.class_rank.max()))
            ax.axhline(floor, color=MUTED, linewidth=1.3,
                       linestyle=(0, (5, 3)), zorder=2)
            if row == 0:
                n_over = int((v.conditional_eb >= floor).sum())
                ax.set_title(f"{name}\n{len(v):,} modes   ·   "
                             f"{n_over:,} score >= {floor:g}",
                             color=INK, fontsize=11, loc="left", pad=8)
            else:
                ax.set_ylim(floor, None)
                ax.set_xlim(0, (xmax * 1.06) if xmax else 1)
                per = {c: int(((v.warhead_class == c)
                               & (v.conditional_eb >= floor)).sum())
                       for c in classes}
                ax.set_title("above the floor only   ·   "
                             + ",  ".join(f"{c} {n:,}" for c, n in per.items()),
                             color=INK, fontsize=10.5, loc="left", pad=8)
                ax.set_xlabel("rank within warhead class", color=INK_2, fontsize=10)
            if col == 0:
                ax.set_ylabel("conditional_eb", color=INK_2, fontsize=10)
    axes[0][0].legend(frameon=False, fontsize=9, labelcolor=INK_2, loc="upper right")
    axes[0][1].annotate(f"budget_floor = {floor:g}",
                        (axes[0][1].get_xlim()[1], floor), color=MUTED,
                        fontsize=9, va="bottom", ha="right",
                        xytext=(-4, 3), textcoords="offset points")

    grow = len(views["all groups (singletons kept)"]) / max(1, len(list(views.values())[0]))
    fig.suptitle("The 12-pose gate, and without it — bdhi only", color=INK,
                 fontsize=13.5, x=0.035, ha="left", y=0.985)
    n_single_over = int(d[d.is_singleton].conditional_eb.ge(floor).sum())
    fig.text(0.035, 0.949,
             f"Same HDBSCAN grouping and the same scores on both sides; only the gate "
             f"differs.  Dropping it multiplies the candidate count by {grow:.1f}× "
             f"and adds {n_single_over} singletons above the floor.",
             color=MUTED, fontsize=9.5, ha="left")
    sv = d[d.is_singleton].conditional_eb.dropna()
    fig.text(0.035, 0.030,
             f"{d.parent_ident.nunique()} bdhi molecules, {skipped} skipped for a "
             f"cloud/table length mismatch. Singletons are "
             f"{d.is_singleton.mean()*100:.0f}% of all groups.",
             color=MUTED, fontsize=8.5, ha="left")
    fig.text(0.035, 0.011,
             f"The flat step on the right is the score's floor of resolution: a "
             f"singleton has 0 or 1 poses in the window, so conditional_eb takes only "
             f"{sv.nunique()} distinct values across {len(sv):,} scorable singletons.",
             color=MUTED, fontsize=8.5, ha="left")
    fig.subplots_adjust(top=0.895, bottom=0.105, left=0.058, right=0.985,
                        wspace=0.13, hspace=0.30)
    dest = t.write("gate_vs_all", ".png")
    fig.savefig(dest, dpi=170, facecolor=SURFACE)

    print("\n" + "=" * 74)
    print("  THE 12-POSE GATE, AND WITHOUT IT — bdhi")
    print("=" * 74)
    for name, v in views.items():
        print(f"\n  {name}")
        print(f"    modes                 {len(v):,}")
        print(f"    scorable (eb notna)   {int(v.conditional_eb.notna().sum()):,}")
        print(f"    score >= {floor:g}          {int((v.conditional_eb >= floor).sum()):,}")
        for c in classes:
            print(f"      {c:9s} {int((v.warhead_class == c).sum()):,}")
    print(f"\n  singletons: {int(d.is_singleton.sum()):,} of {len(d):,} groups "
          f"({d.is_singleton.mean()*100:.0f}%)")
    print(f"  today's shipped ranking carries 1,432 ranked bdhi modes for comparison")
    print(f"\n  wrote {dest}\n")


if __name__ == "__main__":
    main()
