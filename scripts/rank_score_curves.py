#!/usr/bin/env python3
"""
Purpose: score-vs-rank curves for the support-weighted score, mode and molecule level.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-29
Input: rank_score_<topic>/rank_score_<N>.csv
Output: 00_outputs/blacksmith/rank_score_<topic>/rank_curves_<N>.png

@tt8804: "id need to look at the score vs rank curves" -- before choosing the
cutoff that `fraction_above` needs.

THE MOLECULE PANEL IS THE ONE THAT DECIDES ANYTHING. A mode-level curve shows
what the score does; the molecule-level curve at several candidate cutoffs shows
whether any of them SEPARATES molecules, which is what a selector has to do. If
every cutoff gives the same smooth decay, the cutoff is arbitrary and the
aggregation is sampling rather than selecting -- which is what the 1.2 ns sweep
already found for engagement within its own band.

BOTH SCORES ARE DRAWN ON THE MODE PANEL. Support moves the global ordering
almost not at all (rho = 0.9999) while reshuffling 30% of the top 450, and a
single curve cannot show both facts.
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

from shared import engagement_rank as er           # noqa: E402
from shared import outputs as sout                 # noqa: E402
from shared import run_paths as rp                 # noqa: E402
from shared import target_config as tc             # noqa: E402

log = logging.getLogger("rank-curves")
CUTS = [0.05, 0.20, 0.40, 0.60, 0.80]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=None)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    topic = a.topic or rp.topic()
    f = sorted(glob.glob(str(rp.BLACKSMITH / f"rank_score_{topic}/rank_score_*.csv")),
               key=os.path.getmtime)
    if not f:
        raise SystemExit("run scripts/apply_support_score.py first")
    d = pd.read_csv(f[-1])
    d["family"] = d.warhead_class.map(tc.family_of())
    log.info("%s: %d modes, %d molecules", topic, len(d), d.parent_ident.nunique())

    fig, ax = plt.subplots(2, 2, figsize=(14.6, 9.4))
    C = plt.cm.tab10

    # ---- A: mode-level, both scores ------------------------------------- #
    A = ax[0, 0]
    for i, (col, lab, st) in enumerate((("engagement", "engagement (geometry only)", "--"),
                                        ("rank_score", "rank_score (x support)", "-"))):
        y = np.sort(d[col].to_numpy())[::-1]
        A.plot(np.arange(1, len(y) + 1), y, st, lw=1.9, color=C(i),
               label=lab, alpha=.9)
    A.set(xscale="log", xlabel="mode rank (log)", ylabel="score",
          title=f"A · all {len(d):,} modes — support barely moves the curve")
    A.legend(fontsize=9); A.grid(alpha=.25, which="both")

    # ---- B: the reshuffle, where it actually happens ---------------------- #
    B = ax[0, 1]
    re = d.assign(r_e=d.engagement.rank(ascending=False, method="first"),
                  r_s=d.rank_score.rank(ascending=False, method="first"))
    top = re[re.r_e <= 2000]
    B.scatter(top.r_e, top.r_s - top.r_e, s=5, alpha=.35, color=C(3), lw=0)
    B.axhline(0, color="0.4", lw=1)
    B.set(xlabel="rank by geometry alone (top 2,000)",
          ylabel="rank change once support applies",
          title="B · it only reshuffles where geometry already ties")
    B.grid(alpha=.25)

    # ---- C: molecule level, the panel that decides ----------------------- #
    Cx = ax[1, 0]
    for i, c in enumerate(CUTS):
        g = er.rank_ligands(d, how="fraction_above", cutoff=c,
                            ligand_key="parent_ident", score_col="rank_score")
        y = g.ligand_engagement.to_numpy()
        Cx.plot(np.arange(1, len(y) + 1), y * 100, lw=1.8,
                color=C(i), label=f"cutoff {c:.2f}  ({int((y>0).sum())} mols > 0)")
    Cx.set(xlabel="molecule rank", ylabel="% of the molecule's modes above cutoff",
           title="C · molecule level — does any cutoff separate?")
    Cx.legend(fontsize=8.5); Cx.grid(alpha=.25)

    # ---- D: how many molecules survive each cutoff ----------------------- #
    D = ax[1, 1]
    xs = np.arange(0.02, 0.95, 0.02)
    nz, med = [], []
    for c in xs:
        g = er.rank_ligands(d, how="fraction_above", cutoff=float(c),
                            ligand_key="parent_ident", score_col="rank_score")
        nz.append(int((g.ligand_engagement > 0).sum()))
        med.append(float(g.ligand_engagement.median() * 100))
    D.plot(xs, nz, lw=2, color=C(4), label="molecules with ANY mode above cutoff")
    D.set(xlabel="cutoff on rank_score", ylabel="molecules (of 1,684)")
    D.grid(alpha=.25)
    D2 = D.twinx()
    D2.plot(xs, med, lw=2, ls="--", color=C(5))
    D2.set_ylabel("median % of a molecule's modes above cutoff", color=C(5))
    D2.tick_params(axis="y", colors=C(5))
    D.set_title("D · where a cutoff starts to bite")
    D.legend(fontsize=8.5, loc="lower left")

    fig.suptitle(f"{topic} — support-weighted rank score, mode and molecule level",
                 fontsize=13, y=.985)
    fig.tight_layout(rect=(0, 0, 1, .955))
    t = sout.Topic("blacksmith", f"rank_score_{topic}")
    out = t.write("rank_curves", ".png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(out)
    for c in CUTS:
        g = er.rank_ligands(d, how="fraction_above", cutoff=c,
                            ligand_key="parent_ident", score_col="rank_score")
        nzc = int((g.ligand_engagement > 0).sum())
        print(f"  cutoff {c:.2f}: {nzc:5d} molecules have any mode above it "
              f"({nzc/len(g)*100:4.0f}%)   median fraction "
              f"{g.ligand_engagement.median()*100:5.2f}%   top "
              f"{g.ligand_engagement.max()*100:5.1f}%")


if __name__ == "__main__":
    main()
