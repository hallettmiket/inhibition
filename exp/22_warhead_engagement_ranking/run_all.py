#!/usr/bin/env python3
"""
Purpose: which static per-mode metric predicts warhead engagement in MD?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-27
Input: the 147 swept modes with an ok result, the rank_v2 table, the per-pose table
Output: 00_outputs/blacksmith/warhead_engagement/

@tt8804: "now we need to fix ranking accordingly, focusing on warhead optimal
engagement only for now."

THE PROJECT HAS NEVER TESTED ITS RANKING AGAINST AN OUTCOME. `rank_validated` is
False on every shortlist because the enrichment gate fired (D0041), but that gate
asks whether docking separates actives from inactives. This asks something the
project can actually answer with data on disk: does a mode's STATIC score predict
how much of an MD trajectory it spends in attack geometry?

THE TARGET IS `frac_attack_ready` -- the fraction of swept frames in which the
warhead is within the near-attack window. It is the dynamic measurement of
exactly the thing the static metrics estimate, so this is a like-for-like test
rather than a proxy.

THE SAMPLE IS RANGE-RESTRICTED AND THAT IS NOT FIXABLE HERE. These 147 modes were
SELECTED for sweeping by `conditional_eb`, the incumbent metric. So every
correlation below is measured inside the band the incumbent already liked, which
flatters nothing in particular but does mean a metric can look weak here and be
strong across the full 4,432. The comparison between metrics is still valid --
they all face the same restriction -- but no absolute number here is a population
estimate. Stated at the top rather than in a footnote because it is the main
threat to every conclusion in this file.

WHY NOT RANK ON `anchor_quality_max`. It is the obvious "best engagement" choice
and it is measured to be the WORST selector: picking a mode's representative by
argmax anchoring recovers the crystal pose in 6.7% of 15 complexes, against 33.3%
for the medoid of the well-anchored quartile (nac_screen_v2). The maximum of a
noisy score is an outlier, and it grows with mode size for free. Candidates here
are therefore central statistics as well as extremes, and mode size is carried
alongside so a metric that is secretly a size proxy is visible.
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

from shared import nac_criterion as nac            # noqa: E402
from shared import outputs as sout                 # noqa: E402
from shared import run_paths as rp                 # noqa: E402

log = logging.getLogger("engagement")


def newest(pattern: str) -> Path:
    fs = sorted(glob.glob(str(rp.BLACKSMITH / pattern)), key=os.path.getmtime)
    if not fs:
        raise SystemExit(f"nothing matches {pattern}")
    return Path(fs[-1])


def sweep_outcomes() -> pd.DataFrame:
    fs = sorted(glob.glob(str(rp.BLACKSMITH / f"attack_sweep_{rp.topic()}/attack_sweep_*.csv")),
                key=os.path.getmtime)
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d.drop_duplicates("ident", keep="last")
    return d[d.status.astype(str).str.startswith("ok")].copy()


def per_pose() -> pd.DataFrame:
    fs = sorted(glob.glob(str(rp.BLACKSMITH / f"{rp.topic()}/poses_s*.csv")),
                key=os.path.getmtime)
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return d.drop_duplicates(["ident", "pose_idx"])


def engagement_stats(g: pd.DataFrame) -> dict:
    """Candidate per-mode summaries of one group's per-pose anchor quality.

    Several statistics of ONE quantity, because "how well does this group engage"
    has no single obvious summary and the obvious one (max) is measured to be the
    worst. `q75_mean` is the well-anchored-quartile rule that beat argmax 33.3%
    to 6.7% when choosing a representative; it is included here to find out
    whether what selects a good POSE also orders a good MODE.
    """
    a = g["anchor"].to_numpy(dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {}
    top = a[a >= np.percentile(a, 75)] if a.size >= 4 else a
    return {
        "eng_max": float(a.max()),
        "eng_mean": float(a.mean()),
        "eng_median": float(np.median(a)),
        "eng_p90": float(np.percentile(a, 90)),
        "eng_q75_mean": float(top.mean()),
        # frequency of clearing the binary window, for contrast with the
        # continuous statistics above
        "eng_frac_viable": float(g["viable"].mean()),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--top-k", type=int, default=20)
    ap.add_argument("--ready-floor", type=float, default=0.01,
                    help="frac_attack_ready above which a mode counts as engaging")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from scipy.stats import spearmanr

    sw = sweep_outcomes()
    rk = pd.read_csv(newest(f"rank_v2/rank_v2_T4_{rp.topic()}_conditional_eb_*.csv"))
    log.info("%d swept modes, %d ranked modes", len(sw), len(rk))

    # per-pose anchor quality, aggregated per mode
    pp = per_pose()
    pp = pp[pp["mode"] >= 0].copy()
    pp["anchor"] = [nac.anchor_quality(d_, an, me)
                    for d_, an, me in zip(pp.distance, pp.angle, pp.mechanism)]
    eng = (pp.groupby(["ident", "mode"], group_keys=False)
             .apply(lambda g: pd.Series(engagement_stats(g)))
             .reset_index()
             .rename(columns={"ident": "parent_ident"}))
    log.info("engagement statistics for %d modes", len(eng))

    j = sw.merge(rk, on=["parent_ident", "mode"], how="inner", suffixes=("_sw", ""))
    j = j.merge(eng, on=["parent_ident", "mode"], how="left")
    log.info("joined: %d modes carry both an outcome and a score", len(j))

    CANDIDATES = [
        ("conditional_eb", "the INCUMBENT ranking column"),
        ("enrichment", "viable_fraction / isotropic_null"),
        ("viable_fraction", "fraction of poses clearing the window"),
        ("eng_frac_viable", "same, recomputed per pose"),
        ("anchor_quality_max", "max anchor quality (the tempting one)"),
        ("anchor_quality_mean", "mean anchor quality"),
        ("anchor_quality_p90", "90th percentile anchor quality"),
        ("eng_max", "max, recomputed"),
        ("eng_mean", "mean, recomputed"),
        ("eng_median", "median, recomputed"),
        ("eng_p90", "p90, recomputed"),
        ("eng_q75_mean", "mean of the well-anchored quartile"),
        ("start_attack_ready", "was the SWEPT POSE itself attack-ready"),
        ("mean_energy", "AutoDock energy (lower is better)"),
        ("n_poses_mode", "mode size -- the size-proxy control"),
    ]
    y = j["frac_attack_ready"].to_numpy(dtype=float)
    hit = y > a.ready_floor

    rows = []
    for col, note in CANDIDATES:
        if col not in j.columns:
            continue
        x = pd.to_numeric(j[col], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 20:
            continue
        sign = -1.0 if col == "mean_energy" else 1.0
        rho, p = spearmanr(sign * x[ok], y[ok])
        order = np.argsort(-(sign * x[ok]))
        k = min(a.top_k, ok.sum())
        prec = float(hit[ok][order[:k]].mean())
        rows.append(dict(metric=col, note=note, n=int(ok.sum()), rho=float(rho),
                         p=float(p), top_k=k, top_k_precision=prec,
                         size_rho=float(spearmanr(sign * x[ok],
                                                  j["n_poses_mode"].to_numpy()[ok])[0])))
    res = pd.DataFrame(rows).sort_values("rho", ascending=False)

    t = sout.Topic("blacksmith", "warhead_engagement")
    res.to_csv(t.write("metric_comparison", ".csv"), index=False)
    j.to_csv(t.write("joined", ".csv"), index=False)

    base = float(hit.mean())
    P = print
    P("\n" + "=" * 92)
    P("  WHICH STATIC METRIC PREDICTS WARHEAD ENGAGEMENT IN MD?")
    P("=" * 92)
    P(f"\n  {len(j)} swept modes · target = frac_attack_ready > {a.ready_floor} "
      f"· base rate {base * 100:.0f}%\n")
    P(f"  {'metric':<22}{'rho':>8}{'p':>10}{'top-20 hit':>12}{'lift':>7}"
      f"{'rho vs size':>13}   note")
    for _, r in res.iterrows():
        P(f"  {r.metric:<22}{r.rho:+8.3f}{r.p:10.3g}{r.top_k_precision * 100:11.0f}%"
          f"{r.top_k_precision / base if base else float('nan'):6.2f}x"
          f"{r.size_rho:+13.3f}   {r.note}")
    P(f"\n  'lift' is top-20 hit rate over the {base * 100:.0f}% base rate. 1.00x is no better")
    P("  than picking at random from the swept set.")
    P("  'rho vs size' exposes a metric that is secretly ranking by mode size.")
    P("\n  RANGE RESTRICTION: these 147 were SELECTED by conditional_eb, so every")
    P("  number here is measured inside the band the incumbent already liked.")
    P("  Metric-vs-metric comparison is fair; no absolute value is a population estimate.")
    P("\n" + "=" * 92)
    P(f"  written to {t.dir}")
    P("=" * 92 + "\n")


if __name__ == "__main__":
    main()
