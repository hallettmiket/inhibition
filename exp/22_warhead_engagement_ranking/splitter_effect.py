#!/usr/bin/env python3
"""
Purpose: does contact grouping make a mode's engagement score mean anything?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-27
Input: two screens of the SAME molecules, one per splitting method
Output: 00_outputs/blacksmith/warhead_engagement/

run_all.py found that every MODE-LEVEL engagement statistic predicts the MD
outcome at rho ~ 0.11-0.13 while the single SIMULATED POSE predicts it at +0.652,
and offered an explanation: a shipped mode is a mixture. Its poses span a median
0.776 of the anchor-quality scale, which itself only runs 0 to 1, and 93% span
more than half of it. No summary of such a group can predict what one of its
members does.

THAT EXPLANATION MAKES A PREDICTION, WHICH IS WHY IT IS WORTH TESTING RATHER THAN
ASSERTING. If the aggregates fail because the groups are mixtures, then grouping
that produces TIGHT groups should shrink the spread -- and the aggregate should
become a usable score. If the spread does not fall, the explanation is wrong and
something else is limiting the mode-level statistics.

SAME MOLECULES, SAME POSE COUNT, ONE VARIABLE. Both screens run the same ten
idents at the same --nrun through the same receptor and the same seed; the only
difference is `--split-method`. Comparing spread across different molecules would
confound the rule with the chemistry.
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

from shared import engagement_rank as er           # noqa: E402
from shared import outputs as sout                 # noqa: E402
from shared import run_paths as rp                 # noqa: E402

log = logging.getLogger("splitter-effect")


def poses_for(topic: str) -> pd.DataFrame:
    fs = sorted(glob.glob(str(rp.BLACKSMITH / f"{topic}/poses_s*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit(f"no per-pose table for topic {topic!r}")
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return d.drop_duplicates(["ident", "pose_idx"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--contact-topic", default="cmp_contact_linkage")
    ap.add_argument("--dbscan-topic", default="cmp_warhead_dbscan")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    out = {}
    for label, topic in (("contact_linkage", a.contact_topic),
                         ("warhead_dbscan", a.dbscan_topic)):
        d = poses_for(topic)
        d = d[d["mode"] >= 0]
        e = er.mode_engagement(d)
        e["method"] = label
        out[label] = (d, e)
        log.info("%-16s %d molecules, %d poses, %d groups",
                 label, d.ident.nunique(), len(d), len(e))

    shared_ids = set(out["contact_linkage"][0].ident) & set(out["warhead_dbscan"][0].ident)
    P = print
    P("\n" + "=" * 84)
    P("  DOES CONTACT GROUPING MAKE A MODE'S ENGAGEMENT SCORE MEAN ANYTHING?")
    P("=" * 84)
    P(f"\n  {len(shared_ids)} molecules screened both ways, same runs, same seed\n")
    P(f"  {'method':<18}{'groups':>8}{'poses/group':>13}{'SPREAD median':>15}"
      f"{'>0.5 of scale':>15}")
    rows = []
    for label in ("warhead_dbscan", "contact_linkage"):
        d, e = out[label]
        e = e[e.ident.isin(shared_ids)]
        multi = e[e.n_poses_mode > 1]
        P(f"  {label:<18}{len(e):8d}{e.n_poses_mode.mean():13.1f}"
          f"{multi.engagement_spread.median():14.3f}"
          f"{(multi.engagement_spread > 0.5).mean() * 100:14.0f}%")
        rows.append(dict(method=label, groups=len(e),
                         poses_per_group=float(e.n_poses_mode.mean()),
                         spread_median=float(multi.engagement_spread.median()),
                         frac_over_half=float((multi.engagement_spread > 0.5).mean())))
    P("\n  SPREAD is the range of per-pose engagement inside one group, on a scale")
    P("  that runs 0 to 1. A group whose members span most of it has no summary")
    P("  worth ranking -- which is why the shipped rule's aggregates predict the")
    P("  MD outcome at rho ~ 0.11 while one real pose predicts it at +0.652.")

    a_ = out["warhead_dbscan"][1]
    b_ = out["contact_linkage"][1]
    am = a_[(a_.ident.isin(shared_ids)) & (a_.n_poses_mode > 1)].engagement_spread.median()
    bm = b_[(b_.ident.isin(shared_ids)) & (b_.n_poses_mode > 1)].engagement_spread.median()
    P(f"\n  VERDICT: spread {am:.3f} -> {bm:.3f}  "
      f"({'a ' + format(am / bm, '.1f') + 'x reduction' if bm and bm < am else 'NOT reduced'})")
    if bm >= am:
        P("  The explanation in run_all.py is NOT supported: tightening the groups")
        P("  did not make their engagement scores agree, so something else limits")
        P("  the mode-level statistics.")
    else:
        P("  Consistent with the explanation. It does NOT show the aggregate now")
        P("  predicts the outcome -- that needs modes swept under this splitter,")
        P("  which no MD run has been.")

    t = sout.Topic("blacksmith", "warhead_engagement")
    pd.DataFrame(rows).to_csv(t.write("splitter_effect", ".csv"), index=False)
    pd.concat([out[k][1] for k in out]).to_csv(
        t.write("splitter_effect_groups", ".csv"), index=False)
    P("\n" + "=" * 84)
    P(f"  written to {t.dir}")
    P("=" * 84 + "\n")


if __name__ == "__main__":
    main()
