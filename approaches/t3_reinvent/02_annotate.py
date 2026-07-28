"""
Purpose: T_3 step 2 — descriptors, novelty and structural alerts on the CReM frame.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the latest D2 frame (post generation)
Output: D2 with the shared annotation columns

Thin by design: the work is in `shared.annotate`, so T_3 and T_4 cannot drift
into computing "novelty" or "SAscore" differently from each other. The
integration phase pools these axes across all four approaches on one plot, and
that only means anything if one definition produced them.

ALERTS ARE R-GROUP SCOPED, exactly as in T_4 and for the same reason. Every T_3
candidate carries the acrylamide warhead by construction, and an acrylamide is a
Michael acceptor that BRENK flags on sight. Screening whole molecules would
reject the entire approach for having the feature it was designed around, so the
alert screen is scoped to the decoration outside the fixed scaffold.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import annotate as ann                 # noqa: E402
from shared import io as dio                       # noqa: E402

log = logging.getLogger("t3-annotate")

# The fixed scaffold, so alerts are scored on the decoration rather than on the
# warhead every candidate is required to have.
SCAFFOLD_SMARTS = "N(C(=O)C=C)C1CCS(=O)(=O)C1"

EXPERIMENT = "03_t3_reinvent"
APPROACH = "t3"


def main() -> None:
    ap = argparse.ArgumentParser(description="T_3 step 2: annotate.")
    ap.add_argument("--alert-limit", type=int, default=None,
                    help="reject above this many alerts (default: annotate only)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    df, frame_path = dio.latest_frame(EXPERIMENT, APPROACH)
    log.info("loaded %s (%d rows)", frame_path.name, len(df))

    out_df = ann.annotate(df, approach=APPROACH, core_smarts=SCAFFOLD_SMARTS,
                          alert_limit=args.alert_limit)

    out = dio.write_full_frame(
        out_df, approach=APPROACH, experiment=EXPERIMENT, stage="t3_annotate",
        params={"alert_scoping": "R-group outside the fixed scaffold",
                "alert_limit": args.alert_limit,
                "novelty": "1 - max Tanimoto (ECFP4) vs the external set"},
        inputs={"d3_frame": frame_path})

    print(f"\nT_3 annotation -> {out}")
    print(f"  {len(out_df)} candidates, "
          f"{int(out_df['rejected_at'].notna().sum())} stamped rejected and retained")
    for col, label in (("MW", "MW"), ("cLogP", "cLogP"), ("QED", "QED"),
                       ("SAscore", "SAscore"), ("novelty_external", "novelty")):
        if col in out_df:
            s = out_df[col].dropna()
            if len(s):
                print(f"    {label:9s} median {s.median():7.2f}   "
                      f"range {s.min():7.2f} to {s.max():7.2f}")
    if "whole_alert_total" in out_df:
        n = (out_df["whole_alert_total"].fillna(0) > 0).sum()
        print(f"    {int(n)}/{len(out_df)} carry at least one structural alert "
              "(annotated, not gated)")


if __name__ == "__main__":
    main()
