"""
Purpose: T_1 step 2 — descriptors, novelty and structural alerts on the CReM frame.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the latest D2 frame (post generation)
Output: D2 with the shared annotation columns

Thin by design: the work is in `shared.annotate`, so T_2 and T_1 cannot drift
into computing "novelty" or "SAscore" differently from each other. The
integration phase pools these axes across all four approaches on one plot, and
that only means anything if one definition produced them.

ALERTS ARE ANNOTATED, AND HERE THEY COULD REASONABLY BE GATED. T_1 has no seed
and no fixed core, so an alert is not an artifact of the starting material the
way it is in T_2, T_3 and T_4 — it is a property the model chose to generate.
The default is still annotate-only, because a de novo model's whole value is
proposing chemistry nobody would have picked, and filtering it against
medicinal-chemistry intuition before anything is measured discards exactly that.
Pass --alert-limit to gate.
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

log = logging.getLogger("t1-annotate")

EXPERIMENT = "01_t1_de_novo"
APPROACH = "t1"


def main() -> None:
    ap = argparse.ArgumentParser(description="T_1 step 2: annotate.")
    ap.add_argument("--alert-limit", type=int, default=None,
                    help="reject above this many alerts (default: annotate only)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    df, frame_path = dio.latest_frame(EXPERIMENT, APPROACH)
    log.info("loaded %s (%d rows)", frame_path.name, len(df))

    out_df = ann.annotate(df, approach=APPROACH, alert_limit=args.alert_limit)

    out = dio.write_full_frame(
        out_df, approach=APPROACH, experiment=EXPERIMENT, stage="t1_annotate",
        params={"alert_scoping": "whole molecule",
                "alert_limit": args.alert_limit,
                "novelty": "1 - max Tanimoto (ECFP4) vs the external set"},
        inputs={"d1_frame": frame_path})

    print(f"\nT_1 annotation -> {out}")
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
