"""
Purpose: T_2 step 2 — descriptors, novelty and structural alerts on the CReM frame.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the latest D2 frame (post generation)
Output: D2 with the shared annotation columns

Thin by design: the work is in `shared.annotate`, so T_2 and T_1 cannot drift
into computing "novelty" or "SAscore" differently from each other. The
integration phase pools these axes across all four approaches on one plot, and
that only means anything if one definition produced them.

ALERTS ARE ANNOTATED, NOT GATED. T_2 is the derivative neighbourhood of ATRA, a
conjugated polyene that trips BRENK for its own scaffold. Gating on alert count
would reject the seed's entire chemotype for having the feature that defines it —
the same trap T_4 avoids with R-group scoping, arriving here in a different
shape. The counts travel as weighable labels for the panel.
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
from shared import seeds as sd                     # noqa: E402

log = logging.getLogger("t2-annotate")

APPROACH = "t2"
DEFAULT_SEED = "atra"


def main() -> None:
    ap = argparse.ArgumentParser(description="T_2 step 2: annotate.")
    sd.add_seed_argument(ap, APPROACH)
    ap.add_argument("--alert-limit", type=int, default=None,
                    help="reject above this many alerts (default: annotate only)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    seed_name = args.seed or DEFAULT_SEED
    try:
        rec = sd.resolve(APPROACH, seed_name)
    except sd.SeedError as exc:
        raise SystemExit(str(exc)) from exc
    experiment = rec["experiment"]
    log.info("seed %s -> experiment %s", seed_name, experiment)

    df, frame_path = dio.latest_frame(experiment, APPROACH)
    log.info("loaded %s (%d rows)", frame_path.name, len(df))

    out_df = ann.annotate(df, approach=APPROACH, alert_limit=args.alert_limit)

    out = dio.write_full_frame(
        out_df, approach=APPROACH, experiment=experiment, stage="t2_annotate",
        params={"alert_scoping": "whole molecule",
                "alert_limit": args.alert_limit,
                "seed_name": seed_name,
                "novelty": "1 - max Tanimoto (ECFP4) vs the external set"},
        inputs={"d2_frame": frame_path})

    print(f"\nT_2 annotation (seed {seed_name}) -> {out}")
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
