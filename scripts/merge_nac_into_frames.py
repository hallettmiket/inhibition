"""
Purpose: put the near-attack ranking into the frames the GUI reads, with its uncertainty attached.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: 00_outputs/blacksmith/nac_rank/*.csv + the latest T_3/T_4 frames
Output: new versioned D3_/D4_ parquet frames carrying nac_* columns

WHAT THIS ADDS, AND WHAT IT DELIBERATELY DOES NOT TOUCH.

Adds `nac_enrichment`, its 95% confidence interval, the run count that produced
it, and `shortlist_nac`. **It does not overwrite `shortlist` or
`shortlist_synth`.** Those are the existing PI-approved selections (issue #1), and
silently redefining what the GUI calls "the shortlist" would make two different
rankings indistinguishable on screen — the exact confusion this project keeps
paying for.

THE NUMBER COMES WITH ITS LIMITS ATTACHED, BECAUSE IT HAS SERIOUS ONES.

`nac_run_count` is written into every row on purpose. D0068 measured that
enrichment **does not converge**: the same molecules score 2.91x at 200 runs and
0.96x at 2,000, and the crystallographic positives — never selected on score —
fall identically. So an enrichment without its run count is not a quantity, and
the column exists so nobody can quote one.

`nac_enrichment_lo/hi` is a Wilson interval. At 200 runs it is ~1.12x wide
against a median of 1.59x, which is why the ranking is a **filter and not a fine
ordering**: 1,239 of 1,806 molecules have an interval reaching the top-25 band.
The GUI must show the interval, not the point estimate alone.

WHAT IT IS GOOD FOR, MEASURED. The top 300 selected at 200 runs sit at 0.96x on
the converged scale against 0.99x for known crystallographic binders and 0.76x
for random warhead-matched inactives — statistically indistinguishable from known
binders (p = 0.21), and clearly above random (AUC 0.620, p = 0.0017). It
concentrates active-like molecules. That is what a screen is for, and it is what
this column should be read as.

BDHI IS INCLUDED AND WAS RE-SCORED. D0067 found `sn2_ring_opening` applying sp3
backside geometry to an sp2 carbon, scoring 374 candidates at a median of 0.00x.
`load_scored` resolves duplicates newest-first by version, so the corrected rows
supersede.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import io as dio                     # noqa: E402
import nac_rank as nr                            # noqa: E402

log = logging.getLogger("merge-nac")

FRAMES = {"T_3": ("03_t3_reinvent", "D3"), "T_4": ("04_t4_combinatorial", "D4")}
DATA = Path("/data/lab_vm/append_only/inhibition")
SHORTLIST_N = 25          # matches the existing shortlist quota


def latest_frame(subdir: str, stem: str) -> Path:
    fs = list((DATA / subdir).glob(f"{stem}_*.parquet"))
    if not fs:
        raise SystemExit(f"no {stem} frames under {subdir}")
    return max(fs, key=lambda p: int(re.search(r"_(\d+)\.parquet$", p.name).group(1)))


def next_version(subdir: str, stem: str) -> int:
    fs = list((DATA / subdir).glob(f"{stem}_*.parquet"))
    return max(int(re.search(r"_(\d+)\.parquet$", p.name).group(1)) for p in fs) + 1


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    scored = nr.load_scored()
    scored = scored[scored.status == "ok"].copy()
    log.info("scored candidates available: %d", len(scored))

    null = (scored.viable_fraction / scored.enrichment).values
    lo, hi = nr.wilson_enrichment_ci(scored.viable_fraction.values,
                                     scored.n_poses.values.astype(float), null)
    scored["nac_enrichment_lo"], scored["nac_enrichment_hi"] = lo, hi
    scored = scored.rename(columns={"enrichment": "nac_enrichment",
                                    "viable_fraction": "nac_viable_fraction",
                                    "n_poses": "nac_run_count",
                                    "median_angle": "nac_median_angle",
                                    "median_dist": "nac_median_dist"})
    keep = ["ident", "nac_enrichment", "nac_enrichment_lo", "nac_enrichment_hi",
            "nac_viable_fraction", "nac_run_count", "nac_median_angle",
            "nac_median_dist", "mechanism"]
    scored = scored[[c for c in keep if c in scored.columns]]

    for approach, (subdir, stem) in FRAMES.items():
        src = latest_frame(subdir, stem)
        df = pd.read_parquet(src)
        df["_ident"] = df.candidate_id.astype(str)
        merged = df.merge(scored, left_on="_ident", right_on="ident",
                          how="left", suffixes=("", "_nac")).drop(
                              columns=["_ident", "ident"], errors="ignore")

        n_scored = merged.nac_enrichment.notna().sum()
        # Rank WITHIN the approach. Enrichment is comparable across mechanisms
        # only because each is divided by its own isotropic baseline, but the
        # existing shortlists are per-approach and this must line up with them.
        merged["shortlist_nac"] = False
        ranked = merged.dropna(subset=["nac_enrichment"]).nlargest(
            SHORTLIST_N, "nac_enrichment")
        merged.loc[ranked.index, "shortlist_nac"] = True

        log.info("%s: %s -> %d/%d scored, %d flagged shortlist_nac",
                 approach, src.name, n_scored, len(merged),
                 int(merged.shortlist_nac.sum()))
        if args.dry_run:
            continue
        dest = DATA / subdir / f"{stem}_{next_version(subdir, stem)}.parquet"
        merged.to_parquet(dest, index=False)
        log.info("  wrote %s", dest.name)


if __name__ == "__main__":
    main()
