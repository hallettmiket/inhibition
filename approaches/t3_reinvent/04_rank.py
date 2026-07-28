"""
Purpose: T_3 (generative decoration, LibInvent) — rank the docked candidates and build the shortlist.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the latest D3 frame (post docking)
Output: D3 with rank / percentile / ligand-efficiency columns and a shortlist flag

The ranking itself lives in `shared.rank_shortlist`, which all four approaches
call. Each supplies only its own rank metric and grouping; everything else — tie
handling, what a percentile means, what a quota means, and what gets attached to
a ranking before anyone reads it — is identical across approaches, because the
integration phase pools them and four private definitions of "best" is exactly
what that phase cannot absorb.

RANK METRIC: `affinity_kcal` (kcal/mol, LOWER is better).

ONE GROUP. single fixed warhead, so within-class ranking is vacuous and this is a straight ranking.

THIS RANKING IS NOT VALIDATED (D0031). The covalent enrichment gate on this
receptor is indistinguishable from chance, so the shortlist is an ordering the
pipeline produced and not evidence that the molecules at the top bind. The gate
verdict is attached to every ranked row and must be displayed beside the rank.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                     # noqa: E402
from shared import rank_shortlist as rs          # noqa: E402

log = logging.getLogger("t3-rank")

EXPERIMENT = "03_t3_reinvent"
RANK_METRIC = "affinity_kcal"
STRATUM = "covalent"


def main() -> None:
    ap = argparse.ArgumentParser(description="T_3 (generative decoration, LibInvent): rank and shortlist.")
    ap.add_argument("--quota", type=int, default=25,
                    help="how many candidates the shortlist carries forward")
    ap.add_argument("--min-docked", type=int, default=20,
                    help="below this a group's rank is flagged not selective")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    df, frame_path = dio.latest_frame(EXPERIMENT, "t3")
    log.info("loaded %s (%d rows)", frame_path.name, len(df))
    if RANK_METRIC not in df.columns:
        raise SystemExit(f"frame has no {RANK_METRIC!r} — run the docking stage first")

    gated = rs.attach_gate(df, STRATUM, RANK_METRIC)
    ranked = rs.rank(gated, metric=RANK_METRIC, group_col=None,
                     min_docked=args.min_docked)
    final = rs.shortlist(ranked, quota=args.quota)

    n_short = int(final["shortlist"].sum())
    out = dio.write_full_frame(
        final, approach="t3", experiment=EXPERIMENT, stage="t3_rank",
        params={"rank_metric": f"{RANK_METRIC} (lower better)",
                "ranking_scope": "whole approach (single group)",
                "quota": args.quota,
                "min_docked_for_meaningful_rank": args.min_docked,
                "gate_verdict": str(final["gate_verdict"].iloc[0]),
                "rank_validated": bool(final["rank_validated"].iloc[0]),
                "n_shortlisted": n_short},
        inputs={"frame": frame_path})

    print(f"\nT_3 (generative decoration, LibInvent) ranking -> {out}")
    print(f"  shortlisted {n_short} candidates\n")
    print(rs.summarise(final, RANK_METRIC))
    print("\n  shortlist:")
    s = final[final["shortlist"]].sort_values("rank")
    for _, r in s.head(15).iterrows():
        print(f"    #{int(r['rank']):3d}  {r[RANK_METRIC]:7.2f} kcal/mol  "
              f"LE {r['ligand_efficiency'] if r['ligand_efficiency'] is not None else float('nan')}  "
              f"{r['candidate_id']}")
    if n_short > 15:
        print(f"    ... and {n_short - 15} more")


if __name__ == "__main__":
    main()
