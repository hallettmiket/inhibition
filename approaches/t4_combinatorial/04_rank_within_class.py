"""
Purpose: T_4 step 8 — rank docked candidates WITHIN warhead class and build the shortlist.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the latest D4 frame (post covalent docking)
Output: D4 with within-class rank columns and a quota-based shortlist flag

WHY WITHIN CLASS. The rank metric is gnina's Vina-style affinity (D0015), and
that number is comparable only among molecules docked the same way. Each warhead
class is docked against a DIFFERENT reactive-atom constraint, so their pose
ensembles are not samples from one distribution; on top of that, warheads differ
in heavy-atom count and a Vina-style score tracks size. Sorting all 1,683
survivors together would mostly rediscover which warhead is biggest and
greasiest, and would hand the entire shortlist to one chemotype.

So each class is ranked against itself and contributes a fixed quota. The
shortlist is then a designed comparison ACROSS chemotypes — which is the
question T_4 exists to answer, since the core is fixed and the warhead is the
variable.

LIGAND EFFICIENCY IS ADVISORY. LE = -affinity / heavy_atoms is reported because
it is the standard size correction and reviewers will ask. It is NOT the rank
metric: D0015 fixed that on affinity after the enrichment gate measured it, and
no equivalent measurement exists for LE on this target. Swapping the rank metric
for an unmeasured one because it is conventional would undo the only calibration
the choreography actually has.

FLAGS TRAVEL WITH THE ROWS. Classes stamped OUTSIDE_WINDOW by the reactivity
triage (D0019) are ranked and do reach the shortlist, carrying the flag. Classes
with too few successful docks are ranked and flagged too. Nothing is dropped
here — this stage stamps and orders, it does not reject.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                      # noqa: E402

log = logging.getLogger("t4-rank")

EXPERIMENT = "04_t4_combinatorial"
CONFIG = REPO / "config" / "approaches" / "t4_combinatorial.yaml"

RANK_METRIC = "affinity_kcal"      # LOWER is better (D0015)


def load_ranking_config() -> dict:
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    r = cfg.get("ranking") or {}
    return {
        "per_class_quota": int(r.get("per_class_quota", 3)),
        "min_docked": int(r.get("min_docked_for_meaningful_rank", 20)),
        "include_flagged": bool(r.get("include_flagged_classes", True)),
    }


def rank_within_class(df: pd.DataFrame, min_docked: int) -> pd.DataFrame:
    """Add within-class rank, class size, percentile and ligand efficiency.

    Rows without a docking result keep null ranks rather than sorting to the
    bottom — "did not dock" and "docked badly" are different facts and the
    shortlist must not confuse them.
    """
    out = df.copy()
    for col in ("class_rank", "class_n_docked", "class_percentile",
                "ligand_efficiency", "rank_is_selective"):
        out[col] = pd.NA

    docked_mask = out[RANK_METRIC].notna()
    log.info("%d of %d rows carry a docking result", int(docked_mask.sum()), len(out))

    for cls, idx in out[docked_mask].groupby("warhead_class").groups.items():
        sub = out.loc[idx].sort_values(RANK_METRIC, ascending=True)   # lower better
        n = len(sub)
        ranks = range(1, n + 1)
        out.loc[sub.index, "class_rank"] = list(ranks)
        out.loc[sub.index, "class_n_docked"] = n
        # Percentile of the class that this candidate beats; 100 = best in class.
        out.loc[sub.index, "class_percentile"] = [
            round(100.0 * (n - r) / (n - 1), 2) if n > 1 else 100.0 for r in ranks
        ]
        out.loc[sub.index, "rank_is_selective"] = bool(n >= min_docked)
        if n < min_docked:
            log.warning("class %r has only %d successful docks (< %d): its "
                        "within-class rank is not selective and is flagged",
                        cls, n, min_docked)

    # Ligand efficiency, advisory only. HAC comes from the descriptor stage.
    if "HAC" in out.columns:
        hac = pd.to_numeric(out["HAC"], errors="coerce")
        aff = pd.to_numeric(out[RANK_METRIC], errors="coerce")
        out["ligand_efficiency"] = (-aff / hac).where(hac > 0).round(4)
    else:
        log.warning("no HAC column; ligand efficiency not computed")
    return out


def build_shortlist(df: pd.DataFrame, quota: int, include_flagged: bool) -> pd.DataFrame:
    """Flag the top `quota` of each warhead class, recording why each was taken."""
    out = df.copy()
    out["shortlist"] = False
    out["shortlist_reason"] = pd.NA

    ranked = out[out["class_rank"].notna()]
    for cls, sub in ranked.groupby("warhead_class"):
        flag = str(sub["reactivity_flag"].iloc[0]) if "reactivity_flag" in sub else ""
        if flag == "OUTSIDE_WINDOW" and not include_flagged:
            log.info("class %r skipped: OUTSIDE_WINDOW and include_flagged=False", cls)
            continue
        take = sub.sort_values("class_rank").head(quota)
        out.loc[take.index, "shortlist"] = True
        notes = []
        if flag and flag != "IN_WINDOW":
            notes.append(f"reactivity={flag}")
        if not bool(sub["rank_is_selective"].iloc[0]):
            notes.append(f"only {int(sub['class_n_docked'].iloc[0])} docked — "
                         "rank not selective")
        reason = f"top-{quota} of {cls}"
        if notes:
            reason += " (" + "; ".join(notes) + ")"
        out.loc[take.index, "shortlist_reason"] = reason
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="T_4 step 8: rank within warhead class.")
    ap.add_argument("--quota", type=int, default=None,
                    help="override per-class quota from config")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    cfg = load_ranking_config()
    quota = args.quota or cfg["per_class_quota"]

    df, frame_path = dio.latest_frame(EXPERIMENT, "t4")
    log.info("loaded %s (%d rows)", frame_path.name, len(df))
    if RANK_METRIC not in df.columns:
        raise SystemExit(
            f"frame has no {RANK_METRIC!r} column — run 03_covalent_dock.py first")
    if not df[RANK_METRIC].notna().any():
        raise SystemExit(
            f"no row carries a {RANK_METRIC} value; docking produced nothing to rank")

    ranked = rank_within_class(df, cfg["min_docked"])
    final = build_shortlist(ranked, quota, cfg["include_flagged"])

    n_short = int(final["shortlist"].sum())
    n_classes = final.loc[final["shortlist"], "warhead_class"].nunique()

    out = dio.write_full_frame(
        final, approach="t4", experiment=EXPERIMENT, stage="t4_rank_within_class",
        params={"rank_metric": f"{RANK_METRIC} (lower better, D0015)",
                "ranking_scope": "within warhead class",
                "per_class_quota": quota,
                "min_docked_for_meaningful_rank": cfg["min_docked"],
                "include_flagged_classes": cfg["include_flagged"],
                "n_shortlisted": n_short,
                "n_classes_represented": int(n_classes)},
        inputs={"d4_frame": frame_path})

    print(f"\nT_4 within-class ranking -> {out}")
    print(f"  shortlisted {n_short} candidates across {n_classes} warhead classes\n")

    hdr = (f"  {'warhead class':24s} {'docked':>6s} {'best':>7s} {'median':>7s} "
           f"{'best LE':>7s}  flags")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    dk = final[final[RANK_METRIC].notna()]
    for cls, g in dk.groupby("warhead_class"):
        flag = str(g["reactivity_flag"].iloc[0]) if "reactivity_flag" in g else ""
        marks = []
        if flag and flag != "IN_WINDOW":
            marks.append(flag)
        if not bool(g["rank_is_selective"].iloc[0]):
            marks.append("RANK_NOT_SELECTIVE")
        le = pd.to_numeric(g["ligand_efficiency"], errors="coerce").max()
        print(f"  {cls:24s} {len(g):6d} {g[RANK_METRIC].min():7.2f} "
              f"{g[RANK_METRIC].median():7.2f} "
              f"{le:7.3f}  {','.join(marks) if marks else '-'}")

    print("\n  shortlist (ranked within class; affinity NOT comparable across classes):")
    s = final[final["shortlist"]].sort_values(["warhead_class", "class_rank"])
    for _, r in s.iterrows():
        print(f"    {r['warhead_class']:22s} #{int(r['class_rank'])}  "
              f"{r[RANK_METRIC]:6.2f} kcal/mol  {r['candidate_id']}")


if __name__ == "__main__":
    main()
