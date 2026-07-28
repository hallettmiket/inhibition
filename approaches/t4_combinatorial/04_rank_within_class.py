"""
Purpose: T_4 step 8 — rank docked candidates within ADDUCT class and build the shortlist.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27 (moved to the shared ranker, regrouped on adduct class 2026-07-28)
Input: the latest D4 frame (post covalent docking)
Output: D4 with rank columns and a quota-based shortlist flag

WHY WITHIN CLASS. The rank metric is gnina's Vina-style affinity, comparable
only among molecules docked the same way. Each class docks against a different
reactive-atom constraint, and classes differ in heavy-atom count, which a
Vina-style score tracks. Sorting all 1,683 survivors together would mostly
rediscover which warhead is biggest and greasiest and hand the shortlist to one
chemotype. Each class is ranked against itself and contributes a fixed quota, so
the shortlist is a designed comparison ACROSS chemotypes — the question T_4
exists to answer, since the core is fixed and the warhead is the variable.

WHICH CLASS: THE ADDUCT'S, NOT THE WARHEAD'S (D0029). chloroacetamide,
sulfamate_acetamide and sulfonate_acetamide differ only in what leaves, so all
187 R-groups give one IDENTICAL adduct and their affinities are equal to the
last decimal. Quota'ing them as three classes spent 9 shortlist slots on 3
molecules and sent 6 redundant systems through MM-GBSA. Grouping on
`adduct_class` collapses them to one and gives T_4 seven post-reaction
chemotypes, not nine. Warhead class is retained on every row as the synthetic
ROUTE — three routes to one bound molecule is a result worth showing, and the
reactivity triage legitimately distinguishes them on kinetics.

LIGAND EFFICIENCY IS ADVISORY. LE = -affinity / heavy_atoms is reported because
it is the standard size correction and reviewers will ask. It is not the rank
metric.

THE RANKING IS NOT VALIDATED (D0031). On class-matched decoys the covalent gate
reads ROC-AUC 0.537, indistinguishable from chance, so this shortlist is an
ordering the pipeline produced rather than evidence of binding. The verdict
travels on every ranked row.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                     # noqa: E402
from shared import rank_shortlist as rs          # noqa: E402
from shared import warhead_library as wl         # noqa: E402

log = logging.getLogger("t4-rank")

EXPERIMENT = "04_t4_combinatorial"
CONFIG = REPO / "config" / "approaches" / "t4_combinatorial.yaml"
RANK_METRIC = "affinity_kcal"      # LOWER is better
STRATUM = "covalent"


def load_ranking_config() -> dict:
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    r = cfg.get("ranking") or {}
    return {"per_class_quota": int(r.get("per_class_quota", 3)),
            "min_docked": int(r.get("min_docked_for_meaningful_rank", 20))}


def main() -> None:
    ap = argparse.ArgumentParser(description="T_4 step 8: rank within adduct class.")
    ap.add_argument("--quota", type=int, default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    cfg = load_ranking_config()
    quota = args.quota or cfg["per_class_quota"]

    df, frame_path = dio.latest_frame(EXPERIMENT, "t4")
    log.info("loaded %s (%d rows)", frame_path.name, len(df))
    if RANK_METRIC not in df.columns:
        raise SystemExit(f"frame has no {RANK_METRIC!r} — run 03_covalent_dock.py first")
    if "dock_id" not in df.columns:
        raise SystemExit("frame has no `dock_id`; molecule identity is required "
                         "to rank distinct molecules rather than routes (D0029)")

    # Post-reaction identity, read from the library rather than hard-coded here.
    lib = wl.load()
    if "adduct_class" not in lib.columns:
        raise SystemExit("warhead library has no `adduct_class` column (D0029)")
    df = df.copy()
    df["adduct_class"] = df["warhead_class"].map(
        dict(zip(lib["class_id"], lib["adduct_class"])))
    unmapped = df["adduct_class"].isna() & df["warhead_class"].notna()
    if unmapped.any():
        raise SystemExit(
            f"{int(unmapped.sum())} row(s) carry a warhead class the library "
            "does not map to an adduct class")
    n_w = df["warhead_class"].nunique()
    n_a = df["adduct_class"].nunique()
    log.info("%d warhead classes -> %d distinct adduct classes (D0029)", n_w, n_a)

    gated = rs.attach_gate(df, STRATUM, RANK_METRIC)
    # Rank MOLECULES, not rows: three routes reach the same acetamide adduct, so
    # ranking rows gave a top-3 of one molecule listed three times (D0029).
    ranked = rs.rank(gated, metric=RANK_METRIC, group_col="adduct_class",
                     min_docked=cfg["min_docked"], identity_col="dock_id")
    final = rs.shortlist(ranked, quota=quota)

    n_short = int(final["shortlist"].sum())
    n_mols = final.loc[final["shortlist"], "dock_id"].nunique() \
        if "dock_id" in final.columns else n_short

    out = dio.write_full_frame(
        final, approach="t4", experiment=EXPERIMENT, stage="t4_rank_within_class",
        params={"rank_metric": f"{RANK_METRIC} (lower better)",
                "ranking_scope": "within ADDUCT class (D0029)",
                "n_warhead_classes": int(n_w),
                "n_adduct_classes": int(n_a),
                "per_class_quota": quota,
                "min_docked_for_meaningful_rank": cfg["min_docked"],
                "gate_verdict": str(final["gate_verdict"].iloc[0]),
                "rank_validated": bool(final["rank_validated"].iloc[0]),
                "n_shortlisted": n_short,
                "n_distinct_molecules_shortlisted": int(n_mols)},
        inputs={"d4_frame": frame_path})

    print(f"\nT_4 within-adduct-class ranking -> {out}")
    print(f"  shortlisted {n_short} rows = {n_mols} DISTINCT molecules "
          f"across {n_a} adduct classes")
    if n_short != n_mols:
        print(f"  ({n_short - n_mols} extra rows are ALTERNATIVE SYNTHETIC ROUTES "
              "to a shortlisted molecule, not extra candidates — D0029)")
    print()
    print(rs.summarise(final, RANK_METRIC))

    print("\n  shortlist (ranked within adduct class; NOT comparable across classes):")
    s = final[final["shortlist"]].sort_values(["rank_group", "rank"])
    for _, r in s.iterrows():
        route = r["warhead_class"]
        print(f"    {str(r['rank_group']):22s} #{int(r['rank'])}  "
              f"{r[RANK_METRIC]:6.2f} kcal/mol  via {route:22s} {r['candidate_id']}")


if __name__ == "__main__":
    main()
