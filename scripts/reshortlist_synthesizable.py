"""
Purpose: Rebuild each shortlist so its 25 slots hold synthesizable molecules.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: each approach's latest frame
Output: a new frame version with `shortlist_synth` and `rank_synth`

Run:  python scripts/reshortlist_synthesizable.py [--dry-run] [--n 25]

Issue #1, PI decision: "molecules that fail the synthesizability tests should be
removed from the top 25, or at least placed below."

WHY SORTING THE EXISTING SHORTLIST IS NOT ENOUGH. The shortlists are EXACTLY 25
rows. Demoting the failures inside that set only reorders it -- the same 25
molecules remain, 8 of T_1's still among them, because there is nothing below
the cut to promote. To actually fill 25 slots with molecules the lab would make,
the selection has to reach back into the full scored set (3233 candidates for
T_1) and take the next-best passing ones.

THE ORIGINAL SHORTLIST IS NOT OVERWRITTEN. `shortlist` and `rank` stay exactly
as they were; this adds `shortlist_synth` and `rank_synth` beside them. The
original selection is what every earlier result and decision record refers to,
and silently redefining it would make those references wrong without changing
a word of them.

THIS IS A FILTER ON THE CANDIDATE POOL, NOT A NEW SCORE. Order within the
passing set is the approach's own metric, unchanged. The gate has measured that
metric as not demonstrably enriching (D0041) and partly a size ranking (D0043);
removing molecules that cannot be made does not fix that and is not claimed to.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                      # noqa: E402
from shared import synthesizability as syn        # noqa: E402

log = logging.getLogger("reshortlist")

# THE SELECTION RULE BELONGS TO shared.rank_shortlist, NOT TO THIS SCRIPT.
# A first version took a global top-25 by metric and churned 24 of T_4's 25 --
# despite T_4 having ZERO rule failures. The shortlist is a per-`rank_group`
# QUOTA, and T_4 deliberately spends slots across adduct classes so the list is
# a designed comparison rather than whatever the biggest, greasiest warhead
# produced. Re-deriving that rule here would have silently replaced a designed
# comparison with a leaderboard.
#
# So: drop the rule-failures from the pool, then re-run the approach's OWN rank
# and shortlist functions with the same arguments it uses.
APPROACHES = {
    #  approach: (experiment, metric, group_col, identity_col, quota)
    "t1": ("01_t1_de_novo", "vina_affinity", None, None, 25),
    "t2": ("02_t2_atra_crem", "vina_affinity", None, None, 25),
    "t3": ("03_t3_reinvent", "affinity_kcal", None, None, 25),
    "t4": ("04_t4_combinatorial", "affinity_kcal", "adduct_class", "dock_id", 3),
}


def _targets() -> list[tuple[str, str, str, str | None, str | None, int]]:
    """(label, experiment, metric, group_col, identity_col, quota) per run.

    T_2 IS ONE ROW ABOVE AND FIVE RUNS IN PRACTICE. The reseeding gave each
    seed its own experiment directory, and this table named only ATRA's — so
    du_xu and guo_pfizer were ranked and shortlisted but never rebuilt, and the
    GUI showed their RAW top-25 under a banner claiming the synthesizability
    filter was in force. The seed -> experiment lookup comes from
    `shared.seeds`, which is the single source of it; restating it here is the
    duplication that module's own docstring warns about.
    """
    from shared import seeds as sd

    out = []
    for a, (experiment, metric, group_col, identity_col, quota) in APPROACHES.items():
        if a != "t2":
            out.append((a, experiment, metric, group_col, identity_col, quota))
            continue
        for key in sd.declared_for("t2"):
            try:
                rec = sd.resolve("t2", key, require_radius=False)
            except Exception as exc:  # noqa: BLE001
                log.warning("t2 seed %s does not resolve (%s); skipping", key, exc)
                continue
            out.append((f"t2:{key}", rec["experiment"], metric, group_col,
                        identity_col, quota))
        # The degree-2 sample is a derived run of the ATRA seed into its own
        # directory, not a seed in seeds.yaml, and it needs the rebuild too.
        out.append(("t2:atra_degree2", "02_t2_atra_crem_degree2", metric,
                    group_col, identity_col, quota))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--min-docked", type=int, default=1)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    from shared import rank_shortlist as rs

    for label, experiment, metric, group_col, identity_col, quota in _targets():
        a = label.split(":")[0]          # the approach a frame belongs to
        try:
            df, path = dio.latest_frame(experiment, a)
        except Exception as exc:  # noqa: BLE001 - a seed may not have run yet
            log.warning("%s: no frame under %s (%s); skipping",
                        label, experiment, exc)
            continue
        if metric not in df.columns:
            log.warning("%s: no %s column; skipping (still docking?)",
                        label, metric)
            continue

        df = df.copy()
        df["synth_fail"] = [bool(syn.violations(str(x)))
                            for x in df["canonical_smiles"]]
        df["synth_violations"] = [
            "; ".join(r.name for r in syn.violations(str(x))) or pd.NA
            for x in df["canonical_smiles"]]

        # Failures are removed from the POOL, then the approach's own ranking
        # runs on what is left. Their original rank/shortlist columns are
        # untouched -- every earlier result refers to those, and redefining them
        # in place would make those references wrong without changing a word.
        pool = df[~df["synth_fail"]].copy()
        ranked = rs.rank(pool, metric=metric, group_col=group_col,
                         min_docked=args.min_docked, identity_col=identity_col)
        final = rs.shortlist(ranked, quota=quota)

        keep = final[final["shortlist"]]["candidate_id"]
        df["shortlist_synth"] = df["candidate_id"].isin(set(keep))
        rank_map = dict(zip(final["candidate_id"], final["rank"]))
        df["rank_synth"] = df["candidate_id"].map(rank_map)

        old = df[df.get("shortlist", pd.Series(False, index=df.index)).fillna(False)]
        old_ids, new_ids = set(old["candidate_id"]), set(keep)
        log.info("%s: %d of %d candidates fail a rule; old shortlist had %d",
                 label, int(df["synth_fail"].sum()), len(df),
                 int(old["synth_fail"].sum()) if len(old) else 0)
        log.info("%s:   new shortlist %d — %d carried over, %d promoted, %d dropped",
                 label, len(new_ids), len(old_ids & new_ids),
                 len(new_ids - old_ids), len(old_ids - new_ids))

        if args.dry_run:
            continue
        out = dio.write_full_frame(
            df, approach=a, experiment=experiment,
            stage=f"{a}_shortlist_synth",
            params={"decision": "issue #1 PI decision: failures must not hold "
                                "a top-25 slot",
                    "quota": quota, "metric": metric,
                    "group_col": group_col, "identity_col": identity_col,
                    "rule_module": "shared/synthesizability.py",
                    "note": "original shortlist/rank columns untouched"},
            inputs={"previous_frame": path})
        log.info("%s: wrote %s", label, out.name)


if __name__ == "__main__":
    main()
