#!/usr/bin/env python3
"""
Purpose: the modes the screen ranked and never swept, in global rank order.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-11
Input: rank_v2 + attack_sweep, joined on (parent_ident, mode)
Output: 00_outputs/blacksmith/sweep_gaps/sweep_gaps_<N>.csv

FILLING IN THE GAPS (@tt8804, #53). The 2.2.0 sweep took mode 0, once per
molecule -- 233 of the 239 modes it ran. The ranking is per mode. So the modes
that rank highest were, in several cases, never simulated at all: five rank FIRST
in their warhead class. This emits what is missing, best-ranked first, so a run
can work down it.

ALREADY-SWEPT MODES ARE EXCLUDED BY (parent_ident, mode), NEVER BY `ident`. Mode 0
is the bare ident in the sweep table and `_m0` in the rank table, so an `ident`
match would think every simulated mode was missing and re-run all 239 of them
(`shared/mode_key.py`).

THE ORDER IS GLOBAL, AND GLOBAL IS NOT NEUTRAL. Global rank is computed on
`conditional_eb`, and `conditional_eb` exists for **T4 only** -- 0 of 4,607 T3
rows carry it. A global list is therefore a T4 list with T3 sorted to the bottom,
which is a fact about the score's coverage and not a judgement about T3. Ranking
across warhead classes is separately biased: the SN2 angular criterion is
stricter than the perpendicular one (#47). Both are stamped into the output.

EVERY ROW IS CHECKED AGAINST ITS POSE FILE BEFORE IT IS EMITTED. The runner asks
for a pose by `--pose-rank`, so a row whose SDF has no pose at that rank, or whose
pose there carries a different `mode`, is dropped here rather than failing on a
GPU an hour later.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import mode_ranking as mr                     # noqa: E402
from shared import outputs as sout                        # noqa: E402

log = logging.getLogger("sweep-gaps")
POSES = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/nac_v3_poses")


def resolvable(parent: str, mode: int) -> int | None:
    """The `pose_rank` whose pose carries this `mode`, or None.

    Read by identity, both ways: the pose must exist AND its own `mode` property
    must match. `pose_rank - 1 == mode` holds for every pose on disk today and is
    not guaranteed by construction, so it is verified rather than assumed.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    f = POSES / f"{parent}.sdf"
    if not f.is_file():
        return None
    for m in Chem.SDMolSupplier(str(f), removeHs=False, sanitize=False):
        if m is None or not (m.HasProp("pose_rank") and m.HasProp("mode")):
            continue
        if int(m.GetProp("mode")) == mode:
            return int(m.GetProp("pose_rank"))
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--limit", type=int, default=400,
                    help="how many rows to emit; the list is ordered so a "
                         "truncated run is still the best-ranked ones")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    d = mr.gather()
    if d.empty:
        raise SystemExit("no rank tables")

    gap = d[(~d.sent) & d.global_rank.notna()].sort_values("global_rank")
    log.info("%d modes ranked, %d already sent, %d in the gap",
             len(d), int(d.sent.sum()), len(gap))

    rows, checked = [], 0
    for _, x in gap.iterrows():
        if len(rows) >= args.limit:
            break
        checked += 1
        pr = resolvable(str(x.parent_ident), int(x["mode"]))
        if pr is None:
            log.debug("%s: no pose for mode %s", x.ident, x["mode"])
            continue
        rows.append({
            "ident": x.ident, "parent_ident": x.parent_ident,
            "mode": int(x["mode"]), "pose_rank": pr,
            "global_rank": int(x.global_rank),
            "class_rank": int(x.class_rank) if pd.notna(x.class_rank) else None,
            "warhead_class": x.warhead_class,
            "conditional_eb": x.conditional_eb,
            "tier": x.get("tier"),
            "order_basis": "global rank on conditional_eb (T4 only; "
                           "cross-class comparison is biased, #47)",
        })

    out = pd.DataFrame(rows)
    dest = sout.Topic("blacksmith", "sweep_gaps").write("sweep_gaps", ".csv")
    out.to_csv(dest, index=False)
    print(f"\n  {len(out)} runnable of {checked} inspected -> {dest}")
    if len(out):
        print(out.head(10)[["ident", "warhead_class", "global_rank",
                            "class_rank", "pose_rank"]].to_string(index=False))


if __name__ == "__main__":
    main()
