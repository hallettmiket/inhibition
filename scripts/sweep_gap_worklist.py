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
    ap.add_argument("--only-class", default=None,
                    help="restrict to one warhead class")
    ap.add_argument("--strata", default=None,
                    help="depth ladder: 'lo-hi:n,lo-hi:n' over class rank, "
                         "interleaved so an early stop still spans the range")
    ap.add_argument("--per-class", type=int, default=None,
                    help="balanced mode: take this many best-ranked unswept "
                         "modes from EACH warhead class, interleaved")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    d = mr.gather()
    if d.empty:
        raise SystemExit("no rank tables")

    gap = d[(~d.sent) & d.global_rank.notna()].sort_values("global_rank")
    log.info("%d modes ranked, %d already sent, %d in the gap",
             len(d), int(d.sent.sum()), len(gap))

    if args.only_class:
        gap = gap[gap.warhead_class == args.only_class]
        log.info("restricted to %s: %d unswept modes", args.only_class, len(gap))

    if args.strata:
        # A DEPTH LADDER, TO FIND WHERE A CLASS STOPS BEING WORTH SWEEPING.
        #
        # The docking-derived `viable_fraction` decays smoothly with class rank
        # and then collapses -- for bdhi_c5, 37% at rank 1-50 down to 2% by 175
        # and a median of ZERO past 176, where over half the modes have no pose
        # meeting the near-attack criterion at all. Ranks 101-175 have never been
        # swept, so the knee has never been measured, only predicted.
        #
        # If a swept mode with viable_fraction ~ 0 really does score 0, that is a
        # FREE stopping rule for every class: stop where the median crosses the
        # threshold, and never spend a GPU below it. That is what this ladder
        # tests, which is why it is weighted onto the transition rather than
        # spread evenly.
        #
        # Interleaved across strata: the deliverable is a CURVE, so a run that
        # stops early must still span the range rather than fill the top of it.
        buckets = []
        for part in args.strata.split(","):
            span, n = part.split(":")
            lo, hi = (int(v) for v in span.split("-"))
            b = gap[(gap.class_rank >= lo) & (gap.class_rank <= hi)] \
                .sort_values("class_rank")
            take = b.iloc[:: max(1, len(b) // int(n))].head(int(n))
            log.info("  stratum %s: %d of %d unswept taken", span, len(take), len(b))
            buckets.append(take)
        order, i = [], 0
        while any(i < len(b) for b in buckets):
            for b in buckets:
                if i < len(b):
                    order.append(b.iloc[i])
            i += 1
        gap = pd.DataFrame(order)
        log.info("ladder: %d rows, interleaved across %d strata",
                 len(gap), len(buckets))

    if args.per_class:
        # BALANCED, AND ROUND-ROBIN SO ANY STOPPING POINT IS BALANCED TOO.
        #
        # Equal n PER CLASS, not equal coverage. The library is wildly uneven --
        # acrylamide is 4,835 modes and chloroacetamide 200 -- so equalising
        # coverage would spend almost everything on acrylamide, and equalising
        # count gives each class comparable statistical power, which is what a
        # per-class comparison needs. The trade is explicit: chloroacetamide
        # gains far more coverage per sweep than acrylamide does.
        #
        # Within a class the order is CLASS rank, not global. Cross-class
        # ordering compares scores computed under different bars (#47); within a
        # class it is the criterion's own ordering and is the valid one.
        #
        # Interleaved, so a run that stops early has taken roughly the same
        # number from every class rather than all of the first one.
        by = {c: g.sort_values("class_rank") for c, g in gap.groupby("warhead_class")}
        order, i = [], 0
        while i < args.per_class:
            for c in sorted(by):
                if i < len(by[c]):
                    order.append(by[c].iloc[i])
            i += 1
        gap = pd.DataFrame(order)
        log.info("balanced: %d classes x up to %d = %d rows, interleaved",
                 len(by), args.per_class, len(gap))

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
            "order_basis": (f"depth ladder over {args.only_class} class rank, "
                            f"strata {args.strata}, interleaved"
                            if args.strata else
                            "balanced: equal n per warhead class, ordered by "
                            "CLASS rank within each, interleaved"
                            if args.per_class else
                            "global rank on conditional_eb (T4 only; "
                            "cross-class comparison is biased, #47)"),
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
