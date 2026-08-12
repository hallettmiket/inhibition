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
#: Set from `--topic` in `main()`. NOT a hardcoded topic: reading poses from one
#: run while ranking the tables of another is the defect that cost 2.2.0 and
#: nearly cost 3.0.0 (see nac_screen_v2.topic_paths).
POSES = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/nac_v4_poses")


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
    ap.add_argument("--v-strata", default=None,
                    help="strata over VIABLE_FRACTION, 'lo-hi:n,...' — the "
                         "estimand is P(productive | viable_fraction), which is "
                         "dimensionless and comparable across classes and targets")
    ap.add_argument("--strata", default=None,
                    help="depth ladder: 'lo-hi:n,lo-hi:n' over class rank, "
                         "interleaved so an early stop still spans the range")
    ap.add_argument("--per-class", type=int, default=None,
                    help="balanced mode: take this many best-ranked unswept "
                         "modes from EACH warhead class, interleaved")
    ap.add_argument("--by-family", action="store_true",
                    help="the production rule: restrict to the warhead families "
                         "in config/target.yaml and take up to `max_depth` "
                         "best-ranked unswept modes from each, interleaved")
    ap.add_argument("--max-depth", type=int, default=None,
                    help="override the per-family cap; defaults to "
                         "sweep_rule.max_depth in config/target.yaml")
    ap.add_argument("--topic", default="nac_v4",
                    help="which run's representative poses to resolve against; "
                         "must be the run whose tables were ranked")
    args = ap.parse_args()
    global POSES
    POSES = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith"
                 ) / f"{args.topic}_poses"
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

    if args.v_strata:
        # THE ESTIMAND IS P(productive | viable_fraction), NOT P(productive | rank).
        #
        # Rank is a within-class ordinal: it is not comparable across classes, it
        # depends on library size, and it is at chance for predicting sweep
        # outcome (AUC 0.518 over 188 swept modes). `viable_fraction` is the
        # fraction of a mode's docked poses meeting the near-attack criterion --
        # dimensionless, defined identically for every class and every target,
        # and free at ranking time. That makes it the parameter a TARGET-AGNOSTIC
        # tool can carry a threshold on; a rank cut cannot travel.
        #
        # THE ZERO STRATUM IS THE PRIZE. 1,869 modes (28% of the library) have
        # viable_fraction == 0 and FIVE of them have ever been swept. If a mode
        # whose docked ensemble contains no near-attack pose cannot sweep
        # productive, that is a free exclusion of 28% of every future screen. If
        # it can, the docking-derived criterion is not a filter at all. Either
        # answer is worth more than any rank curve, and neither is currently
        # measurable: the stratum is unsampled.
        buckets = []
        for part in args.v_strata.split(","):
            span, n = part.split(":")
            lo, hi = (float(v) for v in span.split("-"))
            b = gap[(gap.viable_fraction >= lo) & (gap.viable_fraction <= hi)]
            # Spread the draw across classes so one class cannot carry a stratum,
            # then across rank within a class -- a stratum drawn from one class is
            # a statement about that class, not about viable_fraction.
            parts = []
            for _c, g in b.sort_values(["warhead_class", "class_rank"]) \
                          .groupby("warhead_class"):
                parts.append(g.iloc[:: max(1, len(g) // 6)].head(6))
            take = pd.concat(parts) if parts else b.iloc[0:0]
            take = take.iloc[:: max(1, len(take) // max(1, int(n)))].head(int(n))
            log.info("  v in [%g,%g]: %d taken of %d unswept, %d classes",
                     lo, hi, len(take), len(b), take.warhead_class.nunique())
            buckets.append(take)
        order, i = [], 0
        while any(i < len(b) for b in buckets):
            for b in buckets:
                if i < len(b):
                    order.append(b.iloc[i])
            i += 1
        gap = pd.DataFrame(order)
        log.info("v-ladder: %d rows across %d strata, interleaved",
                 len(gap), len(buckets))

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

    if args.by_family:
        # THE PRODUCTION SELECTION RULE (@tt8804, 2026-08-12).
        #
        # Two decisions, and they do different jobs. SCOPE says which warhead
        # chemistry the lab will synthesise -- three families of the nine classes
        # screened -- and it is cheap: it removes only 14% of the modes, because
        # acrylamide alone is 76% of what is in scope. DEPTH is what makes the
        # campaign runnable: 29,255 in-scope modes at a measured 20.9 min median
        # per 10 ns sweep is 212 days on two GPUs, and 250 per family is 5.4.
        #
        # ORDERED ON `enrichment`, WITHIN A FAMILY. Ranking across families would
        # be defensible here -- unusually, all five in-scope classes carry
        # isotropic_null = 0.0816, so enrichment is on one scale -- but it would
        # hand every slot to acrylamide, which is 4,223 of 4,779 in-scope
        # molecules. The lopsidedness is an artefact of how the library was
        # generated, not evidence about the chemistry, and the deliverable is a
        # top-5 a chemist can pick ACROSS warhead types. So the quota is per
        # family and the diversity is bought explicitly rather than hoped for.
        #
        # Interleaved, so a run stopped early has gone equally deep in each
        # family instead of exhausting the first one.
        from shared import target_config as tc
        fam_of = tc.family_of()
        depth = args.max_depth if args.max_depth is not None else tc.sweep_max_depth()
        before = len(gap)
        gap = gap[gap.warhead_class.isin(fam_of)].copy()
        gap["family"] = gap.warhead_class.map(fam_of)
        log.info("scope: %d of %d unswept modes in %d families (%s)",
                 len(gap), before, len(set(fam_of.values())),
                 ", ".join(sorted(set(fam_of.values()))))
        # THE FLOOR FIRST, THEN THE CAP, AND THE LOG SAYS WHICH ONE BOUND. A
        # family stopped by the floor ran to its chemistry limit; one stopped by
        # the cap ran out of GPU time with candidates still queued. Those are
        # different claims about the shortlist and they must not look alike.
        bfloor = tc.sweep_budget_floor()
        pre = len(gap)
        gap = gap[gap.enrichment >= bfloor]
        log.info("budget floor: enrichment >= %.1f keeps %d of %d unswept modes",
                 bfloor, len(gap), pre)
        by = {f: g.sort_values("enrichment", ascending=False)
              for f, g in gap.groupby("family")}
        order, i = [], 0
        while i < depth:
            for f in sorted(by):
                if i < len(by[f]):
                    order.append(by[f].iloc[i])
            i += 1
        gap = pd.DataFrame(order)
        for f in sorted(by):
            n = len(by[f])
            log.info("  %-16s %5d above floor, %4d taken%s", f, n, min(n, depth),
                     "  <- CAPPED by max_depth (budget)" if n > depth
                     else "  (exhausted at the floor)")
        log.info("by-family: %d rows, floor %.1f, cap %d/family, interleaved",
                 len(gap), bfloor, depth)

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
            "order_basis": ("per family (config scope), ordered by enrichment "
                            "within family, interleaved, capped at max_depth"
                            if args.by_family else
                            f"depth ladder over {args.only_class} class rank, "
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
