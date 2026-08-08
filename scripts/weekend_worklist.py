"""
Purpose: the priority order for the weekend run — 5 best overall, then acrylamide and the bromines.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: the 2.2.0 mode-based ranking (nac_v3)
Output: 00_outputs/blacksmith/attack_sweep/worklist_<N>.csv + a plain id list

@tt8804: *"wed want to sweep through the 5 best and then as many of the the
acrylamide and the two bromine warheds and then send all survivors to MD."*

THE TWO BROMINE WARHEADS ARE bdhi_c4 AND bdhi_c5 -- bromo-dihydroisoxazoles,
SN2 ring-opening, the only two classes whose reactive SMARTS contains Br. Named
by that property rather than by memory, so the list cannot quietly pick the
wrong classes.

THE ORDER IS THE POINT, NOT THE LENGTH. Those three classes hold ~5,990
molecules and the weekend fits a couple of hundred, so this emits a PRIORITY
ORDER and the workers take as much of it as the time allows. What gets done is
then a function of throughput, not of a number guessed in advance.

  rank 1-5    the five best overall, across ALL classes -- the deliverable's
              own top-5, so it is never crowded out by the class quotas below
  then        acrylamide, bdhi_c4, bdhi_c5, INTERLEAVED by class rank

INTERLEAVED, NOT CONCATENATED. Listing all acrylamides first would spend the
whole weekend on one class if throughput fell short, and the two bromine classes
would get nothing. Round-robin means a shortfall costs each class proportionally
instead of erasing two of them.

DISTINCT MOLECULES ONLY. The ranking is over binding MODES, so a molecule can
appear more than once; it enters this list once, on its best mode, for the same
reason the elevation queue collapses (you synthesise a molecule, not a mode).
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402

log = logging.getLogger("worklist")
RANK = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/rank_v2")
OUT = sout.Topic("blacksmith", "attack_sweep")

N_BEST = 5


def bromine_classes() -> list[str]:
    """The classes whose reactive SMARTS actually contains bromine."""
    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    return sorted(r.class_id for r in wh.itertuples()
                  if "Br" in str(r.reactive_atom_smarts))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--score", default="enrichment_conditional")
    ap.add_argument("--out-list", default=None)
    # ASYMMETRIC BY DESIGN (@tt8804): "we are most interested in the 3 warhead
    # classes I mentioned and then the rest maybe just the top few."
    #
    # A flat per-class number is misleading here because the classes are wildly
    # different sizes: 50 is 1% of acrylamide's 4,216 molecules and 35% of
    # bdhi_c4's 138. A flat 50 would also cost 450 sweeps -- 88 GPU-hours on
    # three cards, against a ~44 hour weekend, before any MD at all.
    ap.add_argument("--per-class-priority", type=int, default=50,
                    help="acrylamide, bdhi_c4, bdhi_c5")
    ap.add_argument("--per-class-other", type=int, default=5,
                    help="the remaining classes — the top few only")
    ap.add_argument("--limit", type=int, default=0,
                    help="hard cap on the whole worklist, applied last")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    frames = []
    for tier in ("T3", "T4"):
        fs = sorted(glob.glob(str(RANK / f"rank_v2_{tier}_{args.score}_*.csv")),
                    key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
        if fs:
            d = pd.read_csv(fs[-1])
            d["tier"] = tier
            frames.append(d)
            log.info("%s: %s (%d rows)", tier, Path(fs[-1]).name, len(d))
    if not frames:
        raise SystemExit(f"no {args.score} ranking under {RANK}")
    d = pd.concat(frames, ignore_index=True)

    # one row per MOLECULE, on its best mode
    if "parent_ident" not in d.columns:
        d["parent_ident"] = d["ident"]
    d["parent_ident"] = d["parent_ident"].fillna(d["ident"])
    d = (d.sort_values("class_rank")
           .groupby("parent_ident", as_index=False).first())

    br = bromine_classes()
    log.info("bromine warhead classes, by SMARTS: %s", br)
    priority = ["acrylamide"] + br

    # "TOP 5 OVERALL" MUST BE BY SCORE, NOT BY class_rank.
    #
    # class_rank is rank WITHIN a warhead class, so nsmallest(5, "class_rank")
    # returns five rows all holding rank 1 -- the winners of five arbitrary
    # classes, chosen by whatever order pandas happened to produce. It is a
    # sample of class winners wearing the name of a top-5.
    score_col = args.score if args.score in d.columns else None
    if score_col is None:
        raise SystemExit(f"{args.score} not in the ranking; cannot take a top-5")
    best = d.nlargest(N_BEST, score_col).copy()
    best["why"] = f"top-{N_BEST} overall by {score_col}"
    taken = set(best.parent_ident)

    others = [c for c in sorted(d.warhead_class.dropna().unique())
              if c not in priority]
    quota = ({c: args.per_class_priority for c in priority} |
             {c: args.per_class_other for c in others})
    log.info("allocation: %s", ", ".join(f"{c}={quota[c]}" for c in priority + others))

    pools = {c: d[(d.warhead_class == c) & (~d.parent_ident.isin(taken))]
                .sort_values("class_rank").head(quota[c]).to_dict("records")
             for c in quota}
    # Round-robin over the PRIORITY classes first so a shortfall costs each of
    # them proportionally rather than erasing two, then the rest.
    rows = []
    for group in (priority, others):
        while any(pools[c] for c in group):
            for c in group:
                if pools[c]:
                    r = pools[c].pop(0)
                    r["why"] = f"{c} rank {int(r['class_rank'])}"
                    rows.append(r)
    rest = pd.DataFrame(rows)

    out = pd.concat([best, rest], ignore_index=True)
    if args.limit:
        out = out.head(args.limit)
    out["priority"] = range(1, len(out) + 1)
    cols = [c for c in ("priority", "parent_ident", "ident", "warhead_class",
                        "tier", "class_rank", args.score, "consensus", "QED", "why")
            if c in out.columns]
    dest = OUT.write("worklist", ".csv")
    out[cols].to_csv(dest, index=False)

    # EXTRA MODES WORTH SWEEPING (@tt8804).
    #
    # The sweep would otherwise take mode 0 only. Measured over 4,482 molecules:
    # 15% are multi-mode, and for 9% of THOSE the dominant mode is NOT the
    # best-anchored one -- with a median anchoring gain of 0.503 when it differs,
    # on a mode holding a median 13% of poses. That is sulfopin's failure in a
    # new place: the reactive mode is the minority one, and consensus will never
    # surface it.
    #
    # So: mode 0 always, PLUS any mode beating it by more than MODE_GAIN on
    # anchoring. The population floor already guarantees such a mode is real,
    # which keeps realism ahead of attack geometry rather than behind it.
    MODE_GAIN = 0.2
    v3 = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(
        "/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/nac_v3/agg_s*_*.csv"))],
        ignore_index=True)
    v3 = v3[v3.status == "ok"]
    pairs = []
    for pid in out.parent_ident:
        g = v3[v3.parent_ident == pid]
        if g.empty:
            pairs.append((pid, 1)); continue
        g = g.sort_values("consensus", ascending=False)
        dom = g.iloc[0]
        pairs.append((pid, 1))                     # mode 0 -> pose_rank 1
        for r in g.iloc[1:].itertuples():
            if (r.anchor_quality_max - dom.anchor_quality_max) > MODE_GAIN:
                pairs.append((pid, int(r.mode) + 1))
    lst = Path(args.out_list) if args.out_list else dest.with_suffix(".txt")
    lst.write_text("\n".join(f"{a} {b}" for a, b in pairs) + "\n")
    extra = len(pairs) - len(out)
    log.info("worklist: %d molecules -> %d sweeps (%d extra modes)",
             len(out), len(pairs), extra)

    print(f"\nWeekend worklist — {len(out):,} molecules in priority order\n")
    print(out[cols].head(12).to_string(index=False))
    print(f"\n  first {N_BEST}: the best overall, whatever class they are in")
    print(f"  then round-robin over {', '.join(priority)}")
    print("\n  class composition:")
    for c, n in out.warhead_class.value_counts().items():
        tag = "  <- priority" if c in priority else ""
        print(f"    {c:<22} {n}{tag}")
    sweeps = len(pairs) if "pairs" in dir() else len(out)
    print(f"\n  {len(out)} molecules; sweep cost ~{len(out)*35/60/3:.0f} h on 3 GPUs")
    print(f"\n  -> {dest}\n  -> {lst}")


if __name__ == "__main__":
    main()
