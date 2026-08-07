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
    ap.add_argument("--limit", type=int, default=0,
                    help="cap the worklist (@tt8804: sweep the top 50)")
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

    best = d.nsmallest(N_BEST, "class_rank").copy()
    best["why"] = "top-5 overall"
    taken = set(best.parent_ident)

    # round-robin the three priority classes so a shortfall costs each of them
    # proportionally rather than erasing two
    pools = {c: d[(d.warhead_class == c) & (~d.parent_ident.isin(taken))]
                .sort_values("class_rank").to_dict("records")
             for c in priority}
    rows, i = [], 0
    while any(pools.values()):
        for c in priority:
            if pools[c]:
                r = pools[c].pop(0)
                r["why"] = f"{c} rank {int(r['class_rank'])}"
                rows.append(r)
        i += 1
    rest = pd.DataFrame(rows)

    out = pd.concat([best, rest], ignore_index=True)
    if args.limit:
        out = out.head(args.limit)
    out["priority"] = range(1, len(out) + 1)
    cols = [c for c in ("priority", "parent_ident", "ident", "warhead_class",
                        "tier", "class_rank", args.score, "QED", "why")
            if c in out.columns]
    dest = OUT.write("worklist", ".csv")
    out[cols].to_csv(dest, index=False)

    lst = Path(args.out_list) if args.out_list else dest.with_suffix(".txt")
    lst.write_text("\n".join(out.parent_ident.astype(str)) + "\n")

    print(f"\nWeekend worklist — {len(out):,} molecules in priority order\n")
    print(out[cols].head(12).to_string(index=False))
    print(f"\n  first {N_BEST}: the best overall, whatever class they are in")
    print(f"  then round-robin over {', '.join(priority)}")
    print("\n  class composition of the first 100:")
    for c, n in out.head(100).warhead_class.value_counts().items():
        print(f"    {c:<22} {n}")
    print(f"\n  -> {dest}\n  -> {lst}")


if __name__ == "__main__":
    main()
