"""
Purpose: re-score named candidates through the near-attack gate at a high run
         count, to test whether their 200-run enrichment converges.
Author: @tt8804 (with Claude Code)
Date: 2026-08-06
Input: candidate idents (as they appear in the T_3/T_4 frames)
Output: <outdir>/nac_converge_one_<N>.csv, one row per ident per --nrun

WHY THIS EXISTS RATHER THAN `nac_rank.py --refine-top`. Two reasons, and the
second is a defect worth writing down.

1. `--refine-top N` re-scores the top N of the whole ranking. Asking about ONE
   molecule should not require re-scoring its neighbours.

2. `nac_rank.refine()` builds its resume set with `_ids_in(REFINE.dir, ...)`,
   which collects EVERY ident already written to a chunk file **regardless of
   status**. A row written as `failed: <reason>` therefore counts as done, and
   the molecule is never retried. Measured here: a first attempt from an env
   without `gemmi` wrote five `failed:` rows, and the immediate re-run in the
   correct env reported "5 assigned, 0 to do" and scored nothing. Append-only
   means those rows cannot be removed, so the molecule is permanently
   un-refinable through that path. This script does not inherit that: it takes
   the idents it is given and always measures them.

THE 200-RUN AND 2,000-RUN NUMBERS ARE DIFFERENT MEASUREMENTS (D0068). They are
written to their own file and must never be concatenated into one column with
the screen's.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import nac_rank as nr          # noqa: E402
import nac_screen as ns        # noqa: E402

log = logging.getLogger("medchem-converge")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--idents", nargs="+", required=True)
    ap.add_argument("--nrun", type=int, nargs="+", default=[2000])
    ap.add_argument("--gpu", default="7")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--version", type=int, default=1)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    out = args.outdir / f"nac_converge_one_{args.version}.csv"
    if out.exists():
        raise SystemExit(f"{out} exists -- bump --version rather than overwriting")

    by_id = {c.ident: c for c in nr.load_candidates()}
    missing = [i for i in args.idents if i not in by_id]
    if missing:
        raise SystemExit(f"not in the candidate frames: {missing}")

    rec = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    rows = []
    for nrun in args.nrun:
        for ident in args.idents:
            r = nr.score_one(by_id[ident], rec, nrun, args.gpu)
            r["nrun_requested"] = nrun
            log.info("%s @ %d runs -> %s", ident, nrun, r.get("status"))
            rows.append(r)

    df = pd.DataFrame(rows)
    # A failed row is reported, not dropped -- and, unlike the refine path, it
    # is not a reason the molecule can never be measured again.
    args.outdir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(df.to_string(index=False))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
