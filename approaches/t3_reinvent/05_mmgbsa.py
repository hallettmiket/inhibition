"""
Purpose: T_3 step 5 — T5 physics rescoring (covalent MM-GBSA) on the shortlist.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the latest D3 frame (post ranking, with `shortlist` and `dock_id`)
Output: D3 with dG columns; per-candidate Amber working directories

The run lives in `shared.mmgbsa_run.run_covalent`, byte-identical to what T_4
executes — same link-atom scheme, same junction parameters, same implicit
solvent, same minimiser. Control S3 applies to the physics tier exactly as it
applies to docking: a within-covalent comparison is only defensible if both
approaches computed dG the same way.

T_3 IS SINGLE-WARHEAD, WHICH HELPS TWICE HERE. Every candidate is an acrylamide,
so the D0020 caveat that dG is comparable only within a warhead class is
satisfied trivially — T_3's whole shortlist is one class and its dG values are
mutually comparable.

AND IT PARAMETERISES ONLY BECAUSE OF D0030. The junction parameters cover an
sp3 attachment carbon. Acrylamide's adduct is sp3 only after its acceptor C=C is
saturated; docked as the alkene it would present a conjugated sp2 carbon and hit
the same missing `2C - S - cc` angle that currently fails T_4's naphthoquinone,
bdhi and sNAr classes. The D0030 fix is what makes this stage runnable.

READ IT WITH D0031. These candidates were selected by docking, which is at
chance on this receptor against class-matched decoys, so this is an independent
estimate rather than confirmation of the ranking.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import compute                      # noqa: E402
from shared import mmgbsa_run as runner          # noqa: E402

log = logging.getLogger("t3-mmgbsa")

EXPERIMENT = "03_t3_reinvent"


def main() -> None:
    ap = argparse.ArgumentParser(description="T_3 step 5: covalent MM-GBSA.")
    ap.add_argument("--selection-col", default="shortlist",
                    help="which boolean column selects candidates; "
                         "use shortlist_synth for the synthesizable list")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=compute.MAX_CPU_WORKERS,
                    help="concurrent candidates; the project budget "
                         "lives in shared/compute.py")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    merged, out, results, failures = runner.run_covalent(
        selection_col=args.selection_col, experiment=EXPERIMENT, approach="t3",
        workers=args.workers, limit=args.limit)

    print(f"\nT_3 covalent MM-GBSA (T5 physics rescoring) -> "
          f"{out if out else '(no frame written — partial run)'}")
    print(f"  scored {len(results)}, failed {len(failures)}\n")
    if results:
        r = pd.DataFrame(results).sort_values("dG_kcal")
        print(f"  {'candidate':22s} {'dG (kcal/mol)':>13s}")
        print("  " + "-" * 38)
        for _, x in r.iterrows():
            print(f"  {x['candidate_id']:22s} {x['dG_kcal']:13.2f}")
        print("\n  All of T_3 is one warhead class, so these dG values are")
        print("  mutually comparable (D0020's caveat is satisfied trivially).")
    for f in failures:
        print(f"  FAILED {f['candidate_id']}: {f['mmgbsa_error'][:120]}")


if __name__ == "__main__":
    main()
