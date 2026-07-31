"""
Purpose: T_2 (ATRA analogues, CReM) — T5 physics rescoring (MM-GBSA) on the shortlist.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the latest frame (post ranking, with `shortlist`)
Output: the frame with dG columns; per-candidate Amber working directories

THE STEP THE PLAN CALLED T5. `docs/approaches/t2.md` step 8: "Physics rescoring
on survivors... computes a higher-fidelity binding estimate on the few
candidates that warrant it... applies no filter of its own." That last clause is
honoured — this stage scores and records, it never stamps `rejected_at`.

The run lives in `shared.mmgbsa_run` so T_1 and T_2 execute one function.

THIS APPROACH: ATRA is a retinoic ACID, so most of this neighbourhood is anionic at pH 7.4.

WHAT IS AND IS NOT DELIVERED. The plan's step 8 lists four things: MM-GBSA,
short explicit-solvent MD, an AI cofold pose with physics relaxation, and
anti-target/selectivity docking. Only the FIRST exists. So "ΔG_bind ±
uncertainty" is reported WITHOUT the uncertainty rather than with an invented
one: a single minimisation has no error bar, and the MD that would supply one is
not in this repo.

COMPARABLE ACROSS THE WHOLE APPROACH. Unlike T_4's covalent dG, which carries a
constant bond term cancelling only within a warhead class (D0020), this dG has
no such term.

READ IT WITH D0031. These candidates were chosen by docking, which is at chance
on this receptor against class-matched decoys. This is therefore an INDEPENDENT
estimate rather than confirmation of the docking ranking; where the two
disagree, the disagreement is the result.
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

log = logging.getLogger("t2-mmgbsa")

EXPERIMENT = "02_t2_atra_crem"


def main() -> None:
    ap = argparse.ArgumentParser(description="T_2 (ATRA analogues, CReM): MM-GBSA rescoring.")
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

    merged, out, results, failures, changed = runner.run(
        selection_col=args.selection_col, experiment=EXPERIMENT, approach="t2",
        workers=args.workers, limit=args.limit)

    print(f"\nT_2 (ATRA analogues, CReM) MM-GBSA (T5 physics rescoring) -> "
          f"{out if out else '(no frame written — partial run)'}")
    print(f"  scored {len(results)}, failed {len(failures)}")
    print(f"  {changed} candidate(s) changed charge at pH 7.4\n")
    if results:
        r = pd.DataFrame(results).sort_values("dG_kcal")
        rank_by_dock = dict(zip(merged["candidate_id"], merged.get("rank", [])))
        print(f"  {'candidate':22s} {'dG (kcal/mol)':>13s} {'dock rank':>10s}")
        print("  " + "-" * 48)
        for _, x in r.iterrows():
            print(f"  {x['candidate_id']:22s} {x['dG_kcal']:13.2f} "
                  f"{str(rank_by_dock.get(x['candidate_id'], '-')):>10s}")
        print("\n  No uncertainty is reported: a single minimisation has none,")
        print("  and the MD the plan pairs with this does not exist here.")
    for f in failures:
        print(f"  FAILED {f['candidate_id']}: {f['mmgbsa_error'][:120]}")


if __name__ == "__main__":
    main()
