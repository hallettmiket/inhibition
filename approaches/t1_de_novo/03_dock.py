"""
Purpose: T_1 (de novo, DiffSBDD) step 3 — non-covalent docking of the surviving candidates.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the latest D1 frame (post annotate-and-gate)
Output: D1 with `vina_affinity`; poses under append_only/

PARITY WITH T_1/T_2's SIBLING. The run lives in `shared.noncovalent_dock_run`
and is byte-identical to what the other reversible approach executes — same
Vina-GPU build, same expanded box, same pinned search_depth 20 (D0017). This
file supplies only this approach's identity, for the same reason T_3 and T_4
share `covalent_dock_run`: two scripts that begin identical do not stay that way.

READ THE RANKING WITH ITS GATE. The non-covalent enrichment gate on this pocket
is ROC-AUC 0.535 with a CI spanning 0.215-0.855 and EF1% of 0.0 (D0016) — on
this shallow, solvent-exposed site, docking barely distinguishes known binders
from property-matched decoys. These scores order the candidates; they do not
establish that the ordering is meaningful, and the GUI must show the verdict
beside them.

THIS APPROACH: seed-free, so its agreement with a seeded approach is not ancestral.

GPU IS PASSED IN, NOT DISCOVERED. Covalent docking may be running concurrently
on other devices, and its allocator cannot see gnina's ~500 MiB footprint. The
device is therefore chosen by the operator.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import noncovalent_dock_run as runner    # noqa: E402

log = logging.getLogger("t1-dock")

EXPERIMENT = "01_t1_de_novo"


def main() -> None:
    ap = argparse.ArgumentParser(description="T_1 (de novo, DiffSBDD) step 3: non-covalent docking.")
    ap.add_argument("--gpu", type=int, required=True,
                    help="GPU device id; must not be one a covalent dock is using")
    ap.add_argument("--limit", type=int, default=None,
                    help="dock only the first N survivors (smoke testing)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    merged, out, survivors, n_docked, elapsed = runner.run(
        experiment=EXPERIMENT, approach="t1", frame_prefix="D1",
        gpu=args.gpu, limit=args.limit)

    print(f"\nT_1 (de novo, DiffSBDD) non-covalent docking -> {out}")
    print(f"  docked successfully {n_docked} / {len(survivors)}")
    print(f"  engine              Vina-GPU 2.1, search_depth {runner.SEARCH_DEPTH} (D0017)")
    print(f"  elapsed             {elapsed:.1f} s on GPU {args.gpu}")
    if n_docked:
        d = merged.dropna(subset=["vina_affinity"])
        print(f"\n  vina_affinity (kcal/mol, lower better) over {len(d)} candidates:")
        print(f"    best   {d['vina_affinity'].min():6.2f}")
        print(f"    median {d['vina_affinity'].median():6.2f}")
        print(f"    worst  {d['vina_affinity'].max():6.2f}")
    print("\n  CAVEAT (D0016): non-covalent enrichment here is ROC-AUC 0.535,")
    print("  CI 0.215-0.855, EF1% 0.0. This ranking is weakly supported.")


if __name__ == "__main__":
    main()
