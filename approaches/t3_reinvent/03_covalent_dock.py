"""
Purpose: T_3 step 3 — covalent docking of the surviving decorations.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the latest D3 frame (post annotate-and-gate)
Output: D3 with dock columns; poses under append_only/

PARITY IS THE POINT (control S3). The run itself lives in
`shared.covalent_dock_run` and is byte-identical to what T_4 executes — same
pinned gnina binary, same parameters, same adduct transform, same protocol
fingerprint. This file supplies only T_3's identity. "We both ran gnina" is not
the claim the integration phase needs; "we ran the same function" is.

T_3 IS SINGLE-WARHEAD. Every candidate is an acrylamide (fixed up front by the
PI, not selected), so the within-class ranking T_4 needs is vacuous here and the
cross-class dedup that merges T_4's three SN2 acetamides has nothing to merge.
What the shared dedup DOES do for T_3 is collapse decorations the generator
proposed more than once.

THE ADDUCT IS SATURATED (D0030). Acrylamide's Michael adduct is
`Cys-S-CH2-CH2-C(=O)NR2` — no re-aromatization is available, so the C=C is gone.
Docking the alkene instead would hand gnina an sp2 carbon and model a planar,
rigid vinyl thioether where the real linker rotates freely. Since acrylamide is
T_3's ONLY warhead, that would have biased the whole approach against exactly
the decorations that need to bend to reach a subsite — which is T_3's question.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import covalent_dock_run as runner    # noqa: E402

log = logging.getLogger("t3-dock")

EXPERIMENT = "03_t3_reinvent"


def main() -> None:
    ap = argparse.ArgumentParser(description="T_3 step 3: covalent docking.")
    ap.add_argument("--limit", type=int, default=None,
                    help="dock only the first N survivors (smoke testing)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    merged, out, proto, survivors, n_docked = runner.run(
        experiment=EXPERIMENT, approach="t3", frame_prefix="D3",
        limit=args.limit)

    print(f"\nT_3 covalent docking -> "
          f"{out if out else '(no frame written — partial run)'}")
    print(f"  docked successfully {n_docked} / {len(survivors)}")
    print(f"  protocol            {proto.version.strip()}")
    print(f"  fingerprint         {proto.fingerprint()[:16]}")
    if n_docked:
        d = merged.dropna(subset=["affinity_kcal"])
        print(f"\n  affinity (kcal/mol, lower better) over {len(d)} candidates:")
        print(f"    best   {d['affinity_kcal'].min():6.2f}")
        print(f"    median {d['affinity_kcal'].median():6.2f}")
        print(f"    worst  {d['affinity_kcal'].max():6.2f}")


if __name__ == "__main__":
    main()
