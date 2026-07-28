"""
Purpose: T_4 step 6 — covalent docking of survivors through the shared protocol.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27 (migrated to the shared runner 2026-07-28)
Input: the latest D4 frame (post reactivity triage)
Output: D4 with dock columns; poses under append_only/

PARITY IS THE POINT (control S3). The run lives in `shared.covalent_dock_run`
and is byte-identical to what T_3 executes — same pinned gnina binary, same
parameters, same adduct transform, same protocol fingerprint. This file supplies
only T_4's identity. Until 2026-07-28 this logic was duplicated here and in
T_3's stage; D0030 collapsed them, because two scripts that begin identical do
not stay that way.

THE ADDUCT IS WHAT GETS DOCKED (D0022, D0030). Chlorine gone, sulfamate's twelve
atoms gone, acrylamide's alkene saturated — the quinones alone pass through
untransformed, because their adduct re-aromatizes onto an sp2 carbon.

ONE DOCK PER DISTINCT ADDUCT (D0029). chloroacetamide, sulfamate_acetamide and
sulfonate_acetamide are SN2 at the same CH2 and differ only in what leaves, so
all three give an IDENTICAL adduct — verified for all 198 R-groups. They are
docked once and mapped back to every class that shares the product, which makes
their scores identical BY CONSTRUCTION. That is the honest representation: what
distinguishes those warheads is kinetics, which is the reactivity window's
business and not docking's. Post-reaction, T_4 covers seven chemotypes, not nine.

RANK METRIC (D0015, re-measured in D0028). gnina's Vina-style `affinity`
(kcal/mol, LOWER better), NOT `CNNaffinity`. On adduct forms affinity gives
ROC-AUC 0.718 against CNNaffinity's 0.392 and EF1% 19.0 against 0.0. Note that
D0028 withdrew D0015's decisive claim — affinity's CI now INCLUDES 0.5 — so the
metric choice stands on the comparison, not on a significant interval.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import covalent_dock_run as runner    # noqa: E402

log = logging.getLogger("t4-dock")

EXPERIMENT = "04_t4_combinatorial"


def main() -> None:
    ap = argparse.ArgumentParser(description="T_4 step 6: covalent docking.")
    ap.add_argument("--limit", type=int, default=None,
                    help="dock only the first N survivors (smoke testing)")
    ap.add_argument("--gpus", default=None,
                    help="comma-separated device ids; pass this whenever another "
                         "docking job is running, since the auto-search cannot "
                         "see gnina's ~500 MiB footprint")
    ap.add_argument("--results-name", default="results_adduct.jsonl",
                    help="resumable JSONL; use a NEW name after a protocol "
                         "change so poses are not resumed across fingerprints")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    gpus = [int(g) for g in args.gpus.split(",")] if args.gpus else None
    merged, out, proto, survivors, n_docked = runner.run(
        experiment=EXPERIMENT, approach="t4", frame_prefix="D4",
        limit=args.limit, results_name=args.results_name, gpus=gpus)

    print(f"\nT_4 covalent docking -> "
          f"{out if out else '(no frame written — partial run)'}")
    print(f"  docked successfully {n_docked} / {len(survivors)}")
    print(f"  protocol            {proto.version.strip()}")
    print(f"  fingerprint         {proto.fingerprint()[:16]}")
    if n_docked:
        print("\n  best affinity per warhead class (kcal/mol, lower better):")
        d = merged.dropna(subset=["affinity_kcal"])
        for cls, g in d.groupby("warhead_class"):
            print(f"    {cls:22s} best {g['affinity_kcal'].min():6.2f}   "
                  f"median {g['affinity_kcal'].median():6.2f}   n={len(g)}")
        print("\n  NOTE (D0029): the three acetamide classes above share one")
        print("  adduct, so their numbers are identical by construction.")


if __name__ == "__main__":
    main()
