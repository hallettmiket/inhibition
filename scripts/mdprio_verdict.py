"""
Purpose: score the MD-priority experiment against the readings fixed before it ran.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: the six completed 100 ns trajectories
Output: 00_outputs/blacksmith/mdprio_reports/mdprio_verdict_<N>.csv + the verdict

`docs/prereg_md_priority.md` fixed ONE readout before any of these ran:

    BPMD occupancy is RANK-CORRELATED with 100 ns RESIDENCE.
    Specifically t4_da2e98512d02 holds and t4_9265b4bff789 does not.

This scores exactly that, on all six, and reports whichever answer comes out.

WHY THIS IS A SEPARATE SCRIPT AND NOT A PARAGRAPH. The readings were fixed in
advance precisely so the answer could not be chosen after seeing the data, and
the most reliable way to break that is to compute the number by hand alongside a
narrative. The prereg's table is transcribed here as a literal, and the verdict
is looked up in it.

AND THE THING NOT TO DO. While checking something else, BPMD occupancy turned out
to predict ATTACK GEOMETRY at rho = +0.900 (p = 0.037) where it predicts
residence at +0.410. That is a better story and it is NOT this experiment. It was
not predicted, it is n = 5, and it gets its own pre-registration
(`docs/prereg_attack_sweep.md`) on its own molecules. This script reports the
residence readout as written -- including as a null.
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                  # noqa: E402

log = logging.getLogger("verdict")
MD = Path("/data/lab_vm/modifiable/inhibition/md_residence_3ikd")
OUT = sout.Topic("blacksmith", "mdprio_reports")

#: The six molecules and the BPMD occupancy each was selected on, transcribed
#: from docs/prereg_md_priority.md. Not re-derived: the prediction was made on
#: THESE numbers, so scoring it against freshly computed ones would be scoring a
#: different prediction.
PREREG = {
    "t4_da2e98512d02": 0.365,
    "t4_7e86b677bb2d": 0.189,
    "t4_9a973be6b946": 0.161,
    "t4_28f5ea16adeb": 0.152,
    "t4_4e608398fd6a": 0.125,
    "t4_9265b4bff789": 0.108,
}
#: The two named predictions, verbatim: "t4_da2e98512d02 holds and
#: t4_9265b4bff789 does not."
NAMED = {"t4_da2e98512d02": True, "t4_9265b4bff789": False}

#: The prereg's reading table, transcribed. The verdict is LOOKED UP here.
READINGS = [
    (0.7, 1.01, "BPMD occupancy is the MD-priority filter. Run it on the 300+, "
                "elevate the top by occupancy (~300 GPU-h against ~1,200)"),
    (0.3, 0.7,  "Useful for the extremes only. Elevate the top and bottom decile "
                "to keep testing it; do NOT rank the middle on it"),
    (-0.3, 0.3, "BPMD does not predict residence either. Fall back to tier-1 "
                "drift on all 300+ (~25 GPU-h) and test THAT per molecule"),
    (-1.01, -0.3, "Something is wrong with the protocol or the readout. Report "
                  "as a failure; do NOT reinterpret"),
]


def _mp():
    spec = importlib.util.spec_from_file_location(
        "mdprio_report", REPO / "scripts" / "mdprio_report.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["mdprio_report"] = m
    spec.loader.exec_module(m)
    return m


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    mp = _mp()

    rows = []
    for cand, occ in PREREG.items():
        rep = MD / cand / "md" / "rep1"
        s = mp.series(rep, mp._er())
        r = mp.residence(s)
        if r.get("status") != "ok":
            log.warning("%s: %s", cand, r.get("status"))
            rows.append({"ident": cand, "bpmd_occupancy": occ,
                         "status": r.get("status")})
            continue
        rows.append({"ident": cand, "bpmd_occupancy": occ, "status": "ok",
                     "residence_frac": r["residence_frac"],
                     "dissociated": r["dissociated"],
                     "left_at_ns": r["left_at_ns"],
                     "length_ns": r["length_ns"]})
    t = pd.DataFrame(rows)
    ok = t[t.status == "ok"]
    dest = OUT.write("mdprio_verdict", ".csv")
    t.to_csv(dest, index=False)

    print("\n" + "=" * 72)
    print("  MD-PRIORITY VERDICT — docs/prereg_md_priority.md, as written")
    print("=" * 72 + "\n")
    print(t.sort_values("bpmd_occupancy", ascending=False)
           [["ident", "bpmd_occupancy", "residence_frac", "dissociated",
             "left_at_ns"]].to_string(index=False))

    if len(ok) < len(PREREG):
        print(f"\n  {len(ok)}/{len(PREREG)} complete — the readout needs all six. "
              "Not scoring a partial set.")
        return

    from scipy.stats import spearmanr
    rho, p = spearmanr(ok.bpmd_occupancy, ok.residence_frac)
    print(f"\n  PRIMARY READOUT: rho(BPMD occupancy, 100 ns residence) = {rho:+.3f}"
          f"   (p = {p:.3f}, n = {len(ok)})")

    print("\n  The two NAMED predictions:")
    for cand, expected_hold in NAMED.items():
        row = ok[ok.ident == cand].iloc[0]
        held = not bool(row.dissociated)
        mark = "CORRECT" if held == expected_hold else "WRONG"
        print(f"    {cand}: predicted {'holds' if expected_hold else 'does not hold'}"
              f" -> {'held' if held else f'left at {row.left_at_ns:.0f} ns'}   {mark}")

    verdict = next(txt for lo, hi, txt in READINGS if lo <= rho < hi)
    print(f"\n  VERDICT, looked up in the pre-registered table:\n    {verdict}")
    print("\n  n = 6 supports only large effects. A null here means \"not "
          "demonstrated at n = 6\",\n  never \"shown to be absent\". No p-value "
          "is quoted as evidence.")
    print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
