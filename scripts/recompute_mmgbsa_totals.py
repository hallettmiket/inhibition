"""
Purpose: Recompute every stored dG with the corrected energy-term list (D0033).
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the existing */min.out files under each approach's mmgbsa workdir
Output: dG_corrected.jsonl per approach + a combined index

Run:  python scripts/recompute_mmgbsa_totals.py [--dry-run]

WHY THIS IS A PARSE, NOT A RERUN. The defect was in reading sander's output,
not in producing it. Every term needed is already sitting in the `.min.out`
files that the original run wrote, so correcting months of dG values costs one
pass over text and no minimisation at all. That is worth stating plainly
because "the energies are wrong" usually implies recomputation, and here it
does not.

WHY THE OLD FILES ARE NOT OVERWRITTEN. `result.json` lives under append_only/,
and the original values are the evidence for what D0032 reported at the time.
Corrected values go to NEW files beside them. A decision record that cites a
number should still be able to find that number.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import mmgbsa  # noqa: E402

log = logging.getLogger("recompute-dg")

DATA = Path("/data/lab_vm/append_only/inhibition")
APPROACH_DIRS = {
    "t1": DATA / "01_t1_de_novo" / "mmgbsa_2",
    "t2": DATA / "02_t2_atra_crem" / "mmgbsa_2",
    "t3": DATA / "03_t3_reinvent" / "mmgbsa",
    "t4": DATA / "04_t4_combinatorial" / "mmgbsa",
    "gate": DATA / "00_shared_substrate" / "mmgbsa_gate",
}


def leg_total(path: Path) -> float:
    """One leg's corrected total, via the same parser the scorer now uses."""
    terms = mmgbsa.parse_energy_block(path.read_text(errors="ignore"))
    return mmgbsa.LegEnergies(terms=terms).total


def recompute(wd: Path) -> dict | None:
    """Corrected dG for one candidate, or None if it was never fully scored."""
    outs = {leg: wd / f"{leg}.min.out" for leg in
            ("complex", "receptor", "ligand")}
    if not all(p.is_file() for p in outs.values()):
        return None
    stored = None
    rj = wd / "result.json"
    if rj.is_file():
        try:
            stored = json.loads(rj.read_text()).get("dG_kcal")
        except Exception:  # noqa: BLE001
            stored = None
    try:
        legs = {leg: leg_total(p) for leg, p in outs.items()}
    except mmgbsa.MMGBSAError as exc:
        return {"id": wd.name, "error": str(exc)[:200]}
    dg = legs["complex"] - legs["receptor"] - legs["ligand"]
    rec = {"id": wd.name,
           "dG_kcal_corrected": round(dg, 3),
           "G_complex": round(legs["complex"], 3),
           "G_receptor": round(legs["receptor"], 3),
           "G_ligand": round(legs["ligand"], 3),
           "energy_terms": list(mmgbsa.ENERGY_TERMS),
           "supersedes_decision": "D0033"}
    if stored is not None:
        rec["dG_kcal_original"] = stored
        rec["shift"] = round(dg - stored, 3)
    return rec


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    combined = []
    for name, root in APPROACH_DIRS.items():
        if not root.is_dir():
            continue
        recs = []
        for wd in sorted(p for p in root.iterdir() if p.is_dir()):
            r = recompute(wd)
            if r is None:
                continue
            r["approach"] = name
            recs.append(r)
            combined.append(r)
        ok = [r for r in recs if "shift" in r]
        if ok:
            shifts = [r["shift"] for r in ok]
            log.info("%-5s %3d candidates  shift mean %+.2f  min %+.2f  "
                     "max %+.2f kcal/mol", name, len(ok),
                     sum(shifts) / len(shifts), min(shifts), max(shifts))
        if not args.dry_run and recs:
            out = root / "dG_corrected.jsonl"
            out.write_text("".join(json.dumps(r) + "\n" for r in recs),
                           encoding="utf-8")

    errs = [r for r in combined if "error" in r]
    log.info("total %d candidates recomputed, %d parse failures",
             len(combined) - len(errs), len(errs))
    for e in errs[:5]:
        log.warning("  %s: %s", e["id"], e["error"][:120])

    if not args.dry_run:
        idx = DATA / "00_shared_substrate" / "dG_corrected_index.jsonl"
        idx.write_text("".join(json.dumps(r) + "\n" for r in combined),
                       encoding="utf-8")
        log.info("wrote %s", idx)


if __name__ == "__main__":
    main()
