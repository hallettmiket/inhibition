"""
Purpose: Merge D0033-corrected dG values onto each approach's latest frame.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: dG_corrected_index.jsonl + each approach's latest D-frame
Output: a new D-frame version per approach carrying corrected energies

Run:  python scripts/merge_corrected_dg.py [--dry-run] [--approach t4]

WHAT MOVES AND WHAT DOES NOT. `dG_kcal` becomes the corrected value, because
every consumer (the GUI, the gate) should read the right number without knowing
this happened. The superseded value is preserved alongside as
`dG_kcal_precorrection` rather than discarded -- D0032 was decided on it, and a
decision record that cites a number should be able to find it.

NO SHORTLIST CHANGES. Every approach ranks on a DOCKING metric (`affinity_kcal`
or `vina_affinity`), never on dG, so the correction cannot move shortlist
membership. It moves the energies carried on those rows and the ordering anyone
would impose using them. Stated explicitly because "the energies were wrong"
sounds like it should invalidate the selection, and here it does not.

A SIDE EFFECT WORTH NAMING. T_4's acrylamide systems were scored under a
`--no-write` run, so three legitimate results never reached the frame at all.
They arrive here, taking T_4 from 9 scored candidates to 12. That is a
merge-coverage gain, not a consequence of the energy fix, and the manifest says
so.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio  # noqa: E402

log = logging.getLogger("merge-corrected-dg")

DATA = Path("/data/lab_vm/append_only/inhibition")
INDEX = DATA / "00_shared_substrate" / "dG_corrected_index.jsonl"
EXPERIMENTS = {
    "t1": "01_t1_de_novo",
    "t2": "02_t2_atra_crem",
    "t3": "03_t3_reinvent",
    "t4": "04_t4_combinatorial",
}
ENERGY_COLS = ("dG_kcal", "G_complex", "G_receptor", "G_ligand")


def load_corrected() -> pd.DataFrame:
    if not INDEX.is_file():
        raise SystemExit(f"no corrected index at {INDEX}; run "
                         "scripts/recompute_mmgbsa_totals.py first")
    recs = [json.loads(l) for l in INDEX.read_text().splitlines() if l.strip()]
    df = pd.DataFrame([r for r in recs if "error" not in r])
    return df.rename(columns={"id": "candidate_id",
                              "dG_kcal_corrected": "dG_kcal"})


def merge_one(approach: str, cor: pd.DataFrame, dry: bool) -> None:
    experiment = EXPERIMENTS[approach]
    df, path = dio.latest_frame(experiment, approach)
    sub = cor[cor.approach == approach]
    if sub.empty:
        log.warning("%s: nothing corrected", approach)
        return

    before = int(df["dG_kcal"].notna().sum()) if "dG_kcal" in df else 0

    # Preserve what D0032 was decided on, before overwriting it.
    if "dG_kcal" in df.columns:
        df = df.rename(columns={"dG_kcal": "dG_kcal_precorrection"})
    for c in ("G_complex", "G_receptor", "G_ligand"):
        if c in df.columns:
            df = df.drop(columns=[c])

    keep = ["candidate_id", *ENERGY_COLS]
    add = sub[[c for c in keep if c in sub.columns]].drop_duplicates(
        "candidate_id")
    merged = df.merge(add, on="candidate_id", how="left")
    if len(merged) != len(df):
        raise SystemExit(f"{approach}: merge changed row count "
                         f"{len(df)} -> {len(merged)}")

    after = int(merged["dG_kcal"].notna().sum())
    moved = merged[merged.dG_kcal.notna()
                   & merged.dG_kcal_precorrection.notna()]
    shift = (moved.dG_kcal - moved.dG_kcal_precorrection)
    log.info("%s: %s  scored %d -> %d  shift mean %+.2f (min %+.2f max %+.2f)",
             approach, path.name, before, after,
             shift.mean(), shift.min(), shift.max())

    if dry:
        return
    out = dio.write_full_frame(
        merged, approach=approach, experiment=experiment,
        stage=f"{approach}_dg_correction",
        params={"decision": "D0033",
                "energy_terms": "BOND ANGLE DIHED VDWAALS EEL EGB ESURF "
                                "1-4 VDW 1-4 EEL CMAP",
                "previously_omitted": "1-4 VDW, 1-4 EEL, CMAP",
                "shortlist_unchanged": True,
                "rank_metric_is_docking_not_dG": True,
                "newly_merged": int(after - before)},
        inputs={"previous_frame": path, "corrected_index": INDEX})
    log.info("%s: wrote %s", approach, out.name)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--approach", action="append", choices=sorted(EXPERIMENTS))
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    cor = load_corrected()
    for a in (args.approach or sorted(EXPERIMENTS)):
        try:
            merge_one(a, cor, args.dry_run)
        except Exception as exc:  # noqa: BLE001
            log.error("%s: %s", a, exc)


if __name__ == "__main__":
    main()
