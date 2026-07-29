"""
Purpose: Merge ensemble MM-GBSA dG and its uncertainty onto each approach's frame.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-29
Input: ensemble_dg_index.jsonl + each approach's latest D-frame
Output: a new D-frame version carrying dG_ensemble, its SEM/SD and provenance

Run:  python scripts/merge_ensemble_dg.py [--dry-run] [--approach t4]

WHY THE FRAME AND NOT THE GUI. D0008 makes the files the source of truth and the
interface a view. Joining the ensemble index inside the Streamlit app would put
a scientific quantity somewhere only the app could see, so the frame carries it
and every consumer -- app, notebook, later analysis -- reads one place.

THE COLUMNS ARE DELIBERATELY NOT CALLED dG_kcal. The single-structure value uses
THREE INDEPENDENT minimisations; the ensemble uses one trajectory. They are
different estimators, not two precisions of the same one (D0036), so they sit in
separate columns and any consumer that wants to compare them has to say which it
means. Overwriting dG_kcal would have silently changed what that column is.
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

log = logging.getLogger("merge-ensemble-dg")

DATA = Path("/data/lab_vm/append_only/inhibition")
INDEX = DATA / "00_shared_substrate" / "ensemble_dg_index.jsonl"
EXPERIMENTS = {
    "t1": "01_t1_de_novo",
    "t2": "02_t2_atra_crem",
    "t3": "03_t3_reinvent",
    "t4": "04_t4_combinatorial",
}
COLS = {
    "dG_kcal_ensemble": "dG_ensemble_kcal",
    "dG_sem_kcal": "dG_ensemble_sem_kcal",
    "dG_sd_kcal": "dG_ensemble_sd_kcal",
    "n_frames": "dG_ensemble_n_frames",
    "statistical_inefficiency": "dG_ensemble_stat_inefficiency",
    "decomposition": "dG_ensemble_decomposition",
}


def load_index() -> pd.DataFrame:
    if not INDEX.is_file():
        raise SystemExit(f"no ensemble index at {INDEX}; run "
                         "scripts/run_mmgbsa_ensemble.py first")
    recs = [json.loads(l) for l in INDEX.read_text().splitlines() if l.strip()]
    df = pd.DataFrame([r for r in recs if "ensemble_error" not in r])
    keep = ["id", "approach", *COLS]
    df = df[[c for c in keep if c in df.columns]].rename(columns=COLS)
    return df.rename(columns={"id": "candidate_id"})


def merge_one(approach: str, ens: pd.DataFrame, dry: bool) -> None:
    experiment = EXPERIMENTS[approach]
    df, path = dio.latest_frame(experiment, approach)
    sub = ens[ens.approach == approach].drop(columns=["approach"])
    if sub.empty:
        log.warning("%s: no ensemble results", approach)
        return

    df = df.drop(columns=[c for c in COLS.values() if c in df.columns])
    merged = df.merge(sub.drop_duplicates("candidate_id"),
                      on="candidate_id", how="left")
    if len(merged) != len(df):
        raise SystemExit(f"{approach}: merge changed row count "
                         f"{len(df)} -> {len(merged)}")

    n = int(merged["dG_ensemble_kcal"].notna().sum())
    both = merged[merged.dG_ensemble_kcal.notna() & merged.dG_kcal.notna()]
    log.info("%s: %s  %d rows carry an ensemble dG (%d also have a "
             "single-structure value)", approach, path.name, n, len(both))
    if not both.empty:
        d = (both.dG_ensemble_kcal - both.dG_kcal)
        log.info("%s:   ensemble - single: mean %+.2f  min %+.2f  max %+.2f "
                 "kcal/mol (different estimators, not an error)",
                 approach, d.mean(), d.min(), d.max())
    if dry:
        return
    out = dio.write_full_frame(
        merged, approach=approach, experiment=experiment,
        stage=f"{approach}_ensemble_dg",
        params={"decision": "D0036",
                "source": INDEX.name,
                "estimator": "single-trajectory 3-leg, sander single points",
                "md": "2 ns GB implicit solvent (GBn2/igb=8), 100 frames",
                "uncertainty": "SEM widened by the statistical inefficiency",
                "distinct_from_dG_kcal": "single-structure dG uses three "
                                         "INDEPENDENT minimisations; these are "
                                         "different estimators (D0036)"},
        inputs={"previous_frame": path, "ensemble_index": INDEX})
    log.info("%s: wrote %s", approach, out.name)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--approach", action="append", choices=sorted(EXPERIMENTS))
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    ens = load_index()
    for a in (args.approach or sorted(EXPERIMENTS)):
        try:
            merge_one(a, ens, args.dry_run)
        except Exception as exc:  # noqa: BLE001
            log.error("%s: %s", a, exc)


if __name__ == "__main__":
    main()
