"""
Purpose: Analyse the explicit-solvent trajectories and merge the result onto frames.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-29
Input: completed GROMACS workdirs under <experiment>/gromacs/<candidate_id>/
Output: a new D-frame per approach carrying the explicit-solvent columns

Run:  /data/lab_vm/envs/dwi_amber_md/bin/python3 scripts/merge_gromacs_results.py \
        [--workers 6] [--dry-run] [--approach t1]

WHY THE FRAME. D0008 makes the files the source of truth and the interface a
view, so a scientific quantity that only the Streamlit app could see would be in
the wrong place. The frame carries it; app, notebook and later analysis all read
one thing.

THE COMPARISON THIS ENABLES, AND ITS ONE TRAP. `explicit_ligand_rmsd_nm_mean`
is directly comparable with the implicit run's `ligand_rmsd_nm_mean`: same
quantity, same units, protein motion removed in both. That comparison is the
point of the whole tier.

`gmx_contacts_mean` is NOT comparable with `mean_contacts`. The implicit metric
counted heavy-atom pairs within 0.45 nm; `gmx mindist -on` uses its own
definition and returns several-fold smaller numbers on the same complex. The
columns are named differently on purpose so nothing can quietly subtract one
from the other.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio  # noqa: E402

log = logging.getLogger("merge-gromacs")

DATA = Path("/data/lab_vm/append_only/inhibition")
EXPERIMENTS = {"t1": "01_t1_de_novo", "t2": "02_t2_atra_crem"}
INDEX = DATA / "00_shared_substrate" / "gromacs_analysis_index.jsonl"

COLS = ("explicit_ligand_rmsd_nm_mean", "explicit_ligand_rmsd_nm_final",
        "explicit_ligand_rmsd_nm_max", "explicit_frac_frames_engaged",
        "gmx_contacts_mean", "ns_analysed", "n_frames_analysed")

# The implicit-solvent counterparts, carried onto the same rows. Without them
# the frame holds half of a comparison, and the comparison IS the result: the
# two solvent models correlate at Spearman -0.102, so a reader seeing only one
# would take a solvent-model artefact for a property of the molecule.
IMPLICIT_INDEX = DATA / "00_shared_substrate" / "md_ensemble_index.jsonl"
IMPLICIT_COLS = ("ligand_rmsd_nm_mean", "ligand_rmsd_nm_max",
                 "frac_frames_engaged", "mean_contacts")


def implicit_residence() -> pd.DataFrame:
    """Implicit-solvent residence per candidate, or an empty frame."""
    if not IMPLICIT_INDEX.is_file():
        log.warning("no implicit index at %s; frames will carry only the "
                    "explicit half of the comparison", IMPLICIT_INDEX)
        return pd.DataFrame()
    recs = [json.loads(l) for l in IMPLICIT_INDEX.read_text().splitlines()
            if l.strip()]
    df = pd.DataFrame([r for r in recs if "md_error" not in r])
    keep = ["candidate_id", *[c for c in IMPLICIT_COLS if c in df.columns]]
    return df[keep].drop_duplicates("candidate_id") if not df.empty else df


def analyse_one(job: dict) -> dict:
    os.nice(19)
    sys.path.insert(0, str(REPO))
    from shared import gromacs_analysis as ga
    wd = Path(job["wd"])
    try:
        return {**job, **ga.analyse(wd)}
    except Exception as exc:  # noqa: BLE001
        (wd / "analysis_traceback.txt").write_text(traceback.format_exc(),
                                                  encoding="utf-8")
        return {**job, "analysis_error": f"{type(exc).__name__}: {str(exc)[:250]}"}


def discover(approaches: list[str]) -> list[dict]:
    """Every trajectory to analyse, INCLUDING the per-replicate ones.

    THE REPLICATE LAYOUT IS NESTED. A candidate directory holds the shared
    solvation and minimisation; each replicate gets its own `repN/` with its own
    `prod.xtc`. This previously looked only at `wd/prod.xtc`, so it would have
    found the 48 original single runs, missed all 240 replicates, and reported
    the old numbers as the new result -- succeeding silently, which is the
    failure mode this project keeps meeting.

    A flat `prod.xtc` is still picked up and labelled `replicate: 0`, so the
    original campaign stays analysable and cannot be mistaken for replicate 1.
    """
    jobs = []
    for a in approaches:
        root = DATA / EXPERIMENTS[a] / "gromacs"
        if not root.is_dir():
            continue
        for wd in sorted(p for p in root.iterdir() if p.is_dir()):
            flat = wd / "prod.xtc"
            if flat.is_file() and flat.stat().st_size:
                jobs.append({"approach": a, "id": wd.name, "wd": str(wd),
                             "replicate": 0})
            for rep in sorted(wd.glob("rep*")):
                xtc = rep / "prod.xtc"
                if not (rep.is_dir() and xtc.is_file() and xtc.stat().st_size):
                    continue
                try:
                    n = int(rep.name[3:])
                except ValueError:
                    continue
                jobs.append({"approach": a, "id": wd.name, "wd": str(rep),
                             "replicate": n})
    return jobs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--approach", action="append", choices=sorted(EXPERIMENTS))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    jobs = discover(args.approach or sorted(EXPERIMENTS))
    if not jobs:
        raise SystemExit("no completed GROMACS trajectory found")
    log.info("analysing %d trajectories on %d workers", len(jobs), args.workers)

    ok, bad = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for fut in as_completed([ex.submit(analyse_one, j) for j in jobs]):
            r = fut.result()
            (bad if "analysis_error" in r else ok).append(r)
            if "analysis_error" in r:
                log.error("%s: %s", r["id"], r["analysis_error"][:120])
    INDEX.write_text("".join(json.dumps(r) + "\n" for r in ok + bad),
                     encoding="utf-8")
    log.info("wrote %s (%d ok, %d failed)", INDEX, len(ok), len(bad))
    if not ok:
        raise SystemExit("no trajectory analysed successfully")

    ens = pd.DataFrame(ok).rename(columns={"id": "candidate_id"})
    for a in sorted({r["approach"] for r in ok}):
        experiment = EXPERIMENTS[a]
        df, path = dio.latest_frame(experiment, a)
        sub = ens[ens.approach == a][["candidate_id",
                                      *[c for c in COLS if c in ens.columns]]]
        df = df.drop(columns=[c for c in (*COLS, *IMPLICIT_COLS)
                              if c in df.columns])
        merged = df.merge(sub.drop_duplicates("candidate_id"),
                          on="candidate_id", how="left")
        imp = implicit_residence()
        if not imp.empty:
            merged = merged.merge(imp, on="candidate_id", how="left")
        if len(merged) != len(df):
            raise SystemExit(f"{a}: merge changed row count")
        n = int(merged["explicit_ligand_rmsd_nm_mean"].notna().sum())
        log.info("%s: %s -> %d rows carry explicit-solvent results",
                 a, path.name, n)

        # The comparison, printed where it is computed rather than left to the
        # interface: does the implicit-solvent verdict survive real water?
        both = merged[merged.explicit_ligand_rmsd_nm_mean.notna()
                      & merged.get("ligand_rmsd_nm_mean",
                                   pd.Series(dtype=float)).notna()] \
            if "ligand_rmsd_nm_mean" in merged.columns else merged.iloc[0:0]
        if not both.empty:
            d = both.explicit_ligand_rmsd_nm_mean - both.ligand_rmsd_nm_mean
            log.info("%s:   explicit - implicit ligand RMSD: mean %+.3f nm "
                     "(min %+.3f, max %+.3f) over %d candidates",
                     a, d.mean(), d.min(), d.max(), len(both))

        if args.dry_run:
            continue
        out = dio.write_full_frame(
            merged, approach=a, experiment=experiment,
            stage=f"{a}_gromacs_explicit",
            params={"solvent": "explicit TIP3P, tleap-solvated, 12 A padding",
                    "engine": "GROMACS 2026.3 CUDA",
                    "rmsd_comparable_with": "ligand_rmsd_nm_mean (implicit tier)",
                    "contacts_NOT_comparable": "gmx mindist -on differs from the "
                                               "implicit tier's heavy-atom pair count",
                    "feeds_no_gate": True},
            inputs={"previous_frame": path, "analysis_index": INDEX})
        log.info("%s: wrote %s", a, out.name)


if __name__ == "__main__":
    main()
