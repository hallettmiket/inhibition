"""
Purpose: T_4 step 9 — MM-GBSA on the true covalent adduct, for the shortlist.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the latest D4 frame (post within-class ranking, with `shortlist`)
Output: D4 with dG columns; per-candidate Amber working directories

The plan's highest-effort stage. Docking answers "can this ligand sit in the
pocket with its warhead on Cys113"; MM-GBSA answers "once the bond is formed, is
the rest of the molecule actually happy there" — with the covalent bond modelled
explicitly rather than imposed as a docking constraint.

SHORTLIST ONLY. Each candidate needs AM1-BCC charges plus three minimisations of
a 2,400-atom system. That is minutes per candidate against milliseconds for
docking, so it runs on the ~27 candidates the class quota selected, not on 1,683.

WITHIN-CLASS ONLY (D0020, and here for a second, independent reason). dG from
the link-atom scheme carries a constant term for the bond formed and the two
C–H/S–H bonds broken. Constant within a warhead class, different between
classes. Ranking chemotypes against each other on this number would be reading
that constant as chemistry.

FAILURES ARE STAMPED, NOT DROPPED. A candidate whose AM1-BCC or tleap build
fails keeps its row with `mmgbsa_error` set. Silent disappearance from a
shortlist is the failure mode that makes a table look cleaner than the run was.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import traceback
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import covalent_adduct as cad         # noqa: E402
from shared import io as dio                      # noqa: E402
from shared import mmgbsa as mg                   # noqa: E402
from shared import warhead_library as wl          # noqa: E402

log = logging.getLogger("t4-mmgbsa")

EXPERIMENT = "04_t4_combinatorial"
WORK_ROOT = Path("/data/lab_vm/append_only/inhibition") / EXPERIMENT / "mmgbsa"
DOCK_ROOT = Path("/data/lab_vm/append_only/inhibition") / EXPERIMENT / "docking"
NICE = 19


def score_one(row: pd.Series, lib, receptor_cache: dict) -> dict:
    """Build, minimise and score one candidate's three legs."""
    cid = row["candidate_id"]
    wd = WORK_ROOT / cid
    wd.mkdir(parents=True, exist_ok=True)

    done = wd / "result.json"
    if done.is_file():                       # resumable: minutes per candidate
        cached = json.loads(done.read_text())
        # Results written before D0029 carry no dock_id, and the merge is now on
        # dock_id — without this they would rejoin as nulls and silently drop a
        # candidate that had in fact been scored.
        cached.setdefault("dock_id", row["dock_id"])
        # A cached value must also have been produced by the CURRENT scorer.
        # This was the fourth copy of this cache in the repo and the one that
        # kept 11 of 27 T_4 rows on pre-D0033 energies — wrong by up to 28
        # kcal/mol, and enough to invert the chloroacetamide ordering — after
        # the other three had been fixed.
        if mg.cached_result_is_current(cached):
            return cached
        log.info("%s: cached result predates the current scorer; rescoring", cid)

    pose = DOCK_ROOT / f"{row['dock_id']}_docked.sdf"
    if not pose.is_file():
        raise mg.MMGBSAError(f"no docked pose at {pose}")

    smarts = cad.adduct_attachment_smarts(row["warhead_class"], library=lib)
    cyx, cys, cyx_idx, n_res = mg.prepare_receptor(wd)
    mol2, frcmod, att, cap, q = mg.parameterize_ligand(pose, wd, smarts, net_charge=0)
    legs = mg.build_topologies(wd, mol2, frcmod, cyx, cys, cyx_idx, n_res + 1,
                               att, cap, q)
    verified = mg.verify_complex(legs["complex"][0], cyx_idx, att)

    energies = {leg: mg.minimize_and_score(wd, leg)
                for leg in ("complex", "receptor", "ligand")}
    result = {"candidate_id": cid, "dock_id": row["dock_id"],
              "warhead_class": row["warhead_class"],
              **mg.delta_g(energies), **{f"topology_{k}": v
                                         for k, v in verified.items()}}
    done.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description="T_4 step 9: covalent MM-GBSA.")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--no-write", action="store_true",
                    help="compute and cache result.json only; write no frame. "
                         "REQUIRED when running several classes in parallel — "
                         "each worker would otherwise write a full frame holding "
                         "only its own class's dG, and the last to finish would "
                         "become the latest frame. Follow a parallel fan-out "
                         "with one plain run, which merges every cached result.")
    ap.add_argument("--classes", default=None,
                    help="comma-separated warhead classes to restrict to")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    os.nice(NICE)
    WORK_ROOT.mkdir(parents=True, exist_ok=True)

    df, frame_path = dio.latest_frame(EXPERIMENT, "t4")
    if "shortlist" not in df.columns:
        raise SystemExit("frame has no shortlist — run 04_rank_within_class.py first")
    if "dock_id" not in df.columns:
        raise SystemExit("frame has no dock_id — re-run 03_covalent_dock.py (D0022)")

    todo = df[df["shortlist"].fillna(False)].copy()
    # ONE SYSTEM PER MOLECULE, NOT PER ROUTE (D0029). Three warhead classes reach
    # an identical adduct, so the shortlist carries the same molecule up to three
    # times. Scoring each row minimised the SAME 2,400-atom system repeatedly:
    # candidates 1 and 4 of the previous run produced byte-identical inputs and a
    # complex energy equal to the decimal. The result is mapped back to every
    # route afterwards.
    n_rows = len(todo)
    todo = todo.drop_duplicates("dock_id")
    if len(todo) != n_rows:
        log.info("%d shortlisted rows -> %d distinct molecules; %d duplicate "
                 "route(s) will reuse their molecule's result (D0029)",
                 n_rows, len(todo), n_rows - len(todo))
    if args.classes:
        keep = {c.strip() for c in args.classes.split(",")}
        todo = todo[todo["warhead_class"].isin(keep)]
    if args.limit:
        todo = todo.head(args.limit)
    log.info("MM-GBSA on %d shortlisted candidates from %s",
             len(todo), frame_path.name)

    lib = wl.load()
    results, failures = [], []
    for i, (_, row) in enumerate(todo.iterrows(), 1):
        log.info("[%d/%d] %s (%s)", i, len(todo), row["candidate_id"],
                 row["warhead_class"])
        try:
            results.append(score_one(row, lib, {}))
        except Exception as exc:  # noqa: BLE001 - one failure must not end the run
            log.error("  FAILED: %s", str(exc)[:200])
            (WORK_ROOT / row["candidate_id"] / "traceback.txt").write_text(
                traceback.format_exc(), encoding="utf-8")
            failures.append({"candidate_id": row["candidate_id"],
                             "dock_id": row["dock_id"],
                             "warhead_class": row["warhead_class"],
                             "mmgbsa_error": str(exc)[:300]})

    if not results and not failures:
        raise SystemExit("nothing to do")

    cols = pd.DataFrame(results + failures)
    keep_cols = [c for c in ("dock_id", "dG_kcal", "G_complex", "G_receptor",
                             "G_ligand", "mmgbsa_error") if c in cols.columns]
    stale = [c for c in keep_cols if c != "dock_id" and c in df.columns]
    if stale:
        df = df.drop(columns=stale)
    # Merge on dock_id: every synthetic route to a scored molecule gets its dG.
    merged = df.merge(cols[keep_cols].drop_duplicates("dock_id"),
                      on="dock_id", how="left")
    if len(merged) != len(df):
        raise RuntimeError(f"merge changed row count {len(df)} -> {len(merged)}")

    if args.no_write:
        print(f"\n--no-write: cached {len(results)} result(s), {len(failures)} "
              "failure(s); no frame written. Run without --no-write to merge.")
        return

    out = dio.write_full_frame(
        merged, approach="t4", experiment=EXPERIMENT, stage="t4_mmgbsa",
        params={"scheme": "link-atom 3-leg, cut at Cys113 SG-C",
                "igb": mg.IGB, "pb_radii": mg.PB_RADII,
                "ensemble_averaged": False,
                "comparable": "within warhead class only (D0020, D0023)",
                "n_scored": len(results), "n_failed": len(failures)},
        inputs={"d4_frame": frame_path})

    print(f"\nT_4 covalent MM-GBSA -> {out}")
    print(f"  scored {len(results)}, failed {len(failures)}\n")
    if results:
        r = pd.DataFrame(results).sort_values(["warhead_class", "dG_kcal"])
        print(f"  {'warhead class':22s} {'candidate':18s} {'dG (kcal/mol)':>13s}")
        print("  " + "-" * 56)
        for _, x in r.iterrows():
            print(f"  {x['warhead_class']:22s} {x['candidate_id']:18s} "
                  f"{x['dG_kcal']:13.2f}")
        print("\n  dG is comparable WITHIN a warhead class only — the link-atom")
        print("  bond term is constant per class and does not cancel across them.")
    for f in failures:
        print(f"  FAILED {f['candidate_id']} ({f['warhead_class']}): "
              f"{f['mmgbsa_error'][:110]}")


if __name__ == "__main__":
    main()
