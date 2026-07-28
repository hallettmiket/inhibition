"""
Purpose: T_2 step 8 — T5 physics rescoring (MM-GBSA) on the shortlist.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the latest D2 frame (post ranking, with `shortlist`)
Output: D2 with dG columns; per-candidate Amber working directories

THE STEP T_2's PLAN CALLED T5. `docs/approaches/t2.md` step 8: "Physics
rescoring on survivors... computes a higher-fidelity binding estimate on the few
candidates that warrant it... applies no filter of its own." That last clause is
honoured here: this stage scores and records, it never stamps `rejected_at`.

SHORTLIST ONLY. Each candidate needs AM1-BCC charges plus three minimisations of
a ~2,400-atom system — minutes each, against milliseconds for docking. It runs
on the 25 the ranking selected, not on 1,882.

WHAT IS AND IS NOT DELIVERED HERE. The plan's step 8 lists four things: MM-GBSA,
short explicit-solvent MD, an AI cofold pose with physics relaxation, and
anti-target/selectivity docking. Only the FIRST is implemented. There is no MD
anywhere in this repo, so the plan's "ΔG_bind ± uncertainty" is reported without
the uncertainty rather than with an invented one — a single minimisation has no
error bar. Cofold poses and anti-target scores are absent entirely.

COMPARABLE ACROSS THE WHOLE APPROACH. Unlike T_4's covalent dG, which carries a
constant bond term that cancels only within a warhead class (D0020), this dG has
no such term. T_2 may rank on it globally.

READ IT WITH D0031. The 25 candidates here were chosen by docking, and on
class-matched decoys docking is at chance on this receptor. This is therefore
not confirmation of the docking ranking but an independent estimate of the same
quantity; where the two disagree, that disagreement is the result.

FAILURES ARE STAMPED, NOT DROPPED. A candidate whose antechamber or tleap build
fails keeps its row with `mmgbsa_error` set.
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

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                          # noqa: E402
from shared import mmgbsa as mg                       # noqa: E402
from shared import mmgbsa_noncovalent as mgn          # noqa: E402

log = logging.getLogger("t2-mmgbsa")

EXPERIMENT = "02_t2_atra_crem"
DATA = Path("/data/lab_vm/append_only/inhibition") / EXPERIMENT
WORK_ROOT = DATA / "mmgbsa"
POSE_ROOT = DATA / "docking" / "poses"
NICE = 19


def score_one(job: dict) -> dict:
    """Build, minimise and score one candidate's three legs."""
    os.nice(NICE)
    cid = job["candidate_id"]
    wd = WORK_ROOT / cid
    wd.mkdir(parents=True, exist_ok=True)

    done = wd / "result.json"
    if done.is_file():
        return json.loads(done.read_text())

    try:
        pose = POSE_ROOT / f"{cid}_out.pdbqt"
        if not pose.is_file():
            pose = POSE_ROOT / f"{cid}.pdbqt"
        if not pose.is_file():
            raise mg.MMGBSAError(f"no docked pose for {cid} under {POSE_ROOT}")

        sdf = mgn.pose_to_sdf(pose, wd, job["smiles"])
        # The apo receptor is what a non-covalent complex needs; prepare_receptor
        # emits it alongside the CYX form, so both approaches share one
        # preparation and cannot drift on histidine states.
        _cyx, cys, _idx, _n = mg.prepare_receptor(wd)
        mol2, frcmod = mgn.parameterize_ligand(sdf, wd,
                                               net_charge=int(job["formal_charge"]))
        legs = mgn.build_topologies(wd, mol2, frcmod, cys)
        verified = mgn.verify_complex(legs)

        energies = {leg: mg.minimize_and_score(wd, leg)
                    for leg in ("complex", "receptor", "ligand")}
        result = {"candidate_id": cid, **mgn.delta_g(energies),
                  **{f"topology_{k}": v for k, v in verified.items()}}
        done.write_text(json.dumps(result, indent=2), encoding="utf-8")
        return result
    except Exception as exc:  # noqa: BLE001 - one failure must not end the run
        (wd / "traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return {"candidate_id": cid, "mmgbsa_error": str(exc)[:300]}


def main() -> None:
    ap = argparse.ArgumentParser(description="T_2 step 8: MM-GBSA rescoring.")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=6,
                    help="candidates scored concurrently. One process writes the "
                         "frame at the end, so parallel workers cannot race on "
                         "frame versions the way a tmux fan-out does.")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    WORK_ROOT.mkdir(parents=True, exist_ok=True)

    df, frame_path = dio.latest_frame(EXPERIMENT, "t2")
    if "shortlist" not in df.columns:
        raise SystemExit("frame has no shortlist — run 04_rank.py first")
    if "formal_charge" not in df.columns:
        raise SystemExit("frame has no formal_charge — antechamber needs the "
                         "ligand's net charge and must not be given a guess")

    todo = df[df["shortlist"].fillna(False)].copy()
    if args.limit:
        todo = todo.head(args.limit)
    log.info("MM-GBSA on %d shortlisted candidates from %s (%d workers)",
             len(todo), frame_path.name, args.workers)

    jobs = [{"candidate_id": r["candidate_id"],
             "smiles": r["canonical_smiles"],
             "formal_charge": r.get("formal_charge", 0) or 0}
            for _, r in todo.iterrows()]

    results, failures = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(score_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            (failures if "mmgbsa_error" in r else results).append(r)
            log.info("[%d/%d] %s -> %s", i, len(jobs), r["candidate_id"],
                     f"dG {r['dG_kcal']:.2f}" if "dG_kcal" in r
                     else f"FAILED {r['mmgbsa_error'][:80]}")

    if not results and not failures:
        raise SystemExit("nothing to do")

    cols = pd.DataFrame(results + failures)
    keep = [c for c in ("candidate_id", "dG_kcal", "G_complex", "G_receptor",
                        "G_ligand", "mmgbsa_error") if c in cols.columns]
    stale = [c for c in keep if c != "candidate_id" and c in df.columns]
    if stale:
        df = df.drop(columns=stale)
    merged = df.merge(cols[keep].drop_duplicates("candidate_id"),
                      on="candidate_id", how="left")
    if len(merged) != len(df):
        raise RuntimeError(f"merge changed row count {len(df)} -> {len(merged)}")

    if args.limit:
        log.warning("--limit %d: NOT writing a frame; a partial run must not "
                    "become the latest frame for the next stage to read.",
                    args.limit)
        print(f"\n--limit {args.limit}: scored {len(results)}, "
              f"failed {len(failures)}; no frame written.")
        return

    out = dio.write_full_frame(
        merged, approach="t2", experiment=EXPERIMENT, stage="t2_mmgbsa",
        params={"tier": "T5 physics rescoring (docs/approaches/t2.md step 8)",
                "scheme": "non-covalent 3-leg",
                "igb": mg.IGB, "pb_radii": mg.PB_RADII,
                "ensemble_averaged": False,
                "uncertainty_reported": False,
                "not_implemented": ["explicit-solvent MD", "AI cofold pose",
                                    "anti-target / selectivity docking"],
                "comparable": "across the whole approach (no link-atom constant)",
                "n_scored": len(results), "n_failed": len(failures)},
        inputs={"frame": frame_path})

    print(f"\nT_2 MM-GBSA (T5 physics rescoring) -> {out}")
    print(f"  scored {len(results)}, failed {len(failures)}\n")
    if results:
        r = pd.DataFrame(results).sort_values("dG_kcal")
        v = merged[merged["candidate_id"].isin(r["candidate_id"])]
        rank_by_dock = dict(zip(v["candidate_id"], v["rank"]))
        print(f"  {'candidate':20s} {'dG (kcal/mol)':>13s} {'dock rank':>10s}")
        print("  " + "-" * 46)
        for _, x in r.iterrows():
            print(f"  {x['candidate_id']:20s} {x['dG_kcal']:13.2f} "
                  f"{str(rank_by_dock.get(x['candidate_id'], '-')):>10s}")
        print("\n  dG is comparable across the whole approach — unlike T_4's")
        print("  covalent dG, there is no per-class constant bond term.")
        print("  No uncertainty is reported: a single minimisation has none,")
        print("  and the MD the plan pairs with this does not exist here.")
    for f in failures:
        print(f"  FAILED {f['candidate_id']}: {f['mmgbsa_error'][:120]}")


if __name__ == "__main__":
    main()
