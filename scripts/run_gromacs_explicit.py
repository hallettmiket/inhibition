"""
Purpose: Explicit-solvent GROMACS MD across T_1 and T_2 shortlisted candidates.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-29
Input: the shortlisted candidates' existing MM-GBSA workdirs
Output: one solvated production trajectory per candidate + a combined index

Run:  /data/lab_vm/envs/dwi_amber_md/bin/python3 \
        scripts/run_gromacs_explicit.py --approach t1 --approach t2 \
        [--gpus 0,2,3,4,5,6] [--production-ps 10000] [--limit N]

WHY ONLY T_1 AND T_2. Both are non-covalent, so their ligands can physically
leave the pocket and the water model governs whether they do. T_1 is the sharper
case: under implicit solvent two of its shortlisted candidates dissociated
outright (RMSD 9.0 and 7.3 nm, engaged 0.07 and 0.14), and GB has no water
structure that could have held them. Whether that is real or an artefact of the
solvent model is a question explicit water can answer.

T_3 and T_4 are excluded on purpose, not for cost. Their ligands are covalently
bonded to Cys113 and cannot dissociate, so residence is trivially 1.0 and the
interesting quantity would be energetic -- and their junction dihedral was
corrected only hours ago (D0037), with the cc/cd term the least certain of the
five. Spending GPU-days on systems whose parameters are still settling would
bake that uncertainty into the result.

THIS FEEDS NO GATE. D0036 showed better sampling did not rescue the ranking.
Nothing here reorders a shortlist or enters the enrichment token.
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

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("gromacs-explicit")

DATA = Path("/data/lab_vm/append_only/inhibition")
APPROACHES = {
    "t1": {"experiment": "01_t1_de_novo", "mmgbsa": "mmgbsa_2"},
    "t2": {"experiment": "02_t2_atra_crem", "mmgbsa": "mmgbsa_2"},
}


def discover(approaches: list[str], limit: int | None) -> list[dict]:
    """Shortlisted candidates that have the Amber inputs solvation needs."""
    from shared import io as dio
    jobs = []
    for a in approaches:
        cfg = APPROACHES[a]
        try:
            df, frame = dio.latest_frame(cfg["experiment"], a)
        except Exception as exc:  # noqa: BLE001
            log.warning("%s: no frame (%s)", a, exc)
            continue
        if "shortlist" not in df.columns:
            log.warning("%s: frame %s has no shortlist column", a, frame.name)
            continue
        s = df[df["shortlist"].fillna(False)].drop_duplicates("candidate_id")
        root = DATA / cfg["experiment"] / cfg["mmgbsa"]
        for _, r in s.iterrows():
            src = root / str(r["candidate_id"])
            if not (src / "ligand.mol2").is_file():
                continue
            jobs.append({"approach": a, "id": str(r["candidate_id"]),
                         "src": str(src),
                         "wd": str(DATA / cfg["experiment"] / "gromacs"
                                   / str(r["candidate_id"]))})
        log.info("%s: %d shortlisted candidates with Amber inputs (from %s)",
                 a, sum(1 for j in jobs if j["approach"] == a), frame.name)
    return jobs[:limit] if limit else jobs


def run_one(job: dict) -> dict:
    os.nice(19)
    sys.path.insert(0, str(REPO))
    from shared import gromacs_explicit as gx

    wd = Path(job["wd"])
    done = wd / "gromacs_result.json"
    if done.is_file() and not job.get("force"):
        r = json.loads(done.read_text())
        # Same protocol check the MD and ensemble caches now carry: a shorter
        # cached run is a different run, not a cheaper one.
        if float(r.get("production_ps", 0)) >= float(job["production_ps"]):
            return {**job, **r, "cached": True}
        log.info("%s: cached run is %.0f ps, %.0f requested; recomputing",
                 job["id"], r.get("production_ps", 0), job["production_ps"])
    try:
        r = gx.run_pipeline(Path(job["src"]), wd, gpu_id=job["gpu"],
                            threads=job["threads"],
                            production_ps=job["production_ps"])
        wd.mkdir(parents=True, exist_ok=True)
        done.write_text(json.dumps(r, indent=2), encoding="utf-8")
        return {**job, **r}
    except Exception as exc:  # noqa: BLE001 - one failure must not end the run
        wd.mkdir(parents=True, exist_ok=True)
        (wd / "traceback.txt").write_text(traceback.format_exc(),
                                          encoding="utf-8")
        return {**job, "gromacs_error": f"{type(exc).__name__}: {str(exc)[:300]}"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--approach", action="append", choices=sorted(APPROACHES))
    ap.add_argument("--gpus", default="0,2,3,4,5,6")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--production-ps", type=float, default=10000.0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--out", default=str(DATA / "00_shared_substrate"
                                         / "gromacs_explicit_index.jsonl"))
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    gpus = [int(g) for g in args.gpus.split(",") if g.strip()]
    jobs = discover(args.approach or ["t1", "t2"], args.limit)
    if not jobs:
        raise SystemExit("no shortlisted candidate has the Amber inputs needed")
    for i, j in enumerate(jobs):
        j.update(gpu=gpus[i % len(gpus)], threads=args.threads,
                 production_ps=args.production_ps, force=args.force)

    log.info("%d candidates, %.1f ns each, GPUs %s, %d threads each",
             len(jobs), args.production_ps / 1000.0, gpus, args.threads)

    ok, failed = [], []
    with ProcessPoolExecutor(max_workers=len(gpus)) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for n, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            (failed if "gromacs_error" in r else ok).append(r)
            if "gromacs_error" in r:
                log.error("[%d/%d] %s %s", n, len(jobs), r["id"],
                          r["gromacs_error"][:120])
            else:
                perf = (r.get("stages", {}).get("prod", {}) or {}).get("ns_per_day")
                log.info("[%d/%d] %s %d atoms, %d waters%s", n, len(jobs),
                         r["id"], r.get("atoms", 0), r.get("waters", 0),
                         f", {perf:.0f} ns/day" if perf else "")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        for r in ok + failed:
            fh.write(json.dumps(r) + "\n")
    log.info("wrote %s (%d ok, %d failed)", out, len(ok), len(failed))

    if ok:
        sizes = sorted(r["atoms"] for r in ok if "atoms" in r)
        perfs = [p for r in ok
                 if (p := (r.get("stages", {}).get("prod", {}) or {})
                     .get("ns_per_day"))]
        print(f"\n=== explicit-solvent MD, {len(ok)} candidates ===")
        print(f"  system size: {sizes[0]}-{sizes[-1]} atoms "
              f"(median {sizes[len(sizes)//2]})")
        if perfs:
            print(f"  throughput:  {min(perfs):.0f}-{max(perfs):.0f} ns/day")
        print("  These trajectories feed no gate and reorder no shortlist. "
              "They exist to show whether a docked pose survives real water.")


if __name__ == "__main__":
    main()
