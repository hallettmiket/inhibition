"""
Purpose: Ensemble MM-GBSA across every candidate with an MD trajectory.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: md/traj_nm.npy + the three prmtops in each MM-GBSA workdir
Output: ensemble/ensemble_dg.json per candidate + a combined index

Run:  /data/lab_vm/envs/dwi_amber_md/bin/python3 scripts/run_mmgbsa_ensemble.py \
          [--workers 50] [--discard-first 10] [--approach t3]

WHAT THIS ANSWERS THAT NOTHING ELSE COULD. Every dG so far is one draw from a
distribution nobody had sampled. The first candidate measured here has a
frame-to-frame standard deviation of 3.26 kcal/mol and a 17 kcal/mol range, so
any ranking that separated candidates by less than a few kcal/mol was reading
noise. This run puts a measured uncertainty on every number.

COST. sander scores 100 frames in ~37 s per leg on one core, so a candidate is
~60 s and the whole set is minutes across 50 workers. The expensive part was
the MD, which already ran.

THREAD PINNING IS NOT OPTIONAL. sander links threaded BLAS; without the pin,
"50 workers" silently becomes the whole machine on a box shared with 60 other
users. See shared/compute.py.
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

from shared import compute  # noqa: E402

log = logging.getLogger("mmgbsa-ensemble")

DATA = Path("/data/lab_vm/append_only/inhibition")
APPROACH_DIRS = {
    "t1": DATA / "01_t1_de_novo" / "mmgbsa_2",
    "t2": DATA / "02_t2_atra_crem" / "mmgbsa_2",
    "t3": DATA / "03_t3_reinvent" / "mmgbsa",
    "t4": DATA / "04_t4_combinatorial" / "mmgbsa",
    "gate": DATA / "00_shared_substrate" / "mmgbsa_gate",
}


def discover(approaches: list[str]) -> list[dict]:
    jobs = []
    for name in approaches:
        root = APPROACH_DIRS[name]
        if not root.is_dir():
            continue
        for wd in sorted(p for p in root.iterdir() if p.is_dir()):
            if (wd / "md" / "traj_nm.npy").is_file():
                jobs.append({"approach": name, "id": wd.name,
                             "workdir": str(wd)})
    return jobs


def run_one(job: dict) -> dict:
    os.nice(compute.NICE)
    compute.pin_to_one_thread()
    sys.path.insert(0, str(REPO))
    from shared import mmgbsa_ensemble as me

    wd = Path(job["workdir"])
    done = wd / "ensemble" / "ensemble_dg.json"
    if done.is_file() and not job.get("force"):
        return {**job, **json.loads(done.read_text()), "cached": True}
    try:
        r = me.ensemble_dg(wd, discard_first=job["discard_first"])
        done.parent.mkdir(parents=True, exist_ok=True)
        done.write_text(json.dumps(r, indent=2), encoding="utf-8")
        return {**job, **r}
    except Exception as exc:  # noqa: BLE001
        (wd / "ensemble_traceback.txt").write_text(
            traceback.format_exc(), encoding="utf-8")
        return {**job, "ensemble_error": f"{type(exc).__name__}: {str(exc)[:300]}"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=compute.MAX_CPU_WORKERS)
    ap.add_argument("--discard-first", type=int, default=10,
                    help="frames dropped as residual equilibration")
    ap.add_argument("--approach", action="append", choices=sorted(APPROACH_DIRS))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--out", default=str(DATA / "00_shared_substrate"
                                         / "ensemble_dg_index.jsonl"))
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    jobs = discover(args.approach or sorted(APPROACH_DIRS))
    if args.limit:
        jobs = jobs[: args.limit]
    if not jobs:
        raise SystemExit("no candidate has md/traj_nm.npy; run the MD tier first")
    for j in jobs:
        j["discard_first"] = args.discard_first
        j["force"] = args.force

    n = compute.available_workers(args.workers)
    log.info("%d candidates with trajectories, %d workers", len(jobs), n)

    ok, failed = [], []
    with ProcessPoolExecutor(max_workers=n) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            (failed if "ensemble_error" in r else ok).append(r)
            if "ensemble_error" in r:
                log.error("[%d/%d] %s %s", i, len(jobs), r["id"],
                          r["ensemble_error"][:110])
            elif i % 10 == 0:
                log.info("[%d/%d] %s dG %.2f +/- %.2f", i, len(jobs), r["id"],
                         r["dG_kcal_ensemble"], r["dG_sem_kcal"])

    out = Path(args.out)
    with out.open("w", encoding="utf-8") as fh:
        for r in ok + failed:
            fh.write(json.dumps(r) + "\n")
    log.info("wrote %s (%d ok, %d failed)", out, len(ok), len(failed))

    if ok:
        sds = sorted(r["dG_sd_kcal"] for r in ok)
        sems = sorted(r["dG_sem_kcal"] for r in ok)
        gs = sorted(r["statistical_inefficiency"] for r in ok)
        mid = len(sds) // 2
        print("\n=== how much of a single-structure dG was noise ===")
        print(f"  frame-to-frame SD (kcal/mol): min {sds[0]:.2f}  "
              f"median {sds[mid]:.2f}  max {sds[-1]:.2f}")
        print(f"  SEM on the mean:              min {sems[0]:.2f}  "
              f"median {sems[mid]:.2f}  max {sems[-1]:.2f}")
        print(f"  statistical inefficiency:     min {gs[0]:.2f}  "
              f"median {gs[mid]:.2f}  max {gs[-1]:.2f}")
        print("  Any ranking that separated candidates by less than the SD was "
              "reading a single draw from this spread.")


if __name__ == "__main__":
    main()
