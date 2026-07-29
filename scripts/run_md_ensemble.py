"""
Purpose: Run the GB implicit-solvent MD tier across every candidate that
         already has an MM-GBSA topology, on all idle GPUs.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: complex.prmtop + complex.min.rst under each approach's mmgbsa workdir
Output: per-candidate trajectory (.npy + .dcd) and residence metrics

Run:  /data/lab_vm/envs/dwi_amber_md/bin/python3 scripts/run_md_ensemble.py \
          [--gpus 0,2,3,4,5] [--production-ps 2000] [--limit N] [--approach t3]

WHAT THIS DELIVERS AND WHAT IT DOES NOT. It delivers step 11's pocket-residence
check and the trajectories that ensemble MM-GBSA needs. It does NOT itself
produce an ensemble dG -- that requires the three-leg decomposition to be
redone per frame, and the link-atom scheme caps the receptor and ligand legs
with hydrogens that do not exist in the complex, so the slice is not a slice.
Doing it wrong would silently produce a confident number, which is the failure
mode this project has already hit four times (D0025, D0028, D0029, D0030). The
trajectories are the expensive part and both uses need them, so they are
produced first and the decomposition is built against real frames.

WHY RESIDENCE IS SAFE TO REPORT NOW. The OpenMM-vs-sander check finds a
consistent -2 to -6 kcal/mol offset on absolute energies (~0.08% of a -7300
kcal/mol total), most likely the GBn2 screening parameters OpenMM warns about
for GAFF ligand atoms. That offset would matter for a dG of ~-15 kcal/mol if it
did not cancel between legs -- an open question deferred to the ensemble stage.
It does not touch residence at all, which is measured from geometry.
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

log = logging.getLogger("md-ensemble")

DATA = Path("/data/lab_vm/append_only/inhibition")

# Where each approach's MM-GBSA workdirs live. The directory name differs by
# approach for historical reasons; it is listed rather than guessed.
APPROACH_DIRS = {
    "t1": DATA / "01_t1_de_novo" / "mmgbsa_2",
    "t2": DATA / "02_t2_atra_crem" / "mmgbsa_2",
    "t3": DATA / "03_t3_reinvent" / "mmgbsa",
    "t4": DATA / "04_t4_combinatorial" / "mmgbsa",
    "gate": DATA / "00_shared_substrate" / "mmgbsa_gate",
}


def discover(approaches: list[str]) -> list[dict]:
    """Candidates with BOTH a topology and minimised coordinates.

    A zero-byte complex.prmtop is a tleap failure that left a file behind
    (act_001 is one), so size is checked rather than existence.
    """
    jobs = []
    for name in approaches:
        root = APPROACH_DIRS[name]
        if not root.is_dir():
            log.warning("%s: no such directory %s", name, root)
            continue
        for wd in sorted(p for p in root.iterdir() if p.is_dir()):
            prm = wd / "complex.prmtop"
            rst = wd / "complex.min.rst"
            if not prm.is_file() or prm.stat().st_size == 0:
                continue
            if not rst.is_file():
                continue
            jobs.append({"approach": name, "id": wd.name, "workdir": str(wd)})
    return jobs


def run_one(job: dict) -> dict:
    """One candidate on one GPU. Runs in its own process."""
    os.nice(19)
    sys.path.insert(0, str(REPO))
    from shared import md_ensemble as md

    wd = Path(job["workdir"])
    done = wd / "md" / "md_result.json"
    if done.is_file() and not job.get("force"):
        r = json.loads(done.read_text())
        # A CACHE HIT MUST MATCH THE PROTOCOL, NOT JUST THE PATH. Checking only
        # for the file's existence meant a 40 ps smoke-test trajectory
        # satisfied a request for 2 ns: four candidates were silently skipped
        # by the production run and only surfaced later when the ensemble
        # rescorer refused to average two frames. A shorter cached run is not
        # a cheaper version of the requested one, it is a different one.
        want_ns = round(job["production_ps"] / 1000.0, 3)
        got_ns = round(float(r.get("ns_simulated", 0.0)), 3)
        if got_ns >= want_ns:
            r["cached"] = True
            return {**job, **r}
        log.info("%s: cached run is %.3f ns but %.3f ns was requested; "
                 "recomputing", job["id"], got_ns, want_ns)
    try:
        res = md.run_md(wd, job["id"], device_index=job["gpu"],
                        production_ps=job["production_ps"],
                        equil_ps=job["equil_ps"])
        return {**job, **res.to_dict()}
    except Exception as exc:  # noqa: BLE001 - one failure must not end the run
        (wd / "md_traceback.txt").write_text(traceback.format_exc(),
                                             encoding="utf-8")
        return {**job, "md_error": f"{type(exc).__name__}: {str(exc)[:300]}"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gpus", default="0,2,3,4,5",
                    help="comma-separated device indices to spread across")
    ap.add_argument("--approach", action="append",
                    choices=sorted(APPROACH_DIRS),
                    help="restrict to one or more approaches (default: all)")
    ap.add_argument("--production-ps", type=float, default=2000.0)
    ap.add_argument("--equil-ps", type=float, default=200.0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true",
                    help="recompute candidates that already have md_result.json")
    ap.add_argument("--out", default=str(DATA / "00_shared_substrate"
                                         / "md_ensemble_index.jsonl"))
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    gpus = [int(g) for g in args.gpus.split(",") if g.strip() != ""]
    approaches = args.approach or ["t1", "t2", "t3", "t4", "gate"]
    jobs = discover(approaches)
    if args.limit:
        jobs = jobs[: args.limit]
    if not jobs:
        raise SystemExit("no candidate has both complex.prmtop and "
                         "complex.min.rst; nothing to simulate")

    for i, j in enumerate(jobs):
        j["gpu"] = gpus[i % len(gpus)]
        j["production_ps"] = args.production_ps
        j["equil_ps"] = args.equil_ps
        j["force"] = args.force

    by_app: dict[str, int] = {}
    for j in jobs:
        by_app[j["approach"]] = by_app.get(j["approach"], 0) + 1
    log.info("%d candidates across %s on GPUs %s, %.1f ns each",
             len(jobs), by_app, gpus, args.production_ps / 1000.0)

    results, failed = [], []
    with ProcessPoolExecutor(max_workers=len(gpus)) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for n, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            (failed if "md_error" in r else results).append(r)
            if "md_error" in r:
                log.error("[%d/%d] %s FAILED %s", n, len(jobs), r["id"],
                          r["md_error"][:100])
            else:
                log.info("[%d/%d] %s rmsd %.3f nm  engaged %.2f  %s",
                         n, len(jobs), r["id"],
                         r.get("ligand_rmsd_nm_mean", float("nan")),
                         r.get("frac_frames_engaged", float("nan")),
                         "(cached)" if r.get("cached") else "")

    out = Path(args.out)
    with out.open("w", encoding="utf-8") as fh:
        for r in results + failed:
            fh.write(json.dumps(r) + "\n")
    log.info("wrote %s (%d ok, %d failed)", out, len(results), len(failed))

    # The energy cross-check is reported in aggregate: if OpenMM and sander
    # disagree everywhere, that is a property of the setup, not of a candidate.
    checked = [r for r in results if r.get("valid_checked")]
    if checked:
        deltas = [r["valid_delta_kcal"] for r in checked
                  if "valid_delta_kcal" in r]
        agree = sum(1 for r in checked if r.get("valid_agrees"))
        if deltas:
            print(f"\n=== OpenMM vs sander on the minimised structure ===")
            print(f"  {agree}/{len(checked)} within tolerance; "
                  f"delta min {min(deltas):+.2f} max {max(deltas):+.2f} "
                  f"mean {sum(deltas)/len(deltas):+.2f} kcal/mol")
            print("  A consistent same-signed offset largely cancels in a "
                  "three-leg dG; an inconsistent one does not. That is "
                  "measured at the ensemble stage, not assumed here.")

    if results:
        eng = [r["frac_frames_engaged"] for r in results
               if "frac_frames_engaged" in r]
        rms = [r["ligand_rmsd_nm_mean"] for r in results
               if "ligand_rmsd_nm_mean" in r]
        if eng:
            print(f"\n=== pocket residence ({len(eng)} candidates) ===")
            print(f"  frac_frames_engaged: min {min(eng):.2f} "
                  f"median {sorted(eng)[len(eng)//2]:.2f} max {max(eng):.2f}")
            print(f"  ligand RMSD (nm):    min {min(rms):.3f} "
                  f"median {sorted(rms)[len(rms)//2]:.3f} max {max(rms):.3f}")


if __name__ == "__main__":
    main()
