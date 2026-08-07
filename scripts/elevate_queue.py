"""
Purpose: launch the elevation suite on the selection queue — one 100 ns trajectory per GPU, in tmux.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: 00_outputs/blacksmith/elevation_queue/queue_<N>.csv
Output: tmux session `elevate100`, one window per molecule; MD under modifiable/

@tt8804, #22: "100 ns non Cov MD is run along with BPMD ... if the 100 ns MD
results show that there is not only 90%+ residence but also majority of time in
attack distance and angle for the warhead then we do covalent MD at 50ns."

WHAT THIS LAUNCHES AND WHAT IT DOES NOT. It starts the 100 ns leg only. The gate
and the covalent leg are deliberately NOT chained behind it: the gate is a
judgement about whether a molecule is worth another 50 ns, and firing it
automatically at 4am on a single replicate would commit compute on a measurement
that carries ~100% relative standard error. The trajectories land; the gate is
applied in the morning against them.

IT STARTS FROM THE POSE SELECTION CHOSE, NOT POSE 1. `select_elevate` walks the
persisted poses in energy order and picks the first in attack geometry, because
the lowest-energy pose frequently is not -- two of the top three T_3 molecules
had a rank-1 pose outside the window, one 8.03 A from the sulfur. The chosen
rank travels in the queue and is passed through, so the trajectory starts where
the chemistry says it should.

ONE REPLICATE PER MOLECULE, AND IT IS A SCREEN NOT A MEASUREMENT. One
dissociation event is one draw from an exponential; the residence estimate that
comes out has roughly 100% relative standard error. Three replicates is the
measurement, and at 4 GPU-hours each that is a 2.5-day queue for 20 molecules --
a decision about throughput that belongs to @tt8804, not to this script.

Fair use: GPUs are taken from --gpus only, one molecule per card, nice -n 19.
"""

from __future__ import annotations

import argparse
import glob
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("elevate-queue")
DATA = Path("/data/lab_vm/append_only/inhibition")
QUEUE = DATA / "00_outputs/blacksmith/elevation_queue"
LOGS = Path("/data/lab_vm/modifiable/inhibition/elevate100_logs")
PY = Path.home() / ".micromamba/envs/dwi_reactive/bin/python"
SESSION = "elevate100"
FORBIDDEN = {0, 4, 7}


def tmux(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["tmux", *args], capture_output=True, text=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--n", type=int, default=4, help="molecules to launch")
    ap.add_argument("--gpus", type=int, nargs="+", default=[1, 2, 3, 5])
    ap.add_argument("--production-ps", type=float, default=100000.0)
    ap.add_argument("--queue", default=None)
    ap.add_argument("--winners", default=None,
                    help="pose_rank_bpmd CSV; the elevated pose is taken from "
                         "its is_winner rows rather than from geometry alone")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    bad = FORBIDDEN & set(args.gpus)
    if bad:
        raise SystemExit(f"GPUs {sorted(bad)} are spoken for; refusing")

    q = Path(args.queue) if args.queue else None
    if q is None:
        fs = sorted(glob.glob(str(QUEUE / "queue_*.csv")),
                    key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
        if not fs:
            raise SystemExit(f"no queue under {QUEUE}")
        q = Path(fs[-1])
    df = pd.read_csv(q)
    log.info("queue %s: %d molecules", q.name, len(df))

    # only molecules whose elevated pose is actually in attack geometry. A
    # trajectory started from a pose 8 A off the sulfur measures nothing about
    # the reaction, so it is not worth 4 GPU-hours.
    ok = df[df.geometry_ok].copy() if "geometry_ok" in df.columns else df.copy()
    skipped = len(df) - len(ok)
    if skipped:
        log.info("%d queued molecules have no reaction-competent pose; not launched",
                 skipped)
    if ok.empty:
        raise SystemExit("nothing in the queue has a viable pose to start from")

    # The pose to elevate comes from the BPMD pose ranking when it exists.
    # Geometry says a pose COULD react; BPMD says it is physically there, and
    # spending 4 GPU-hours on the former when the latter is available is the
    # thing this stage exists to prevent.
    if args.winners:
        w = pd.read_csv(args.winners)
        w = w[w.get("is_winner", False) == True]                 # noqa: E712
        if len(w):
            best = w.set_index("ident").pose_rank.to_dict()
            occ = w.set_index("ident").get("frac_in_window", pd.Series(dtype=float)).to_dict()
            before = ok.pose_rank.copy()
            ok["pose_rank"] = ok.ident.map(best).fillna(ok.pose_rank)
            ok["bpmd_occupancy"] = ok.ident.map(occ)
            moved = int((before != ok.pose_rank).sum())
            log.info("BPMD chose a different pose for %d of %d molecules",
                     moved, len(ok))
            # a molecule with no BPMD winner has not been pose-ranked; say so
            missing = sorted(set(ok.ident) - set(best))
            if missing:
                log.warning("no BPMD winner for %s — elevating the "
                            "geometry-chosen pose", ", ".join(missing[:4]))
        else:
            log.warning("%s has no is_winner rows; using geometry-chosen poses",
                        args.winners)

    ok = ok.sort_values(["class_rank", "warhead_class"]).head(min(args.n, len(args.gpus)))
    LOGS.mkdir(parents=True, exist_ok=True)
    if tmux("has-session", "-t", SESSION).returncode != 0:
        tmux("new-session", "-d", "-s", SESSION, "-n", "idle")

    launched = []
    for (gpu, r) in zip(args.gpus, ok.itertuples()):
        pose_rank = int(r.pose_rank) if r.pose_rank == r.pose_rank else 1
        logf = LOGS / f"{r.ident}.log"
        cmd = (f"cd {REPO} && nice -n 19 {PY} scripts/md_residence_3ikd.py "
               f"--candidate {r.ident} --pose {r.sdf} --pose-rank {pose_rank} "
               f"--production-ps {args.production_ps:.0f} --gpu {gpu} --keep "
               f"--tag elevate 2>&1 | tee {logf}")
        log.info("%s -> GPU %d, pose rank %d (d=%.2f A), log %s",
                 r.ident, gpu, pose_rank, r.pose_distance_A, logf.name)
        if args.dry_run:
            print(f"  would run: {cmd}\n")
            continue
        win = r.ident[-8:]
        tmux("kill-window", "-t", f"{SESSION}:{win}")
        res = tmux("new-window", "-t", SESSION, "-n", win,
                   f"bash -lc '{cmd}; exec bash'")
        if res.returncode != 0:
            log.error("tmux failed for %s: %s", r.ident, res.stderr[:200])
            continue
        launched.append((r.ident, gpu, pose_rank))

    print(f"\n=== launched {len(launched)} × {args.production_ps/1000:.0f} ns ===")
    for i, g, pr in launched:
        print(f"  {i:<22} GPU {g}   pose rank {pr}")
    print(f"\n  watch:  tmux attach -t {SESSION}")
    print(f"  logs:   {LOGS}")
    print("\n  ONE replicate each — a screen, not a residence measurement. The "
          "#22 gate\n  (>=90% residence AND majority in attack geometry) is "
          "applied to these in the\n  morning; the covalent leg is deliberately "
          "not chained behind it unattended.")


if __name__ == "__main__":
    main()
