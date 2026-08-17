#!/usr/bin/env python3
"""
Purpose: run one pipeline stage to completion — the fan-out that used to be shell.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: <stage> (screen | sweep | production | bpmd)
Output: whatever the stage writes; this owns only scheduling and resume

WHY THIS REPLACES THE SHELL WORKERS. Each of these stages fans one script out
over many molecules and several GPUs, and that fan-out was written as bash every
time it was needed. The scripts were never the problem; the bash was:

  * the sweep worker read the worklist POSITIONALLY, taking field 5 as
    `pose_rank` when it holds `global_rank`, so 24 jobs asked for ranks in the
    hundreds. Here the worklist is a DataFrame and columns are named.
  * the 100 ns worker asked for an RMSD column `attack_sweep` does not write,
    and reported "no survivors" for eight hours. Here the survivor set comes
    from `pipeline.survivors()`, which raises rather than returning empty when
    it cannot evaluate.
  * killing a worker left `md_residence_3ikd` orphaned and relaunching `gmx`.
    Here children are tracked and torn down parents-first.

RESUME IS PER ITEM, ALWAYS. Every stage below asks "is this one already done?"
from artefacts on disk before spending a GPU on it, so stopping and restarting
costs only what was in flight -- which is what makes the dashboard's stop button
safe to press.
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import pipeline as pl                     # noqa: E402
from shared import run_paths as rp                    # noqa: E402
from shared import target_config as tc                # noqa: E402

log = logging.getLogger("stage")
PY = pl.PY


def _run(cmd: list[str], timeout: int) -> int:
    """One job, niced. The box is shared, so every job yields."""
    full = ["nice", "-n", "19"] + cmd
    try:
        p = subprocess.run(full, cwd=REPO, capture_output=True, text=True,
                           timeout=timeout)
        if p.returncode != 0:
            log.warning("rc=%d %s", p.returncode,
                        (p.stderr or p.stdout or "")[-300:].replace("\n", " "))
        return p.returncode
    except subprocess.TimeoutExpired:
        log.warning("timeout: %s", " ".join(cmd[-4:]))
        return 124


# ---------------------------------------------------------------------------

def stage_screen(gpus: list[int]) -> None:
    """Dock + NAC over the in-scope molecules, one shard per GPU."""
    scope = pl.scope_idents()
    f = pl.run_dir() / "scope.txt"
    f.write_text("\n".join(scope) + "\n")
    log.info("screen: %d molecules in scope over %d GPUs", len(scope), len(gpus))
    n = len(gpus)
    with ThreadPoolExecutor(max_workers=n) as ex:
        for i, g in enumerate(gpus):
            ex.submit(_run, [PY, str(REPO / "scripts/nac_screen_v2.py"),
                             "--shard", str(i), "--n-shards", str(n),
                             "--gpu", str(g), "--only", str(f), "--no-gnina"],
                      86_400)
    log.info("screen finished")


def stage_sweep(gpus: list[int]) -> None:
    """8 ns triage over the worklist, one job per GPU at a time.

    The worklist is read BY COLUMN NAME. Reading it positionally is what sent
    every job a `global_rank` in place of a `pose_rank`.
    """
    wl = pl.worklist_path()
    if wl is None:
        raise SystemExit("no worklist for this run's ranking")
    d = pd.read_csv(wl)
    need = {"parent_ident", "pose_rank", "ident"}
    if not need <= set(d.columns):
        raise SystemExit(f"worklist lacks {sorted(need - set(d.columns))}")
    jobs = [(str(r.ident), str(r.parent_ident), int(r.pose_rank))
            for r in d.itertuples()]
    log.info("sweep: %d modes over %d GPUs, %.0f ns each",
             len(jobs), len(gpus), tc.md_sweep_ps() / 1000)

    def one(job, gpu):
        ident, parent, prank = job
        if _swept(parent, prank):
            log.info("skip %s (already swept)", ident)
            return
        log.info("gpu%d %s pose_rank %d", gpu, ident, prank)
        # --sweep-ps omitted: it defaults to md.sweep_ps from config.
        _run([PY, str(REPO / "scripts/attack_sweep.py"),
              "--candidates", parent, "--pose-rank", str(prank),
              "--gpu", str(gpu)], 21_600)

    _fan(jobs, gpus, one)
    log.info("sweep finished")


def _swept(parent: str, prank: int) -> bool:
    import glob
    for f in glob.glob(str(rp.sweep_dir() / "*.csv")):
        try:
            d = pd.read_csv(f)
        except Exception:                                      # noqa: BLE001
            continue
        if not {"parent_ident", "pose_rank", "status"} <= set(d.columns):
            continue
        m = d[(d.parent_ident.astype(str) == parent) & (d.pose_rank == prank)
              & d.status.astype(str).str.startswith("ok")]
        if len(m):
            return True
    return False


def stage_production(gpus: list[int], poll: int = 600) -> None:
    """100 ns on sweep survivors, as they qualify.

    A WATCHER, NOT A BATCH: survivors appear across the whole sweep, so this
    polls and starts each as soon as it clears the bar. It exits when the sweep
    is over and the queue is empty.
    """
    poses = rp.poses_dir()
    while True:
        surv = pl.survivors()                    # raises if it cannot evaluate
        md = pl._cat(str(rp.residence_dir() / "*.csv"))
        done = set()
        if not md.empty and {"ident", "production_ps"} <= set(md.columns):
            done = set(md[md.production_ps >= 90_000].ident.astype(str))
        todo = [r for r in surv.itertuples() if str(r.ident) not in done]
        if todo:
            log.info("production: %d survivor(s) queued", len(todo))

            def one(r, gpu):
                mode = str(r.ident).rsplit("_m", 1)[-1]
                log.info("gpu%d 100 ns %s (sweep max %.3f nm)", gpu, r.ident, r.rmsd_max)
                _run([PY, str(REPO / "scripts/md_residence_3ikd.py"),
                      "--candidate", str(r.parent_ident),
                      "--pose", str(poses / f"{r.parent_ident}.sdf"),
                      "--pose-rank", str(int(r.pose_rank)),
                      "--mode", mode, "--gpu", str(gpu), "--keep",
                      "--tag", f"md100_{r.ident}"], 86_400)

            _fan(todo, gpus, one)
        elif not pl.running(pl.BY_NAME["sweep"]):
            log.info("sweep over and queue empty — production done")
            return
        else:
            log.info("no new survivors")
        time.sleep(poll)


def stage_bpmd(gpus: list[int], poll: int = 900) -> None:
    """BPMD on molecules that held for the whole 100 ns."""
    bar = tc.md_survivor_rmsd_nm()
    while True:
        md = pl._cat(str(rp.residence_dir() / "*.csv"))
        todo = []
        if not md.empty and {"production_ps", "explicit_ligand_rmsd_nm_max"} <= set(md.columns):
            held = md[(md.production_ps >= 90_000)
                      & (md.explicit_ligand_rmsd_nm_max < bar)]
            bp = pl._cat(str(rp.bpmd_dir() / "*.csv"))
            seen = set(bp.ident.astype(str)) if not bp.empty and "ident" in bp.columns else set()
            todo = [i for i in held.ident.astype(str).unique() if i not in seen]
        if todo:
            def one(ident, gpu):
                log.info("gpu%d bpmd %s", gpu, ident)
                _run([PY, str(REPO / "scripts/bpmd_run.py"), "--pose", str(ident),
                      "--gpu", str(gpu)], 86_400)
            _fan(todo, gpus, one)
        elif not pl.running(pl.BY_NAME["production"]):
            log.info("production over and queue empty — bpmd done")
            return
        else:
            log.info("no new 100 ns holders")
        time.sleep(poll)


def _fan(items, gpus, fn) -> None:
    """Run `fn(item, gpu)` over items, at most one per GPU concurrently."""
    if not items:
        return
    with ThreadPoolExecutor(max_workers=len(gpus)) as ex:
        futs = []
        for i, item in enumerate(items):
            futs.append(ex.submit(fn, item, gpus[i % len(gpus)]))
        for f in futs:
            f.result()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("stage", choices=("screen", "sweep", "production", "bpmd"))
    ap.add_argument("--gpus", default=None,
                    help="comma-separated GPU ids; defaults to the stage's "
                         "declared budget starting at --gpu0")
    ap.add_argument("--gpu0", type=int, default=0)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    if args.gpus:
        gpus = [int(x) for x in args.gpus.split(",") if x.strip() != ""]
    else:
        n = pl.BY_NAME[args.stage].gpus
        gpus = list(range(args.gpu0, args.gpu0 + max(1, n)))
    log.info("stage %s on GPUs %s (topic %s)", args.stage, gpus, rp.topic())
    {"screen": stage_screen, "sweep": stage_sweep,
     "production": stage_production, "bpmd": stage_bpmd}[args.stage](gpus)


if __name__ == "__main__":
    main()
