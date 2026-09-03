#!/usr/bin/env python3
"""
Purpose: keep the SWEEP page current as results land, and rebuild nothing else.
Author: @twu383 (with Claude Code)
Date: 2026-09-02
Input: the run's sweep rows + its worklist
Output: sweep_assets/, sweep_pages/, sweep.html — rebuilt only when a new sweep finishes

@twu383, 2026-09-02: *"build the sweep page update as we go not the ranking
page"*.

The 1.2 ns campaign delivers a result every few minutes for ~11 days. Rebuilding
by hand means the page is stale the moment you stop typing; rebuilding on a timer
regardless of whether anything changed burns CPU on a shared box for nothing.
This watches the completed-sweep COUNT and rebuilds only when it moves.

WHAT IT DELIBERATELY DOES NOT TOUCH
-----------------------------------
**The ranking page.** `ranking_page.py` renders 132,027 modes and their pose
assets; it takes minutes and its inputs do not change while the sweep runs --
the screen is finished, the ranking is fixed, and the only thing moving is which
modes have a result. Rebuilding it on every sweep would be minutes of CPU to
change a handful of badges. Run it by hand when the ranking itself changes.

ORDER MATTERS, AND GETTING IT WRONG IS INVISIBLE
------------------------------------------------
1. `sweep_assets`  — movies and plots for the new results
2. `build_gui`     — the counts and `sweep_state.json`, WITH the worklist
3. `sweep_combine` — the reports, LAST

`build_gui` and `sweep_combine` both write `sweep.html`. Running them the other
way round is what produced "No sweep has finished yet" on a run with 17
finished sweeps: `build_gui`, given no worklist, rendered the empty state over
the top of `sweep_combine`'s reports. And building the pages before the assets
is what left 16 of 17 reports without their movie. Both look like a rendering
bug and are an ordering bug.

`build_gui` is given `--worklist` explicitly. Without it `sweep_state.json`
carries `worklist: null`, the page cannot say how many modes are queued, and the
pending count falls back to a number from another campaign.

THIS OUTLIVES THE SESSION THAT STARTED IT, like the supervisor beside it:

    ps -eo pid,etime,cmd | grep [s]weep_gui_refresh
    touch <state>/STOP_GUI       # stops it within one poll
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import run_paths as rp                        # noqa: E402

PY = Path.home() / ".micromamba/envs/dwi_reactive/bin/python"
log = logging.getLogger("sweep-gui")


def n_finished(topic: str) -> int:
    """Completed sweeps for this topic. The trigger, and only `ok` counts.

    A `stage0 only` placeholder or an `aborted` give-up does not change what the
    reports show, so rebuilding on them would be work for no visible difference.
    """
    # DISTINCT MODES, not rows. Summing per file counted a mode once per file it
    # appears in, and a mode legitimately appears in several: a re-scored table
    # (`recompute_attack_ready`) rewrites every finished row into a new
    # versioned file, which took the count from 98 to 196 with nothing new run.
    # Only a trigger, so it cost nothing here -- but it is a completion count
    # that can read as progress, and it was wrong by a factor of two.
    seen = set()
    d = rp.BLACKSMITH / rp.sweep_topic(topic)
    for f in d.glob("*.csv"):
        try:
            x = pd.read_csv(f, usecols=["ident", "status"])
        except Exception:                                 # noqa: BLE001
            continue
        seen |= set(x.loc[x.status.astype(str) == "ok", "ident"].astype(str))
    return len(seen)


def newest_worklist(topic: str) -> Path | None:
    fs = sorted((rp.BLACKSMITH / f"sweep_worklist_{topic}").glob("worklist_*.csv"),
                key=lambda f: int(f.stem.split("_")[-1]))
    return fs[-1] if fs else None


def run(step: str, args: list[str], timeout: int) -> bool:
    t0 = time.time()
    r = subprocess.run([str(PY)] + args, cwd=str(REPO),
                       capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        log.warning("%s failed rc=%d: %s", step, r.returncode,
                    (r.stderr or "").strip()[-200:])
        return False
    log.info("%s ok (%.0fs)", step, time.time() - t0)
    return True


def _results_only_worklist(topic: str, wl: Path, state: Path) -> Path:
    """The worklist narrowed to modes that HAVE a result.

    `sweep_assets` walks every row it is given and asks the filesystem whether a
    finished trajectory exists. Handed the full 4,296-row campaign worklist it
    does that 4,296 times to build assets for the one mode that just finished --
    on a rebuild that fires every few minutes for a week.

    Narrowing is safe because an asset can only be built from a finished run:
    a row with no result contributes nothing but a stat() call. The full list is
    still what `build_gui` and `sweep_combine` see, so the QUEUE is unaffected --
    only the asset builder's input shrinks.
    """
    try:
        d = pd.concat([pd.read_csv(f) for f in
                       (rp.BLACKSMITH / rp.sweep_topic(topic)).glob("*.csv")],
                      ignore_index=True)
        done = set(d.loc[d.status.astype(str) == "ok", "ident"].astype(str))
        w = pd.read_csv(wl)
        key = "task_id" if "task_id" in w.columns else "ident"
        sub = w[w[key].astype(str).isin(done)]
        if sub.empty:
            return wl
        dest = state / "worklist_with_results.csv"
        sub.to_csv(dest, index=False)
        log.info("assets over %d finished mode(s), not the whole %d-row list",
                 len(sub), len(w))
        return dest
    except Exception as exc:                              # noqa: BLE001
        log.warning("could not narrow the worklist (%s); using the full one", exc)
        return wl


def _finished_idents(topic: str) -> list[str]:
    out = []
    for f in (rp.BLACKSMITH / rp.sweep_topic(topic)).glob("*.csv"):
        try:
            x = pd.read_csv(f, usecols=["ident", "status"])
        except Exception:                                 # noqa: BLE001
            continue
        out += list(x.loc[x.status.astype(str) == "ok", "ident"].astype(str))
    return sorted(set(out))


def _build_reports(topic: str) -> int:
    """One `sweep_report` per finished mode whose page is missing or stale.

    THE STEP THAT WAS MISSING, and its absence is why no movie ever appeared.
    `sweep_combine` only LINKS `sweep_pages/<ident>.html`; it builds nothing and
    silently skips any ident whose page is absent. The pages that existed were
    40 KB stubs from an earlier stage, and because they existed they were never
    replaced. A real report is ~9.5 MB -- 99 trajectory frames and the RMSD
    plot -- so the size difference alone said the page was empty.

    Built only when missing or OLDER THAN ITS ASSETS: each page is ~9.5 MB and
    takes seconds, so rebuilding all of them on every new result would be
    minutes of work and a GB of writes to add one report.
    """
    pages = rp.reports_dir(topic) / "sweep_pages"
    assets = rp.reports_dir(topic) / "sweep_assets"
    built = 0
    for ident in _finished_idents(topic):
        pg = pages / f"{ident}.html"
        src = _asset_mtime(topic, ident)
        fresh = (pg.is_file() and src is not None
                 and pg.stat().st_mtime >= src
                 and pg.stat().st_size > 1_000_000)     # a stub is ~40 KB
        if fresh:
            continue
        if run(f"sweep_report {ident}",
               ["scripts/sweep_report.py", "--ident", ident], 600):
            built += 1
    return built


def _asset_mtime(topic: str, ident: str) -> float | None:
    """Newest mtime across EVERY asset a per-mode page embeds.

    THE PAGE EMBEDS BOTH, so staleness has to consider both. This compared the
    page against `<ident>.pdb` alone, and the RMSD/distance figure is a
    base64-embedded `<ident>.png` -- so regenerating every figure with a
    corrected green zone (D0111's 3.5 A band, not the screen's 4.2 A window)
    changed no `.pdb`, and 114 pages would have gone on serving the old band
    with nothing reporting anything wrong. A staleness check scoped to one of
    two inputs is a guard that cannot fail for the other.
    """
    assets = rp.reports_dir(topic) / "sweep_assets"
    ts = [f.stat().st_mtime for f in
          (assets / f"{ident}.pdb", assets / f"{ident}.png") if f.is_file()]
    return max(ts) if ts else None


def n_stale(topic: str) -> int:
    """Finished modes whose per-mode page is missing or older than its asset.

    THE COUNT IS NOT THE ONLY REASON TO REBUILD. The loop watched `n_finished`
    alone, so a change that invalidates every page but adds no result triggered
    nothing -- which is what happened when the movie's warhead atom was fixed
    (D0112): 105 pages held a readout measured from the wrong atom and the
    refresher had no reason to notice. Touching the assets marks them stale;
    this is what makes stale mean something.
    """
    pages = rp.reports_dir(topic) / "sweep_pages"
    n = 0
    for ident in _finished_idents(topic):
        src = _asset_mtime(topic, ident)
        if src is None:
            continue
        pg = pages / f"{ident}.html"
        if (not pg.is_file() or pg.stat().st_mtime < src
                or pg.stat().st_size <= 1_000_000):
            n += 1
    return n


def repo_mid_operation() -> str | None:
    """Is the repo mid-rebase / merge / bisect? Returns the operation, else None.

    THE REFRESHER READS CODE FROM A LIVE WORKING TREE. During a rebase that tree
    is at an INTERMEDIATE commit -- older code, briefly, with no warning -- and a
    rebuild that lands in that window regenerates every page from it. Measured
    2026-09-02: a `git pull --rebase` mid-session had the refresher rebuild
    `sweep.html` from the pre-change `sweep_combine`, and the GUI went to
    "No sweep has finished yet · 0 modes" on a run with 139 finished sweeps.

    Nothing was lost -- the next rebuild after the rebase restored it -- but the
    page said the campaign had produced nothing, which is the worst thing this
    page can say incorrectly. Skipping the rebuild leaves the LAST GOOD page in
    place, which is always a better answer than one built from code that is
    halfway between two commits.
    """
    g = REPO / ".git"
    for name, op in (("rebase-merge", "rebase"), ("rebase-apply", "rebase"),
                     ("MERGE_HEAD", "merge"), ("CHERRY_PICK_HEAD", "cherry-pick"),
                     ("BISECT_LOG", "bisect")):
        if (g / name).exists():
            return op
    return None


def rebuild(topic: str, wl: Path, state: Path) -> None:
    """Assets, reports, counts, then the index. Order is load-bearing.

    `build_gui` and `sweep_combine` both write `sweep.html`; `sweep_combine`
    runs LAST so the reports are what survives. And the per-mode reports are
    built BEFORE `sweep_combine`, because it only links pages that already
    exist.
    """
    run("sweep_assets", ["scripts/sweep_assets.py", "--worklist",
                         str(_results_only_worklist(topic, wl, state))], 3600)
    nb = _build_reports(topic)
    if nb:
        log.info("built %d per-mode report(s)", nb)
    run("build_gui", ["scripts/build_gui.py", "--worklist", str(wl)], 900)
    run("sweep_combine", ["scripts/sweep_combine.py", "--worklist", str(wl)], 900)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=None,
                    help="resolved ONCE at startup, like sweep_supervisor")
    ap.add_argument("--poll", type=int, default=180)
    ap.add_argument("--once", action="store_true", help="rebuild and exit")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    try:
        os.nice(19)
    except OSError:
        pass

    topic = args.topic or rp.topic()
    state = Path("/data/lab_vm/modifiable/inhibition/sweep_supervisor") / topic
    state.mkdir(parents=True, exist_ok=True)
    stop = state / "STOP_GUI"

    # ONE REBUILD AT A TIME. `build_gui` and `sweep_combine` write the same
    # file, so two overlapping rebuilds -- two refreshers, or a refresher and a
    # hand-run rebuild -- can finish in an order that leaves `build_gui`'s
    # "awaiting stage" placeholder on top of real reports. That is exactly what
    # was on screen while the reports sat complete on disk.
    lock = state / "refresh.lock"
    try:
        fd = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, str(os.getpid()).encode()); os.close(fd)
    except FileExistsError:
        try:
            other = int(lock.read_text().strip())
            os.kill(other, 0)
            raise SystemExit(f"another refresher (pid {other}) holds {lock}")
        except (ValueError, ProcessLookupError):
            log.warning("stale lock from a dead process; taking it")
            lock.write_text(str(os.getpid()))
    import atexit
    atexit.register(lambda: lock.unlink(missing_ok=True))

    wl = newest_worklist(topic)
    if wl is None:
        raise SystemExit(f"no worklist for topic {topic!r}")
    log.info("topic %s, worklist %s, poll %ds", topic, wl.name, args.poll)
    log.info("STOP: touch %s", stop)

    if args.once:
        rebuild(topic, wl, state)
        return

    last = -1
    while True:
        if stop.exists():
            log.info("STOP seen; exiting")
            (state / "gui_status.json").write_text(json.dumps(
                {"state": "stopped", "updated": datetime.now().isoformat()}))
            return
        # THE WORKLIST IS RE-RESOLVED TOO. It can be re-scoped mid-campaign
        # (an angle filter, a different band), and rebuilding against the old
        # one would report a queue that no longer exists.
        wl = newest_worklist(topic) or wl
        n = n_finished(topic)
        stale = n_stale(topic)
        op = repo_mid_operation()
        if op:
            # Do NOT touch `last`: the rebuild still owes to be done once the
            # tree settles, and clearing it here would skip it for good.
            log.warning("repo is mid-%s — skipping this rebuild; the current "
                        "page stays up rather than being rebuilt from a "
                        "half-applied tree", op)
        elif n != last or stale:
            why = (f"{n} finished sweeps (was {'?' if last < 0 else last})"
                   if n != last else f"{stale} page(s) stale")
            log.info("%s — rebuilding", why)
            rebuild(topic, wl, state)
            last = n
        (state / "gui_status.json").write_text(json.dumps({
            "state": "watching", "pid": os.getpid(), "topic": topic,
            "worklist": wl.name, "n_finished": n, "n_stale": stale,
            "updated": datetime.now().isoformat(),
            "stop_file": str(stop),
            "note": "rebuilds sweep assets/pages only; never the ranking page",
        }, indent=2))
        time.sleep(args.poll)


if __name__ == "__main__":
    main()
