#!/usr/bin/env python3
"""
Purpose: work a sweep worklist for weeks, sizing GPU use to what nobody else is using.
Author: @twu383 (with Claude Code)
Date: 2026-09-02
Input: 00_outputs/blacksmith/sweep_worklist_<topic>/worklist_<N>.csv
Output: rows into attack_sweep_<topic>/ via attack_sweep.py, one per mode; plus
        a heartbeat at modifiable/inhibition/sweep_supervisor/<topic>/status.json

@twu383, 2026-09-02: *"we will not be able to just hog 4 gpus for 12 days ...
ramp up at night if noone is using and ramp down during the days, even pausing
our worklist for a bit if needed ... make sure this happens in a stable and
stealthy manner."*

8,488 modes x ~8 min is ~1,130 GPU-hours. There is no schedule that makes that
polite AND fast, so the design gives up on fast: it takes whatever the box is
not using, and it is correct to run for weeks at one worker or to sit at zero
for a whole afternoon.

WHY A SUPERVISOR AND NOT A BATCH SUBMISSION. There is no scheduler on this box.
The alternative is a fixed fan-out, which is exactly the thing being asked
against -- a run sized for an idle night is antisocial at 10am and a run sized
for 10am wastes the night.

------------------------------------------------------------------------------
THIS PROCESS OUTLIVES THE SESSION THAT STARTED IT, AND CLAUDE.md SAYS WHAT THAT
COSTS. Two `overnight.sh` supervisors once survived 14 and 9 days from a session
in another repo, kept polling, and rebuilt a GUI for a topic that had no data
the moment the topic was bumped. So:

  * the topic is RESOLVED ONCE, at startup, and passed to every child
    explicitly. Bumping `run.topic` mid-run cannot redirect work in flight --
    which is the specific failure that record describes.
  * `status.json` names the PID, the topic, the worklist file and the stop file,
    so whoever finds this can see what it is without reading the code.
  * `STOP` in the state directory halts it cleanly within one poll. Killing a
    worker only makes the supervisor start another; stop the SUPERVISOR.
  * discover it with:  ps -eo pid,etime,cmd | grep [s]weep_supervisor
------------------------------------------------------------------------------

THE STEALTH POLICY, and every number in it is a deliberate choice:

  * GPUs 0, 4 and 7 are NEVER touched, matching `overnight.sh` and
    `elevate_queue.FORBIDDEN`. Candidates are 1, 2, 3, 5, 6.
  * A GPU is available only if NO process belonging to another user is on it.
    Ours are identified by PID ownership, not by name, so another copy of this
    pipeline run by someone else correctly counts as foreign.
  * Concurrency is capped by the clock: quiet hours get more workers. A worker
    already running is never killed to meet a lower cap -- it finishes. Ramping
    down means declining to START work, which is the difference between polite
    and destructive.
  * A CPU-load ceiling, because each sweep takes ~8 threads and the GPU being
    free does not mean the box is. Load is read as a 1-minute average, so a
    transient spike does not stall the queue.
  * One claim per task, by `os.mkdir` -- atomic on POSIX, so two workers (or two
    supervisors) cannot take the same mode. A claim older than STALE_H is
    reclaimed, because a worker killed by a reboot leaves its claim behind.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
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

#: Never ours to take. Matches overnight.sh and elevate_queue.FORBIDDEN.
FORBIDDEN_GPUS = {0, 4, 7}
CANDIDATE_GPUS = [1, 2, 3, 5, 6]

#: Workers and CPU-load ceiling by hour of day, local: (from, to, workers, load).
#:
#: THE BUILT-IN VALUES ARE THE CONSERVATIVE FALLBACK, not the operating policy --
#: `policy.json` in the state directory overrides them and is re-read every
#: poll, so tuning a run that is weeks long does not need a restart. A malformed
#: or missing file falls back to THESE numbers, never to the aggressive ones: a
#: policy file that fails to parse must not quietly turn the campaign loose.
#:
#: Night takes all five candidate cards (@twu383, 2026-09-02: "we can be a bit
#: more aggressive at night"). It still never touches 0, 4 or 7, so three of
#: eight are untouched even at full tilt.
SCHEDULE = {
    "night": (22, 8, 4, 170.0),     # 22:00-07:59
    "evening": (19, 22, 2, 150.0),  # 19:00-21:59
    "day": (8, 19, 1, 130.0),       # 08:00-18:59  -- deliberately ONE
}
#: Fallback ceiling if a window carries none. 224 cores, ~8 threads per sweep.
LOAD_CEILING = 130.0
POLL_S = 120
STALE_H = 6.0


def _now() -> datetime:
    return datetime.now()


def load_policy(state: Path | None = None) -> dict:
    """`policy.json` if it parses, else the built-in SCHEDULE.

    FAILS CLOSED. A truncated or hand-edited file that does not parse gets the
    conservative built-in, and says so -- the alternative is a typo silently
    setting the night cap to something nobody chose, on a process that runs for
    weeks unattended.
    """
    if state is None:
        return dict(SCHEDULE)
    f = state / "policy.json"
    if not f.is_file():
        return dict(SCHEDULE)
    try:
        raw = json.loads(f.read_text())
        out = {}
        for name, v in raw["schedule"].items():
            lo, hi, w, ld = int(v[0]), int(v[1]), int(v[2]), float(v[3])
            if not (0 <= lo <= 23 and 0 <= hi <= 24 and 0 <= w <= len(CANDIDATE_GPUS)
                    and 0 < ld <= (os.cpu_count() or 1)):
                raise ValueError(f"window {name!r} out of range: {v}")
            out[name] = (lo, hi, w, ld)
        if not out:
            raise ValueError("empty schedule")
        return out
    except Exception as exc:                              # noqa: BLE001
        print(f"WARNING policy.json unusable ({exc}); using the built-in "
              f"conservative schedule")
        return dict(SCHEDULE)


def target_workers(when: datetime | None = None,
                   schedule: dict | None = None) -> tuple[int, str, float]:
    """(workers, window name, load ceiling) for the clock."""
    h = (when or _now()).hour
    for name, v in (schedule or SCHEDULE).items():
        lo, hi, n = v[0], v[1], v[2]
        ld = v[3] if len(v) > 3 else LOAD_CEILING
        inside = (lo <= h < hi) if lo < hi else (h >= lo or h < hi)
        if inside:
            return n, name, ld
    return 1, "day", LOAD_CEILING


def gpu_owners() -> dict[int, set[str]]:
    """{gpu index: {usernames with a process on it}} -- ours included.

    By PID OWNERSHIP, never by process name: another user running this same
    pipeline must count as foreign, and a name match would call it ours.
    """
    out: dict[int, set[str]] = {g: set() for g in CANDIDATE_GPUS}
    try:
        q = subprocess.run(
            ["nvidia-smi", "--query-compute-apps=gpu_bus_id,pid",
             "--format=csv,noheader"],
            capture_output=True, text=True, timeout=30)
        idx = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,gpu_bus_id",
             "--format=csv,noheader"],
            capture_output=True, text=True, timeout=30)
    except Exception:                                     # noqa: BLE001
        # A FAILED PROBE MEANS "ASSUME BUSY", never "assume free". An
        # unreadable nvidia-smi that defaulted to free would fan out onto
        # somebody else's cards precisely when we cannot see them.
        return {g: {"__probe_failed__"} for g in CANDIDATE_GPUS}

    bus2idx = {}
    for ln in idx.stdout.strip().splitlines():
        parts = [x.strip() for x in ln.split(",")]
        if len(parts) == 2:
            bus2idx[parts[1]] = int(parts[0])
    for ln in q.stdout.strip().splitlines():
        parts = [x.strip() for x in ln.split(",")]
        if len(parts) != 2:
            continue
        g = bus2idx.get(parts[0])
        if g is None or g not in out:
            continue
        who = subprocess.run(["ps", "-o", "user=", "-p", parts[1]],
                             capture_output=True, text=True).stdout.strip()
        if who:
            out[g].add(who)
    return out


def free_gpus(me: str) -> list[int]:
    """Candidate GPUs with no OTHER user's process on them."""
    owners = gpu_owners()
    return [g for g in CANDIDATE_GPUS
            if g not in FORBIDDEN_GPUS and not (owners[g] - {me})]


def load1() -> float:
    return os.getloadavg()[0]


class Claims:
    """Atomic one-worker-per-task claiming, by mkdir."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)

    def take(self, task_id: str) -> bool:
        d = self.root / task_id
        try:
            os.mkdir(d)                                   # atomic on POSIX
        except FileExistsError:
            # RECLAIM A STALE CLAIM. A worker killed by a reboot leaves its
            # directory behind, and without this the task is lost for good.
            try:
                age_h = (time.time() - d.stat().st_mtime) / 3600.0
            except OSError:
                return False
            if age_h > STALE_H and not (d / "done").exists():
                try:
                    shutil.rmtree(d)
                    os.mkdir(d)
                except OSError:
                    return False
            else:
                return False
        (d / "claimed_at").write_text(datetime.now().isoformat())
        return True

    def finish(self, task_id: str) -> None:
        try:
            (self.root / task_id / "done").write_text(datetime.now().isoformat())
        except OSError:
            pass


def done_tasks(topic: str) -> set[str]:
    """Task ids already on disk, from the sweep table itself.

    Read from the RESULTS, not from the claim directory: a claim says a worker
    started, and only a row says it finished. Recomputing this every poll costs
    a few file reads and means a manual re-run or a hand-added row is honoured.
    """
    out: set[str] = set()
    d = rp.BLACKSMITH / rp.sweep_topic(topic)
    if not d.is_dir():
        return out
    for f in d.glob("*.csv"):
        try:
            df = pd.read_csv(f, usecols=["ident", "status"])
        except Exception:                                 # noqa: BLE001
            continue
        # ONLY A COMPLETED SWEEP COUNTS. `attack_sweep --stage0-only` writes a
        # row with `status = "stage0 only"` and no measurements, and a plain
        # ident match treated those as finished -- 12 modes were marked done by
        # a free geometry probe and would never have been simulated.
        #
        # `ok` and nothing else, matching `nac_screen_v2`'s resume rule and for
        # its stated reason: counting `failed:` rows as done means a transient
        # failure is never retried. A mode that fails REPEATEDLY is held off by
        # its claim directory rather than by this set.
        out.update(df.loc[df.status.astype(str) == "ok", "ident"].astype(str))
    return out


def newest_worklist(topic: str, pinned: Path | None = None) -> Path:
    """The highest-numbered worklist for this topic, or the pinned one.

    RE-RESOLVED EVERY POLL, so the campaign can be re-scoped -- an angle
    filter, a different distance band -- by dropping in `worklist_<N+1>.csv`,
    with no restart. Tasks already finished are skipped because `done` comes
    from the RESULTS table, so a narrower list never re-runs anything and a
    wider one just adds work.
    """
    if pinned is not None:
        return pinned
    fs = sorted((rp.BLACKSMITH / f"sweep_worklist_{topic}").glob("worklist_*.csv"),
                key=lambda f: int(f.stem.split("_")[-1]))
    if not fs:
        raise SystemExit(f"no worklist for topic {topic!r}")
    return fs[-1]


def read_worklist(f: Path, fallback: pd.DataFrame | None) -> pd.DataFrame:
    """Load a worklist, FAILING CLOSED onto the one already in use.

    A half-written CSV -- dropped in while being generated -- must not empty
    the queue or crash a run that is weeks long. Required columns are checked
    rather than assumed, because a file missing `pose_rank` would launch every
    sweep on pose 1 (D0105's defect, arriving by a new route).
    """
    need = {"task_id", "ident", "pose_rank", "priority"}
    try:
        d = pd.read_csv(f)
        missing = need - set(d.columns)
        if missing:
            raise ValueError(f"missing columns {sorted(missing)}")
        if d.empty:
            raise ValueError("empty worklist")
        if d.pose_rank.isna().any():
            raise ValueError("some rows have no pose_rank")
        return d
    except Exception as exc:                              # noqa: BLE001
        print(f"WARNING worklist {f.name} unusable ({exc}); keeping the "
              f"previous list" if fallback is not None else
              f"ERROR worklist {f.name} unusable ({exc}) and no previous list")
        if fallback is None:
            raise
        return fallback


def launch(row, gpu: int, topic: str, logs: Path):
    """One sweep, detached, nice 19. Returns the Popen."""
    logs.mkdir(parents=True, exist_ok=True)
    log = logs / f"{row.task_id}.log"
    cmd = ["nice", "-n", "19", str(PY), "scripts/attack_sweep.py",
           "--candidates", row.ident,
           # THE POSE SET IS NAMED, not defaulted. `attack_sweep` resolves
           # `--pose-dir` from `run.topic` when omitted, so a topic bump
           # mid-campaign would silently start sweeping another run's poses.
           "--pose-dir", str(rp.poses_dir(topic)),
           "--pose-rank", str(int(row.pose_rank)),
           "--sweep-ps", "1200",
           # EARLY GIVE-UP, passed EXPLICITLY rather than left to the script's
           # default. `attack_sweep`'s default is 0 (off) because an abort
           # discards work and must be asked for; the campaign asks for it here,
           # so the policy lives with the campaign and not in a library default
           # that some other caller inherits by accident.
           #
           # 6.0 A is LOOSE and set from few observations: it is past the 4.2 A
           # window and past every pose that has scored above zero so far, so it
           # catches only the unambiguously departed. `start_dist_a` accumulates
           # on every completed row, so this can be re-derived from hundreds of
           # sweeps rather than the handful it rests on now.
           "--abort-above-a", "6.0",
           "--gpu", str(gpu)]
    return subprocess.Popen(cmd, cwd=str(REPO),
                            stdout=log.open("w"), stderr=subprocess.STDOUT,
                            start_new_session=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=None,
                    help="resolved ONCE here; children are told explicitly")
    ap.add_argument("--worklist", default=None)
    ap.add_argument("--max-workers", type=int, default=None,
                    help="hard ceiling on top of the schedule")
    ap.add_argument("--poll", type=int, default=POLL_S)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # NICE OURSELVES, not just the children. The workers were launched with
    # `nice -n 19` from the start, but this process sat at 0 -- and its
    # per-poll work GROWS with the campaign: `done_tasks` re-reads every CSV in
    # the sweep directory, which is a handful today and thousands of files in
    # three weeks. A supervisor that starts negligible and ends up scanning a
    # large directory every two minutes at priority 0 is exactly the kind of
    # slow creep nobody notices.
    #
    # It also makes the children's nice unconditional: `nice -n 19` is RELATIVE
    # to the parent, so this guarantees 19 even if someone launches the
    # supervisor from an already-niced shell.
    try:
        os.nice(19)
    except OSError:
        pass                                              # already at the floor
    topic = args.topic or rp.topic()
    pinned = Path(args.worklist) if args.worklist else None
    wl = newest_worklist(topic, pinned)
    tasks = read_worklist(wl, None)

    state = Path("/data/lab_vm/modifiable/inhibition/sweep_supervisor") / topic
    state.mkdir(parents=True, exist_ok=True)
    stop = state / "STOP"
    claims = Claims(state / "claims")
    logs = state / "logs"
    me = subprocess.run(["id", "-un"], capture_output=True, text=True).stdout.strip()

    def status(**kw):
        (state / "status.json").write_text(json.dumps({
            "pid": os.getpid(), "topic": topic, "worklist": str(wl),
            "stop_file": str(stop), "state_dir": str(state),
            "updated": datetime.now().isoformat(),
            "n_tasks": int(len(tasks)), "user": me,
            "worklist_live": "newest worklist_<N>.csv is re-read every poll",
            "discover": "ps -eo pid,etime,cmd | grep [s]weep_supervisor",
            **kw}, indent=2))

    print(f"supervisor pid {os.getpid()}  topic {topic}")
    print(f"  worklist {wl.name}: {len(tasks):,} tasks "
          f"({int((tasks.priority == 0).sum())} front-loaded)")
    print(f"  state    {state}")
    print(f"  STOP     touch {stop}")
    if args.dry_run:
        n, win, ld = target_workers(schedule=load_policy(state))
        status(state="dry-run", window=win, target=n, load_ceiling=ld,
               free_gpus=free_gpus(me))
        print(f"  window {win}: target {n} workers, load ceiling {ld}; "
              f"free GPUs {free_gpus(me)}")
        return

    running: dict[str, tuple[subprocess.Popen, int, float]] = {}
    stopping = False

    def on_term(*_):
        nonlocal stopping
        stopping = True
    signal.signal(signal.SIGTERM, on_term)
    signal.signal(signal.SIGINT, on_term)

    while True:
        # reap
        for tid in list(running):
            proc, gpu, t0 = running[tid]
            if proc.poll() is not None:
                claims.finish(tid)
                print(f"{datetime.now():%H:%M} done {tid} rc={proc.returncode} "
                      f"{(time.time()-t0)/60:.1f} min gpu{gpu}")
                del running[tid]

        if stop.exists() or stopping:
            status(state="stopping", running=list(running))
            if not running:
                print("stopped cleanly")
                status(state="stopped", running=[])
                return
            time.sleep(10)
            continue

        # RE-READ EVERY POLL, so the policy can be tuned without a restart.
        n, window, ceiling = target_workers(schedule=load_policy(state))
        if args.max_workers is not None:
            n = min(n, args.max_workers)
        ld = load1()
        gpus = free_gpus(me)
        # A RUNNING WORKER IS NEVER KILLED to meet a lower cap. Ramping down
        # means declining to start, which is the whole difference between
        # polite and destructive.
        if ld > ceiling:
            n = 0
        n = min(n, len(gpus))

        # RE-RESOLVE THE WORKLIST TOO. Dropping in worklist_<N+1>.csv re-scopes
        # the campaign live; nothing finished is re-run because `done` comes
        # from the results table rather than from the list.
        newest = newest_worklist(topic, pinned)
        if newest != wl:
            fresh = read_worklist(newest, tasks)
            if fresh is not tasks:
                print(f"{datetime.now():%H:%M} worklist -> {newest.name} "
                      f"({len(fresh):,} tasks, was {len(tasks):,})")
                tasks, wl = fresh, newest
        done = done_tasks(topic)
        todo = tasks[~tasks.task_id.isin(done) & ~tasks.task_id.isin(running)]
        status(state="running", window=window, target=n, load1=round(ld, 1),
               load_ceiling=ceiling, free_gpus=gpus, running=list(running),
               n_done=len(done & set(tasks.task_id)), n_todo=int(len(todo)))

        if todo.empty and not running:
            print("worklist complete")
            status(state="complete", running=[])
            return

        busy = {g for _, g, _ in running.values()}
        for _, row in todo.iterrows():
            if len(running) >= n:
                break
            spare = [g for g in gpus if g not in busy]
            if not spare:
                break
            if not claims.take(row.task_id):
                continue
            g = spare[0]
            running[row.task_id] = (launch(row, g, topic, logs), g, time.time())
            busy.add(g)
            print(f"{datetime.now():%H:%M} start {row.task_id} gpu{g} "
                  f"({window}, target {n}, load {ld:.0f})")

        time.sleep(args.poll)


if __name__ == "__main__":
    main()
