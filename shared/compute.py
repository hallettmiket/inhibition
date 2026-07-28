"""
Purpose: The project's CPU/GPU concurrency budget, in one place.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: none
Output: the constants every parallel stage sizes itself from

WHY A CONSTANT RATHER THAN A DEFAULT PER SCRIPT. Concurrency was previously a
literal in each stage (6 here, 20 there), so raising the budget meant finding
every one of them and the effective limit was whichever script you happened to
run. It lives here now, and a stage that wants a different number has to say so
explicitly.

THE NUMBER IS THE PI's, NOT A GUESS. Set to 50 by @mhallet on 2026-07-28,
raising the lab's standing ≤20-core default for this project. The machine has
224 cores and was measured sitting at 18 in use while a 4-worker job crawled.

CONCURRENCY IS WORKERS, NOT THREADS. Every Amber worker must also be pinned to a
single thread. antechamber's AM1-BCC step shells out to `sqm`, which is
OpenMP/MKL-threaded and takes every core it can see — one `sqm` was measured at
7,543% CPU (~75 cores), and six nominal "workers" consumed 204 of 224 cores
while 63 users were logged in. `SINGLE_THREAD_ENV` is what keeps `--workers N`
meaning N cores.

YIELD WHEN THE MACHINE IS BUSY. This is a shared server. `available_workers()`
reports what is actually free so a long run can size itself down rather than
assuming the budget is always spendable.
"""

from __future__ import annotations

import logging
import os
import subprocess

log = logging.getLogger(__name__)

# The CPU budget for this project's parallel stages.
MAX_CPU_WORKERS = 50

# GPU policy: use up to 6 idle devices, dropping to 4 when other jobs are
# present. gnina occupies only ~500 MiB, which a memory-threshold check cannot
# see — prefer explicit device allocation when another job is already running.
MAX_GPUS_IDLE = 6
MAX_GPUS_SHARED = 4

NICE = 19

SINGLE_THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


def pin_to_one_thread() -> None:
    """Call at the top of every worker process, before any Amber tool runs."""
    os.environ.update(SINGLE_THREAD_ENV)


def available_workers(requested: int | None = None, *, reserve: int = 24) -> int:
    """How many workers to actually start, given what the machine is doing.

    `reserve` leaves headroom for other researchers rather than filling the box.
    A shared server that one job saturates is a job that gets killed.
    """
    want = requested or MAX_CPU_WORKERS
    try:
        total = os.cpu_count() or 1
        used = float(subprocess.check_output(
            "ps -eo pcpu= | awk '{s+=$1} END {print s/100}'",
            shell=True, text=True, timeout=30).strip() or 0)
        free = max(1, int(total - used - reserve))
        n = max(1, min(want, free))
        if n < want:
            log.info("requested %d workers; %.0f of %d cores are in use, so "
                     "starting %d and leaving %d spare", want, used, total, n,
                     reserve)
        return n
    except Exception as exc:  # noqa: BLE001 - never let sizing block a run
        log.warning("could not size the pool (%s); using %d", exc, want)
        return want
