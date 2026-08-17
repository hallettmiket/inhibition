"""Every directory a RUN owns, derived from one topic.

WHY THIS EXISTS. D0080 made a run's topic a directory rather than a behaviour
flag, and it was applied to docking and ranking and stopped there:

    docking / NAC     <topic>/, <topic>_poses/        topic-scoped
    ranking           rank_v2_<tier>_<topic>_*        topic-scoped
    triage sweep      attack_sweep/                   FLAT
    100 ns            md_residence/                   FLAT
    reports           mdprio_reports/                 FLAT
    workdirs          attack_sweep_10ns/, md_residence_3ikd/   FLAT

So bumping `run.topic` emptied one page of the GUI out of three, and the Sweep
and MD results pages went on showing every screen ever run -- 554 sweep rows and
647 residence rows accumulated across four topics, joined by nothing except
sharing a folder. @tt8804, looking at a freshly bumped run: "the gui is not
updated fresh".

The flat layout also makes a re-run destructive in a way the append-only root is
supposed to prevent: a second screen's sweep CSVs land beside the first's, and
`rank_v2`-style globs concatenate them into one frame with no column saying
which run each row came from.

THE RULE: nothing downstream of docking may name a bare directory. Ask here, and
the answer moves when the topic moves.

LEGACY READS. `mdprio_reports` and friends still exist unscoped on disk with
3.0.0's contents. They are not migrated -- the root is append-only and 2.5 GB of
them is archived browsable under `gui_archive/` -- they simply stop being read.
"""

from __future__ import annotations

from pathlib import Path

from . import target_config as tc

#: Governed outputs (append-only, hook-protected).
BLACKSMITH = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")

#: Scratch: trajectories and MD workdirs. Freely deletable, and deleted between
#: runs -- 1.5 TB of it at the end of 3.0.0.
SCRATCH = Path("/data/lab_vm/modifiable/inhibition")


def topic() -> str:
    return tc.topic()


# -- governed output topics -------------------------------------------------

def sweep_topic(t: str | None = None) -> str:
    """Topic name for the triage sweep's tables."""
    return f"attack_sweep_{t or topic()}"


def residence_topic(t: str | None = None) -> str:
    """Topic name for the production-MD rows."""
    return f"md_residence_{t or topic()}"


def reports_topic(t: str | None = None) -> str:
    """Topic name for the report HTML the GUI browses."""
    return f"mdprio_reports_{t or topic()}"


def bpmd_topic(t: str | None = None) -> str:
    return f"bpmd_{t or topic()}"


def sweep_dir(t: str | None = None) -> Path:
    return BLACKSMITH / sweep_topic(t)


def residence_dir(t: str | None = None) -> Path:
    return BLACKSMITH / residence_topic(t)


def reports_dir(t: str | None = None) -> Path:
    d = BLACKSMITH / reports_topic(t)
    d.mkdir(parents=True, exist_ok=True)
    return d


def bpmd_dir(t: str | None = None) -> Path:
    return BLACKSMITH / bpmd_topic(t)


def poses_dir(t: str | None = None) -> Path:
    """Representative poses -- already topic-scoped before this module."""
    return BLACKSMITH / f"{t or topic()}_poses"


def allposes_dir(t: str | None = None) -> Path:
    return BLACKSMITH / f"{t or topic()}_allposes"


# -- scratch workdir roots --------------------------------------------------

def sweep_work(t: str | None = None) -> Path:
    """Triage-sweep workdirs.

    SEPARATE FROM THE PRODUCTION ROOT, and that is not cosmetic: the workdir is
    <root>/<ident>/md/rep1 regardless of tag, so a sweep and a production run of
    the same mode would collide and the second would find a finished prod.xtc
    sitting there and skip itself.
    """
    return SCRATCH / f"attack_sweep_{t or topic()}"


def residence_work(t: str | None = None) -> Path:
    return SCRATCH / f"md_residence_{t or topic()}"


def bpmd_work(t: str | None = None) -> Path:
    return SCRATCH / f"bpmd_{t or topic()}"
