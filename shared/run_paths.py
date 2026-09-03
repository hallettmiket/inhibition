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

def _p(key: str, default: str) -> Path:
    """A root from `config/target.yaml`, so it can point at another dataset.

    THE TOOL IS MEANT TO TRAVEL. A literal `/data/lab_vm/.../inhibition` is this
    dataset on this filesystem forever, and 212 of them were spread across the
    scripts. The defaults preserve the current deployment so nothing moves
    today; the point is that they are defaults.
    """
    try:
        return Path(str(tc.get(f"paths.{key}", default=default)))
    except Exception:                                          # noqa: BLE001
        return Path(default)


#: The governed dataset root (append-only, hook-protected).
DATA = _p("governed", "/data/lab_vm/append_only/inhibition")

#: Which agent's outputs: `<DATA>/00_outputs/<agent>/<topic>/`.
AGENT = str(tc.get("paths.agent", default="blacksmith"))

#: Where this agent's topics live.
BLACKSMITH = DATA / "00_outputs" / AGENT

#: Scratch: trajectories and MD workdirs. Freely deletable, and deleted between
#: runs -- 1.5 TB of it at the end of 3.0.0.
SCRATCH = _p("scratch", "/data/lab_vm/modifiable/inhibition")


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


def controls_topic(t: str | None = None) -> str:
    """Topic for the crystal / reactant controls.

    THE LAST UNSCOPED TOPIC (@tt8804: "get rid of the yellow controls, they
    dont belong in this version"). `mdprio_combine` read a flat
    `crystal_controls` directory, so every run showed whatever controls any
    previous run had left there -- including `xtal_6VAJ` and `rx_6VAJ`, built
    against the receptor 3IKD replaced. A control is evidence about a specific
    screen against a specific receptor; carried across runs it is decoration
    that reads as evidence. Scoped like every other topic, so a run shows the
    controls it produced and no others.
    """
    return f"crystal_controls_{t or topic()}"


def controls_dir(t: str | None = None) -> Path:
    return BLACKSMITH / controls_topic(t)


def worklist_topic(t: str | None = None) -> str:
    """Topic for the sweep worklists.

    THE LAST FLAT DIRECTORY, and the one that already cost a live incident: a
    chained script took the newest `sweep_gaps_*.csv` and got the PREVIOUS
    screen's, five days old, then reported its 170 modes as this run's
    selection. `pipeline.worklist_path()` guards it with an mtime check, but a
    guard downstream of a shared directory is a patch, not a fix.
    """
    return f"sweep_gaps_{t or topic()}"


def sweep_dir(t: str | None = None) -> Path:
    return BLACKSMITH / sweep_topic(t)


def sweep_result_files(t: str | None = None) -> list[Path]:
    """This run's sweep tables, OLDEST FIRST, so `keep="last"` means newest.

    ONE RESOLVER, because ten readers had four different answers and two of
    them were wrong:

    * mtime-sorted (right): `sweep_state.results`, `mode_ranking.gather`,
      `mdprio_report`
    * `int(stem.split("_")[-1])` -- **raises** on any stem that is not a bare
      integer. `attack_sweep_21_corrected.csv` (a deliberately superseded row,
      written under the append-only rule because the original could not be
      deleted) crashed `mdprio_combine` outright, so the MD results page
      stopped building and said nothing about why.
    * plain lexicographic `sorted()` -- no crash, but `_10` sorts before `_9`,
      so with `keep="last"` the OLDER measurement of a mode wins. That is the
      defect `mode_ranking`'s own comment warns about.
    * unsorted `glob` -- order is filesystem-dependent, i.e. undefined.

    Sorting on mtime rather than on the version integer is deliberate: it needs
    nothing from the filename, so a corrected, re-scored or otherwise
    non-conforming stem orders correctly instead of crashing or silently
    landing in the wrong place.
    """
    import glob as _glob
    import os as _os
    fs = _glob.glob(str(sweep_dir(t) / "attack_sweep_*.csv"))
    return [Path(f) for f in sorted(fs, key=_os.path.getmtime)]


def residence_dir(t: str | None = None) -> Path:
    return BLACKSMITH / residence_topic(t)


def reports_dir(t: str | None = None) -> Path:
    d = BLACKSMITH / reports_topic(t)
    d.mkdir(parents=True, exist_ok=True)
    return d


def bpmd_dir(t: str | None = None) -> Path:
    return BLACKSMITH / bpmd_topic(t)


def worklist_dir(t: str | None = None) -> Path:
    d = BLACKSMITH / worklist_topic(t)
    d.mkdir(parents=True, exist_ok=True)
    return d


def poses_dir(t: str | None = None) -> Path:
    """Representative poses -- already topic-scoped before this module.

    RESOLVED BY GLOB, NEVER PINNED. `<topic>_poses` is append-only, so a
    correction cannot overwrite it -- it is written beside the original as
    `<topic>_poses_2`, and the original stays because every artefact built
    before the correction was read from it. The highest integer wins, which is
    the same rule `reference_set.latest_reference` follows and for the same
    reason: a literal directory name is a pin, and a pin cannot announce that a
    newer version exists (`how_this_project_breaks.md`, disguise #3).

    WHY THERE IS A `_2` AT ALL (2026-08-29): `write_sdf` stamped each
    representative's `mode` from its WRITE POSITION rather than from the pose,
    so for 128 of 1,684 molecules every representative carried another mode's
    label. The scores were never affected -- they are computed from the cloud --
    but the viewer drew one pose beside another pose's number, and any join on
    `mode` into the representative file paired the wrong rows.
    """
    base = f"{t or topic()}_poses"
    cands = [d for d in BLACKSMITH.glob(f"{base}_*")
             if d.is_dir() and d.name[len(base) + 1:].isdigit()]
    if cands:
        return max(cands, key=lambda d: int(d.name[len(base) + 1:]))
    return BLACKSMITH / base


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


# -- the dataset's own inputs ------------------------------------------------

def frames(tier: str) -> list[Path]:
    """Every version of one tier's candidate library, oldest first.

    ACCEPTS THE TIER OR THE FRAME STEM, and RAISES on anything else. The config
    keys these by tier (`T4`), but three call sites had always spoken in frame
    stems (`D4`) because that is what the filenames say. After the roots moved
    into config, `frames("D4")` matched no key and returned an empty list --
    so the SMILES lookup behind every structure depiction and every rail
    thumbnail silently found nothing, and the panels vanished from the reports
    with no error anywhere. @tt8804: "how many times do I need to ask for the
    viewer to show the structure".

    An unknown key is now a failure, not an empty result: returning [] for a
    plausible-looking name is how a query that cannot be answered gets mistaken
    for a question whose answer is nothing.

    Sorted by the integer version suffix, never lexicographically: `_10` sorts
    before `_9` as a string, and callers take the last element.
    """
    cfg = tc.get("paths.frames", default={}) or {}
    key = str(tier).upper()
    stem = cfg.get(key)
    if stem is None:
        # by frame stem: `D4` -> the entry whose path ends in `/D4`
        for k, v in cfg.items():
            if str(v).rsplit("/", 1)[-1].upper() == key:
                stem, key = v, k
                break
    if stem is None:
        raise KeyError(
            f"no frame library named {tier!r}; config/target.yaml paths.frames "
            f"has {sorted(cfg)} with stems "
            f"{sorted(str(v).rsplit('/', 1)[-1] for v in cfg.values())}")
    fs = list(DATA.glob(f"{stem}_*.parquet"))

    def _v(p: Path) -> int:
        tail = p.stem.rsplit("_", 1)[-1]
        return int(tail) if tail.isdigit() else -1

    return sorted([f for f in fs if _v(f) >= 0], key=_v)


def latest_frame(tier: str) -> Path | None:
    fs = frames(tier)
    return fs[-1] if fs else None


def receptor_prep() -> Path:
    """The prepared receptor, named from `target.pdb`.

    Every script held `receptor_3ikd_prep/3IKD_noligand.pdb` as a literal, which
    is the current TARGET rather than a property of the tool.
    """
    pdb = str(tc.get("target.pdb", default="3IKD"))
    tpl = str(tc.get("paths.receptor_prep",
                     default="receptor_{pdb}_prep/{PDB}_noligand.pdb"))
    return SCRATCH / tpl.format(pdb=pdb.lower(), PDB=pdb.upper())


def receptor_plain() -> Path:
    pdb = str(tc.get("target.pdb", default="3IKD"))
    tpl = str(tc.get("paths.receptor_plain", default="receptor_{pdb}_plain"))
    return SCRATCH / tpl.format(pdb=pdb.lower(), PDB=pdb.upper())


def sidecars() -> Path:
    """SMILES + charge for poses that are not library rows (controls, refs)."""
    d = BLACKSMITH / "pose_sidecars"
    return d
