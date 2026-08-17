"""A run must be standalone: nothing may reach a previous run's data.

WHY THIS FILE EXISTS. @tt8804, after deleting 1.5 TB to start clean:

    "I am very annoyed that things keep carrying over from previous runs when I
     have gone through the trouble of deleting past runs to start fresh. stop
     contaminating things with shoddy builds."

The complaint is exact. `run.topic` was bumped, the data root was cleared, and a
fresh screen was launched -- and the GUI still showed the previous run, five
separate times, each found by the user rather than by a test:

  1. the report SERVER was started with a literal path and served the old topic
  2. `sweep_state` globbed the unscoped `attack_sweep/`
  3. `mode_ranking._step_counts` read the unscoped `sweep_state.json` -> "447 ok"
  4. `mode_ranking.gather` joined status from unscoped `attack_sweep/` AND
     `md_residence/`, so every badge on the ranking page described 3.0.0
  5. a chained script took the newest `sweep_gaps_*.csv` and got 3.0.0's

Each was fixed on discovery, which is the wrong process: five instances of one
defect means the defect is structural. This test is the structural fix. It reads
the SOURCE of every module on the pipeline path and fails if any of them names a
shared run directory, a topic, or a dataset root.

It is deliberately a source scan rather than a behavioural test. The failure mode
is "reads the wrong directory and finds real data there", which cannot be
detected by running the code -- both paths exist and both return rows.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

#: Every module the pipeline invokes or imports to produce a run. If a stage is
#: added to `shared/pipeline.py` it belongs here too.
PIPELINE_PATH = [
    "scripts/nac_screen_v2.py", "scripts/rank_v2.py",
    "scripts/sweep_gap_worklist.py", "scripts/attack_sweep.py",
    "scripts/md_residence_3ikd.py", "scripts/bpmd_run.py",
    "scripts/pipeline.py", "scripts/pipeline_stage.py",
    "scripts/build_gui.py", "scripts/ranking_page.py",
    "scripts/sweep_combine.py", "scripts/sweep_report.py",
    "scripts/sweep_assets.py", "scripts/mdprio_combine.py",
    "scripts/mdprio_report.py", "scripts/promote_to_bpmd.py",
    "scripts/serve_reports.py",
    "shared/pipeline.py", "shared/run_paths.py", "shared/mode_ranking.py",
    "shared/sweep_state.py", "shared/mode_assets.py",
]

#: READERS, not producers (#74). These do not create a run; they display one.
#: They are held to the CROSS-RUN CONTAMINATION rule only -- no bare
#: `attack_sweep/` or `md_residence/` -- because that is the half that produces a
#: wrong answer rather than an inconvenient one: unscoped, those directories are
#: every screen that has ever run, so the page renders and describes the wrong
#: campaign with nothing saying so.
#:
#: They are NOT yet held to the topic and dataset-root rules, and that is a
#: statement about work outstanding rather than about what is acceptable. Each
#: still carries stale topic literals -- `integration/app/app.py` reads
#: `nac_v2_poses` and `nac_v3_poses` while the run is nac_v5, two topics behind.
#: Promote a file into PIPELINE_PATH when those are fixed too; #74 is the list.
RUN_READERS = [
    "integration/app/app.py",          # the Streamlit GUI members are told to open
    "scripts/pose_modes_report.py",
    "scripts/shortlist_report.py",
    "shared/pipeline_schematic.py",
]

#: Directories every run would share if they were not scoped.
SHARED_DIRS = ("attack_sweep", "md_residence", "mdprio_reports", "bpmd",
               "sweep_gaps")


def _code(path: str) -> list[tuple[int, str]]:
    """Source lines with comments and docstring bodies removed.

    The prose in this project NAMES the defect in order to explain it -- every
    fix cites the directory it stopped reading. Matching that would make these
    tests unfailable, which is the mistake an earlier test in this suite had to
    correct for.
    """
    out, in_doc = [], False
    for i, l in enumerate((REPO / path).read_text(errors="replace").splitlines(), 1):
        t = l.strip()
        if t.count('"""') == 1:
            in_doc = not in_doc
            continue
        if in_doc or t.startswith("#") or t.startswith('"""') or t.startswith("#:"):
            continue
        out.append((i, l))
    return out


@pytest.mark.parametrize("path", PIPELINE_PATH + RUN_READERS)
def test_no_module_names_a_shared_run_directory(path):
    """A bare `attack_sweep/` is every screen's; `attack_sweep_<topic>/` is one.

    MATCHED ON THE JOINED SOURCE, not line by line. The first version of this
    test scanned lines, and two leaks in `mdprio_report` survived it because
    their paths were split across concatenated string literals -- one reading
    every screen's `rank_v2` tables, one reading every screen's `attack_sweep/`.
    A guard that a line break defeats is not a guard.
    """
    joined = " ".join(l for _i, l in _code(path))
    joined = re.sub(r'"\s*\n?\s*"', "", joined)      # glue "a" "b" into "ab"
    bad = []
    for d in SHARED_DIRS:
        if re.search(rf'Topic\(\s*[^,]+,\s*["\']{d}["\']\s*\)', joined) \
           or re.search(rf'/\s*["\']{d}(/|["\'])', joined) \
           or re.search(rf'["\']{d}/', joined):
            bad.append(d)
    assert not bad, f"{path} names shared run directory(ies) directly: {bad}"


@pytest.mark.parametrize("path", PIPELINE_PATH)
def test_no_module_reads_a_rank_table_without_the_topic(path):
    """`rank_v2_<tier>_<topic>_<score>_*.csv` -- the topic is what makes it this
    run's. A glob without it matched every screen's tables and `[-1]` chose by
    string order."""
    joined = re.sub(r'"\s*\n?\s*"', "", " ".join(l for _i, l in _code(path)))
    for m in re.finditer(r'rank_v2_\{?[^"\']{0,24}?\*\.csv', joined):
        frag = m.group(0)
        assert "topic" in frag or "{rp.topic" in joined[max(0, m.start()-60):m.start()], \
            f"{path}: rank table glob without a topic: {frag}"


@pytest.mark.parametrize("path", PIPELINE_PATH)
def test_no_module_hardcodes_a_topic(path):
    """`run.topic` is the whole ceremony for a fresh screen (D0080). A literal
    `nac_v4` in code is a run that cannot be started over."""
    bad = [f"{path}:{i}: {l.strip()[:80]}" for i, l in _code(path)
           if re.search(r'["\'][a-z_]*nac_v[0-9]', l)]
    assert not bad, "hardcoded topic:\n" + "\n".join(bad)


@pytest.mark.parametrize("path", PIPELINE_PATH)
def test_no_module_hardcodes_the_dataset_root(path):
    """The tool is meant to screen OTHER datasets. A literal
    `/data/lab_vm/.../inhibition` is this dataset, on this filesystem, forever.
    `shared/run_paths` is where the root is resolved."""
    if path in ("shared/run_paths.py",):
        return                     # this is where the root is allowed to live
    bad = [f"{path}:{i}: {l.strip()[:80]}" for i, l in _code(path)
           if "/data/lab_vm" in l and "inhibition" in l]
    assert not bad, "hardcoded dataset root:\n" + "\n".join(bad)


def test_run_paths_scopes_every_directory_a_run_owns():
    from shared import run_paths as rp
    t = rp.topic()
    for fn in ("sweep_dir", "residence_dir", "reports_dir", "bpmd_dir",
               "worklist_dir", "poses_dir", "allposes_dir",
               "sweep_work", "residence_work", "bpmd_work"):
        assert t in str(getattr(rp, fn)()), f"{fn}() is not scoped to {t}"


def test_the_pipeline_path_list_covers_every_declared_stage():
    """A stage added without adding its script here would be unguarded."""
    from shared import pipeline as pl
    for s in pl.STAGES:
        cmd = " ".join(s.launch())
        named = [p for p in PIPELINE_PATH if Path(p).name in cmd]
        assert named, f"stage {s.name} launches {cmd}, which no test guards"


def test_both_servers_are_threaded():
    """`HTTPServer` serves ONE request at a time. A browser holding a connection
    open, or a slow transfer of a 100 MB report, blocks every other request and
    the whole GUI stops answering -- it looks down while the process is fine.

    `python -m http.server` has used ThreadingHTTPServer since 3.7, so replacing
    it with the plain class to add no-store headers silently downgraded that."""
    for f in ("scripts/serve_reports.py", "scripts/pipeline.py"):
        src = (REPO / f).read_text()
        assert "ThreadingHTTPServer" in src, f
        code = "\n".join(l for l in src.splitlines()
                         if not l.strip().startswith("#"))
        assert re.search(r"(?<!Threading)HTTPServer\(", code) is None, \
            f"{f} instantiates a serial HTTPServer"
