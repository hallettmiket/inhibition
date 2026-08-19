"""The run spec lives in config/target.yaml, and every stage reads it there.

WHY THIS FILE EXISTS. Before a full re-run the pipeline was audited for
cohesion and the spec turned out to be scattered:

  * `--topic` carried FIVE different literals -- nac_screen_v2 and
    score_selection said "nac_v3", rank_v2 said "nac_v2", sweep_gap_worklist
    said "nac_v4" -- while config said something else again. A run launched
    without an explicit flag would have written the SCREEN to one topic and the
    WORKLIST to another, silently. D0080 exists because that already happened.
  * `md.production_ps: 100000` was read by `pipeline_schematic` -- the DIAGRAM
    -- and by nothing that runs. The runner's own default was 100.0 ps, a
    thousandth of it. The page stated 100 ns; the code would have produced
    0.1 ns; the real number came from a flag in a scratch shell script.
  * `SWEEP_PS` stayed a 10 ns literal after D0085 moved the decision to 8 ns:
    25% more GPU time per mode than the experiment concluded was needed.

Each of these is the same failure -- a value displayed in one place and applied
from another -- and each is invisible in the output, because a run at the wrong
length or the wrong topic still produces a complete-looking result.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import target_config as tc                     # noqa: E402

#: Every stage that takes --topic. A new one added without config wiring is
#: exactly the drift this file is here to stop.
TOPIC_STAGES = [
    "nac_screen_v2", "rank_v2", "score_selection",
    "backfill_anchor_p90", "sweep_gap_worklist", "rebuild_representatives",
]


@pytest.mark.parametrize("stage", TOPIC_STAGES)
def test_no_stage_carries_its_own_topic_literal(stage):
    src = (REPO / "scripts" / f"{stage}.py").read_text()
    line = next((l for l in src.splitlines() if '"--topic"' in l), None)
    assert line, f"{stage}: no --topic argument found"
    assert "default=tc.topic()" in line, (
        f"{stage} defaults --topic from `{line.strip()}` rather than config; "
        "a default in six places is six defaults")


@pytest.mark.parametrize("stage", TOPIC_STAGES + ["attack_sweep",
                                                 "md_residence_3ikd"])
def test_the_stage_can_actually_build_its_parser(stage):
    """TEXT IS NOT BEHAVIOUR. The string check above passed while
    `nac_screen_v2` and `sweep_gap_worklist` died on launch with
    `NameError: name 'tc' is not defined` -- both already imported
    target_config INSIDE a function, so my "is it imported?" guard saw the name
    and skipped adding it at module scope. The screen exited in under a second
    and the only evidence was in a log.

    Running `--help` builds every argparse default for real. A stage that
    cannot construct its own parser cannot run, whatever its source says."""
    import subprocess
    r = subprocess.run([sys.executable, str(REPO / "scripts" / f"{stage}.py"),
                        "--help"], capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, f"{stage} --help failed:\n{r.stderr[-800:]}"


def test_the_screen_and_the_worklist_cannot_disagree():
    """The specific pairing that would split a run in half: the screen writes
    poses under one topic, the worklist reads representatives from another."""
    a = (REPO / "scripts" / "nac_screen_v2.py").read_text()
    b = (REPO / "scripts" / "sweep_gap_worklist.py").read_text()
    assert "tc.topic()" in a and "tc.topic()" in b


def test_the_runner_length_is_the_configured_length():
    """100.0 vs 100_000 -- the sort of default that produces a finished-looking
    0.1 ns trajectory nobody looks at twice."""
    src = (REPO / "scripts" / "md_residence_3ikd.py").read_text()
    line = next((l for l in src.splitlines() if '"--production-ps"' in l), None)
    assert line and "default=tc.md_production_ps()" in line, \
        f"the 100 ns runner does not read md.production_ps: {line}"


def test_the_sweep_length_is_the_configured_length():
    src = (REPO / "scripts" / "attack_sweep.py").read_text()
    assert re.search(r"^SWEEP_PS\s*=\s*tc\.md_sweep_ps\(\)", src, re.M), \
        "attack_sweep carries its own sweep length"


def test_one_bar_gates_sweep_production_and_promotion():
    """The same question at three timescales, so the same number -- a second
    constant drifts the moment either is retuned."""
    for f in ("scripts/promote_to_bpmd.py", "scripts/sweep_combine.py",
              "scripts/mdprio_combine.py", "scripts/mdprio_report.py"):
        assert "sweep_survivor_rmsd_nm" in (REPO / f).read_text(), f


def test_the_spec_resolves_to_the_decisions_on_record():
    """D0087 (5 ns), D0085 (0.35 nm) and D0081 (T_4 only), read back through the
    API the runners now use. A config edit that contradicts a decision fails
    here.

    The sweep length was 8 ns under D0085 and is 5 ns under D0087, which amends
    it on length alone: the optimum was measured as a plateau (bootstrap 95% CI
    4.3-9.5 ns) and truncation is one-sided, so a shorter window can only admit
    extras, never drop a survivor. The BAR is untouched.
    """
    assert tc.md_sweep_ps() == 5_000.0
    assert tc.md_survivor_rmsd_nm() == 0.35
    assert tc.md_production_ps() == 100_000.0
    assert list(tc.get("run.tiers")) == ["T4"]
    assert set(tc.sweep_families()) == {"acrylamide", "bdhi_c4", "bdhi_c5"}


# --------------------------------------------------------------------------
# a run's topic reaches every stage, not just the first two
# --------------------------------------------------------------------------

def test_every_run_directory_moves_with_the_topic():
    """D0080 was applied to docking and ranking and stopped there, so the sweep,
    the 100 ns stage and the reports all lived in flat directories shared by
    every screen ever run. Bumping run.topic emptied one GUI page out of three
    and the other two went on showing 554 sweep rows and 647 residence rows from
    superseded runs. @tt8804: "the gui is not updated fresh"."""
    from shared import run_paths as rp
    t = rp.topic()
    for fn in ("sweep_dir", "residence_dir", "reports_dir", "bpmd_dir",
               "poses_dir", "sweep_work", "residence_work"):
        p = str(getattr(rp, fn)())
        assert t in p, f"{fn}() -> {p} does not carry the topic {t}"


def test_the_sweep_and_production_workdirs_are_different_roots():
    """The workdir is <root>/<ident>/md/rep1 regardless of tag, so sharing a
    root means a production run finds the triage run's finished prod.xtc and
    skips itself."""
    from shared import run_paths as rp
    assert rp.sweep_work() != rp.residence_work()


def test_no_gui_builder_reads_a_flat_run_directory():
    """The readers are what the user actually sees. A builder still globbing
    `attack_sweep/` or `md_residence/` shows the previous screen's results under
    the current screen's title, which is worse than showing nothing."""
    import re
    bad = re.compile(r'["\'](?:attack_sweep|md_residence|mdprio_reports|bpmd)(?:/|["\'])')
    for f in ("scripts/mdprio_combine.py", "scripts/sweep_report.py",
              "scripts/sweep_combine.py", "scripts/build_gui.py",
              "scripts/sweep_assets.py", "scripts/promote_to_bpmd.py",
              "shared/sweep_state.py"):
        src = (REPO / f).read_text()
        hits = [l.strip() for l in src.splitlines()
                if bad.search(l) and not l.lstrip().startswith("#")]
        assert not hits, f"{f} still reads a flat run directory: {hits[:2]}"
