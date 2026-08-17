"""The pipeline: stages that know their own state, and the rules they enforce.

WHY THIS EXISTS. Running 3.1.0 by hand needed about a dozen bespoke shell
launchers, and the stage scripts were never what broke -- the glue was,
differently each time. Every failure had the same shape: a value taken by
POSITION, by RECENCY, by DEFAULT, or by a NAME THAT DOES NOT EXIST, failing
silently because an empty result and a broken query are indistinguishable when
nothing asserts the difference.

Each test below pins one of those, so the pipeline cannot reacquire it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import pipeline as pl                     # noqa: E402


def test_every_stage_declares_what_it_needs_and_the_graph_is_acyclic():
    seen = set()
    for s in pl.STAGES:
        for n in s.needs:
            assert n in pl.BY_NAME, f"{s.name} needs unknown stage {n}"
            assert n in seen, f"{s.name} needs {n}, which is declared after it"
        seen.add(s.name)


def test_the_stage_order_is_the_cascade():
    assert [s.name for s in pl.STAGES] == [
        "screen", "rank", "worklist", "sweep", "production", "bpmd"]


def test_unknown_is_a_state_distinct_from_zero():
    """THE ONE THAT COST THE MOST. A survivor query looked for an RMSD column
    `attack_sweep` does not write, found nothing, and reported "no survivors"
    for eight hours while two sat under the bar. A probe that cannot be
    evaluated must not present as a probe that found nothing."""
    class Boom:
        name = "boom"; title = "Boom"; needs = (); gpus = 0; proc = ""; note = ""
        @staticmethod
        def probe():
            raise pl.StageError("cannot evaluate")
    st = pl.state(Boom())
    assert st["state"] == "unknown"
    assert st["done"] is None and st["total"] is None
    assert "cannot evaluate" in st["error"]


def test_survivors_raises_rather_than_returning_empty_when_it_cannot_read():
    """Checked at the source, since the failure only reproduces with a real
    sweep tree: the function must have a raise on the unresolved path, not a
    return."""
    src = (REPO / "shared" / "pipeline.py").read_text()
    body = src[src.index("def survivors("):src.index("def p_production(")]
    assert "raise StageError" in body
    assert "unresolved" in body


def test_the_worklist_is_chosen_by_the_run_not_by_recency():
    """`ls -t | head -1` returned the PREVIOUS screen's worklist when this
    run's ranking had failed, and 170 of its modes were reported as this run's
    selection. A worklist only belongs to this run if it postdates this run's
    ranking table."""
    src = (REPO / "shared" / "pipeline.py").read_text()
    body = src[src.index("def worklist_path("):src.index("@dataclass")]
    assert "getmtime" in body and "newest_rank" in body
    assert rp_topic_in(body), "the glob is not scoped to this run's topic"


def rp_topic_in(body: str) -> bool:
    return "rp.topic()" in body


def test_children_are_killed_parents_first():
    """`md_residence_3ikd` spawns `gmx` and relaunches it for the next stage of
    its own chain, so killing `gmx` first makes a new one appear ~13 s later --
    a "stopped" fleet held eight GPUs at 90% that way."""
    order = list(pl._KILL_ORDER)
    assert order.index("md_residence_3ikd.py") < order.index("gmx mdrun")
    assert order.index("attack_sweep.py") < order.index("md_residence_3ikd.py")


def test_stopping_never_touches_another_users_gmx():
    """The box is shared and other people's jobs sit on the same cards."""
    src = (REPO / "shared" / "pipeline.py").read_text()
    body = src[src.index("def stop("):]
    assert "/proc/" in body and "cwd" in body
    assert "rp.topic() not in cwd" in body


def test_process_matching_excludes_the_querying_shell():
    """`pgrep -f` matches the shell running the query as often as the target,
    which made every count off by one and made a pkill kill its own caller."""
    src = (REPO / "shared" / "pipeline.py").read_text()
    body = src[src.index("def _procs("):src.index("def running(")]
    assert "os.getpid()" in body and "bash -c" in body


def test_a_stage_refuses_to_start_before_its_inputs_are_done():
    src = (REPO / "shared" / "pipeline.py").read_text()
    body = src[src.index("def start("):src.index("_KILL_ORDER")]
    assert "already running" in body
    assert "needs" in body and "raise StageError" in body


def test_the_sweep_stage_reads_the_worklist_by_column_name():
    """Positional parsing took field 5 as pose_rank when it holds global_rank,
    so all 24 launched jobs asked for ranks in the hundreds."""
    src = (REPO / "scripts" / "pipeline_stage.py").read_text()
    body = src[src.index("def stage_sweep("):src.index("def _swept(")]
    assert "r.pose_rank" in body and "r.parent_ident" in body
    assert "split(\",\")" not in body


def test_resume_is_keyed_on_the_pair_not_the_molecule():
    """One swept mode must not stand for every mode of its molecule."""
    src = (REPO / "scripts" / "pipeline_stage.py").read_text()
    body = src[src.index("def _swept("):src.index("def stage_production(")]
    assert "parent_ident" in body and "pose_rank" in body


def test_status_serialises_for_the_dashboard():
    """The page renders from this JSON; a value it cannot encode is a blank
    dashboard."""
    import json
    st = pl.status()
    json.loads(json.dumps(st, default=str))
    assert {"topic", "stages", "spec"} <= set(st)
    for s in st["stages"]:
        assert {"name", "title", "state", "done", "total", "ready"} <= set(s)


@pytest.mark.parametrize("stage", [s.name for s in pl.STAGES])
def test_every_stage_has_a_launch_command(stage):
    s = pl.BY_NAME[stage]
    assert s.launch is not None, f"{stage} cannot be started"
    cmd = s.launch()
    assert cmd and all(isinstance(x, str) for x in cmd)
