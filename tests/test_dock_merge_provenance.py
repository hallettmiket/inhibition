"""
Purpose: the docking merge must name the frame it consumed, on every call path.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: shared.noncovalent_dock_run
Output: pass/fail

WHY THIS EXISTS. `merge_poses_onto_frame` was split out of `run()` in 2a22970
so that a chunked run and a single-GPU run share one merge. The reference to
`frame_path` moved with it; the BINDING did not. `frame_path` is assigned only
inside the `if df is None:` branch, and the manifest call at the bottom reads it
unconditionally -- so:

* the CHUNKED path passes no `df`, binds `frame_path`, and works;
* the `run()` path passes `df`, never binds it, and raised
  `UnboundLocalError` at the manifest call.

That is the whole of T_1's and T_2's ordinary `03_dock.py` route, and it failed
AFTER the GPU run was already spent -- ~1.3 h for atra, ~10 h for du_xu. It went
unseen because the only pool docked after the refactor was liu_2024_c3, which is
large enough to require chunking, i.e. the one path that works.

WHY THE FIX IS A PARAMETER AND NOT A `dio.latest` CALL. Re-resolving the newest
frame here would make the crash go away and record whichever frame happened to
be newest at merge time rather than the one the caller actually read. That is a
provenance lie in a manifest whose entire purpose is to record the SHA-256 of
what a run consumed, and it would be silent. The caller knows; the caller passes
it. A missing `frame_path` alongside `df` now raises.

WHAT MAKES THIS TEST NON-VACUOUS. `test_stale_guard`'s first version counted
only dotted attribute reads, and a defensive `getattr` left it zero accesses to
check -- it passed for free. So the runtime test here asserts that
`write_full_frame` was ACTUALLY REACHED and that the recorded input is the frame
that was passed in. If the function returned early, `captured` stays empty and
the assertion fails rather than passing by omission.

The AST test is the class-level guard: it fails for any FUTURE caller that
supplies `df` without `frame_path`, which is the mistake itself rather than this
one instance of it.
"""

from __future__ import annotations

import ast
import pathlib
import tempfile

import pandas as pd
import pytest

from shared import io as dio
from shared import noncovalent_dock_run as ncd

REPO = pathlib.Path(__file__).resolve().parent.parent
SCAN_DIRS = ("shared", "scripts", "approaches", "integration", "tests")

FRAME = pathlib.Path(
    "/data/lab_vm/append_only/inhibition/02_t2_atra_crem/D2_32.parquet")


@pytest.fixture
def captured(monkeypatch):
    """Intercept the frame write so the test never touches the data root."""
    seen: dict = {}

    def fake_write(merged, **kwargs):
        seen.update(kwargs)
        seen["n_rows"] = len(merged)
        return pathlib.Path("/dev/null")

    monkeypatch.setattr(dio, "write_full_frame", fake_write)
    return seen


def _frame() -> pd.DataFrame:
    return pd.DataFrame({"candidate_id": ["t2_aaa", "t2_bbb"],
                         "rejected_at": [None, None]})


def test_run_path_records_the_frame_it_was_given(captured):
    """`df` + `frame_path` -> the manifest names that exact frame.

    This is the call shape `run()` uses and the one that was broken.
    """
    df = _frame()
    ncd.merge_poses_onto_frame(
        experiment="02_t2_atra_crem", approach="t2", frame_prefix="D2",
        out_dir=pathlib.Path(tempfile.mkdtemp()), elapsed=1.0, gpu=0,
        df=df, frame_path=FRAME, survivors=df)

    # Reached the manifest at all -- see "what makes this test non-vacuous".
    assert captured, ("write_full_frame was never reached; the merge returned "
                      "early and this test would otherwise pass for free")
    assert captured["inputs"]["frame"] == FRAME, (
        "the manifest must record the frame the caller read, not one resolved "
        "at merge time")
    assert captured["n_rows"] == len(df)


def test_df_without_frame_path_is_refused(captured):
    """A merge that cannot name its input frame must not write one."""
    df = _frame()
    with pytest.raises(ValueError, match="frame_path"):
        ncd.merge_poses_onto_frame(
            experiment="02_t2_atra_crem", approach="t2", frame_prefix="D2",
            out_dir=pathlib.Path(tempfile.mkdtemp()), elapsed=1.0, gpu=0,
            df=df, survivors=df)
    assert not captured, "refused merge must not have written a frame"


def _calls_to_merge(path: pathlib.Path) -> list[ast.Call]:
    try:
        tree = ast.parse(path.read_text(errors="replace"))
    except SyntaxError:
        return []
    return [n for n in ast.walk(tree)
            if isinstance(n, ast.Call)
            and isinstance(n.func, ast.Attribute)
            and n.func.attr == "merge_poses_onto_frame"]


def test_every_caller_supplying_df_also_supplies_frame_path():
    """The class-level guard: catches the NEXT caller, not just this one."""
    offenders = []
    for d in SCAN_DIRS:
        for py in (REPO / d).rglob("*.py"):
            # This file deliberately contains the defect, in
            # `test_df_without_frame_path_is_refused`, to prove the runtime
            # guard rejects it. Excluded BY NAME rather than by dropping
            # `tests/` from the scan, so a bad caller in any other test is
            # still caught.
            if py.name == pathlib.Path(__file__).name:
                continue
            for call in _calls_to_merge(py):
                kw = {k.arg for k in call.keywords if k.arg}
                if "df" in kw and "frame_path" not in kw:
                    offenders.append(f"{py.relative_to(REPO)}:{call.lineno}")

    assert not offenders, (
        "these calls pass `df` without `frame_path`, so the manifest cannot "
        f"name the frame the run consumed: {offenders}")


def test_the_ast_guard_can_actually_fail(tmp_path):
    """Check the guard can fail -- rule 3 of `how_this_project_breaks.md`.

    A guard that inspects something absent passes for free. Feed it the exact
    defect and assert it is seen.
    """
    bad = tmp_path / "bad_caller.py"
    bad.write_text(
        "ncd.merge_poses_onto_frame(experiment='e', approach='t2',\n"
        "                           frame_prefix='D2', out_dir=None,\n"
        "                           elapsed=1.0, gpu=0, df=df)\n")
    calls = _calls_to_merge(bad)
    assert len(calls) == 1, "the AST matcher did not see the call at all"
    kw = {k.arg for k in calls[0].keywords if k.arg}
    assert "df" in kw and "frame_path" not in kw, (
        "the matcher did not register the defect it exists to catch")
