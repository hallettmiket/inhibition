"""
Purpose: the orientation doc's measured numbers must match the data, or the suite fails.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: docs/state_of_the_project.md, the latest frame per experiment
Output: pass/fail

ISSUE #11, and the load-bearing half of it.

`state_of_the_project.md` is this project's tier-II memory: it is what a fresh
Claude Code session reads to learn what is true before touching anything. It
says of itself that it *"drifted badly within 24 h of being written and a new
maintainer read it as fact."* A stale number here is therefore not a
documentation nit -- it is bad context handed to whoever works next.

WHY THE REFRESH SCRIPT ALONE IS NOT ENOUGH. A generator only helps if somebody
runs it, and nobody runs it exactly when the numbers are moving fastest, which
is when drift happens. So this test fails the suite when the document no longer
matches the frames. The script makes the fix a one-liner; this makes skipping
it impossible.

PROOF THAT IT CATCHES SOMETHING REAL: on the day it was written, the hand-
maintained line claimed **59,323** molecules across the six T_2 variants. The
measured total was **60,123**. Eight hundred molecules of drift, in the
document whose stated purpose is to stop people being misled.

WHAT IT DOES NOT CHECK. Only the generated AUTO blocks. The prose sections --
"what is established", "what is ruled out", "what to do next" -- are judgements
and are deliberately outside its scope; asserting on them would either be
vacuous or would freeze conclusions that are supposed to change.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DOC = REPO / "docs" / "state_of_the_project.md"
SCRIPT = REPO / "scripts" / "refresh_orientation.py"

pytestmark = pytest.mark.skipif(
    not Path("/data/lab_vm/append_only/inhibition").is_dir(),
    reason="frames are not on this machine")


def _auto_keys() -> list[str]:
    text = DOC.read_text(encoding="utf-8")
    import re
    return re.findall(r"<!-- AUTO:([a-z0-9_]+):BEGIN -->", text)


def test_the_doc_has_the_auto_blocks_it_is_supposed_to_have():
    """Guard against the guard being vacuous.

    If the markers were removed, `--check` would raise rather than pass -- but
    a future refactor that renamed them could otherwise leave this file
    checking nothing at all. Name what must be present.
    """
    keys = _auto_keys()
    assert set(keys) >= {"arms", "t2", "decisions"}, (
        f"orientation doc is missing AUTO blocks; found {keys}")


def test_the_measured_numbers_in_the_doc_are_current():
    """The whole point. Fails when a frame moved and the doc did not."""
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--check"],
        capture_output=True, text=True, cwd=REPO)
    assert proc.returncode == 0, (
        "docs/state_of_the_project.md no longer matches the data on disk.\n"
        "Run:  python3 scripts/refresh_orientation.py\n\n"
        f"{proc.stdout}{proc.stderr}")


def test_the_check_can_actually_fail(tmp_path, monkeypatch):
    """Rule 3 of `how_this_project_breaks.md`: check the guard can fail.

    Corrupt a copy of the doc and assert `--check` rejects it. Without this,
    a `--check` that silently returned 0 -- for instance because the marker
    regex stopped matching -- would pass forever and look like a healthy guard.
    """
    doc_copy = tmp_path / "state_of_the_project.md"
    text = DOC.read_text(encoding="utf-8")
    # Replace the generated arms table with something plainly wrong.
    import re
    broken = re.sub(r"(<!-- AUTO:arms:BEGIN -->).*?(<!-- AUTO:arms:END -->)",
                    r"\1\nnot the real numbers\n\2", text, flags=re.DOTALL)
    assert broken != text, "failed to corrupt the doc; the test is vacuous"
    doc_copy.write_text(broken, encoding="utf-8")

    import importlib.util
    spec = importlib.util.spec_from_file_location("refresh", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setattr(mod, "STATE_DOC", doc_copy)
    monkeypatch.setattr(sys, "argv", ["refresh_orientation.py", "--check"])

    with pytest.raises(SystemExit) as exc:
        mod.main()
    assert "STALE" in str(exc.value)
