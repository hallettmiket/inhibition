"""
Purpose: a frame produced by code that has since changed must have a MEASURED impact recorded.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: scripts/check_frame_code_currency.py
Output: pass/fail

WHY THIS EXISTS. A fix lands in a generation stage, the stage is never re-run,
and the fix is inert on the data everyone is looking at. Nothing announces it --
the frame is present, populated and plausible, and the code reads as though the
fix is in force. Found twice by hand on 2026-08-04 (T_3's pocket ceiling, T_2
ATRA's), and the check then surfaced two more nobody had looked at.

The test does NOT demand that stale stages be re-run. Re-running a 16,806-
molecule CReM expansion to pick up a comment change would be absurd, and three
of the four flags dissolved once measured. It demands that staleness be
MEASURED and the number written down -- "we checked and it does not matter" is
worth nothing without the impact attached.

Allowlist, not denylist (D0051): a stage nobody has looked at fails, rather
than passing because it was not anticipated.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "check_frame_code_currency.py"

pytestmark = pytest.mark.skipif(
    not Path("/data/lab_vm/append_only/inhibition").is_dir(),
    reason="frames are not on this machine")


def test_no_generation_frame_is_stale_without_a_measured_impact():
    proc = subprocess.run([sys.executable, str(SCRIPT), "--check"],
                          capture_output=True, text=True, cwd=REPO)
    assert proc.returncode == 0, (
        "a generation stage's code changed after its frame was produced, and "
        "nobody has recorded what that changed in the data.\n"
        "Re-run the stage, or measure the impact and add it to ACKNOWLEDGED "
        "in scripts/check_frame_code_currency.py.\n\n"
        f"{proc.stdout}{proc.stderr}")


def test_every_acknowledgement_carries_a_measurement():
    """Check the guard can't be silenced with a shrug.

    An ACKNOWLEDGED entry that just says "fine" would pass the test above while
    recording nothing. Require the word MEASURED and a digit -- the number is
    the whole point.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("cfc", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    assert mod.ACKNOWLEDGED, "no acknowledgements to check; test would be vacuous"
    for stage, note in mod.ACKNOWLEDGED.items():
        assert "MEASURED" in note, f"{stage}: acknowledgement records no measurement"
        assert any(c.isdigit() for c in note), f"{stage}: no number in the impact"
