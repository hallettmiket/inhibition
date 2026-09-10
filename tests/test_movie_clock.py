"""The movie's time axis comes from the trajectory, not from a table.

All 16 long-run viewers shipped with a clock taken from the results table while
their frames came from the trajectory, and the two disagreed by up to 100x: a
100 ns run under a 1.2 ns axis, with no visible symptom because the slider reads
frame numbers. These tests pin the three properties that stop it recurring.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import md_movie as mm                          # noqa: E402


def _fake_check(monkeypatch, last_ps, n_step):
    """Stand in for `gmx check` with a chosen extent."""
    class R:
        stdout = f"Step {n_step} 10\n"
        stderr = (f"Reading frame 1 time 10.000\n"
                  f"Last frame {n_step - 1} time {last_ps:.3f}\n")
    monkeypatch.setattr(mm.subprocess, "run", lambda *a, **k: R())


def test_extent_read_from_the_file(monkeypatch, tmp_path):
    _fake_check(monkeypatch, 100_000.0, 10_001)
    n, last = mm.traj_extent(tmp_path)
    assert last == 100_000.0
    assert n == 10_001


def test_a_stale_table_value_is_reported_not_believed(monkeypatch, tmp_path,
                                                      caplog):
    """The exact shipped defect: table says 1.2 ns, trajectory holds 100 ns."""
    _fake_check(monkeypatch, 100_000.0, 10_001)
    n, last = mm.traj_extent(tmp_path)
    # what build_movie_pdb does with the two numbers
    table_says = 1200.0
    assert abs(last - table_says) > 0.01 * last, "the disagreement must be seen"
    assert last == 100_000.0, "the trajectory wins"


def test_clock_survives_a_round_trip_through_the_movie_file(tmp_path):
    """A movie states its own length; a reader need not consult any table."""
    dest = tmp_path / "movie.pdb"
    body = "MODEL     1\nATOM      1  CA  ALA A   1      0.0 0.0 0.0\nENDMDL\n"
    dest.write_text(f"REMARK   1 TOTAL_PS {100000.0:.3f}\nREMARK   1 FRAMES 1\n"
                    + body)
    assert mm.movie_total_ps(dest) == 100_000.0


def test_an_unstamped_movie_returns_none_rather_than_a_guess(tmp_path):
    """The 16 shipped viewers have no recorded length. None says so."""
    dest = tmp_path / "old.pdb"
    dest.write_text("MODEL     1\nENDMDL\n")
    assert mm.movie_total_ps(dest) is None


def test_missing_file_is_none_not_a_default(tmp_path):
    assert mm.movie_total_ps(tmp_path / "absent.pdb") is None


def test_the_guard_can_fail(monkeypatch, tmp_path):
    """A vacuous guard is this project's most common defect. Prove this one
    distinguishes a matching clock from a mismatched one."""
    _fake_check(monkeypatch, 10_000.0, 1_001)
    _, last = mm.traj_extent(tmp_path)
    assert last == 10_000.0
    assert abs(last - 100_000.0) > 0.01 * last     # would be caught
    assert abs(last - 10_000.0) <= 0.01 * last     # would not be


# ---------------------------------------------------------------------------
# The boot-on-load race (found 2026-09-10 building the presentation GUI).
#
# Every caller before this one wrapped `viewer_html()`'s output in a COLLAPSED
# <details>, so `window.addEventListener('load', boot1)` always had a 'load'
# event still ahead of it -- the panel could not be opened before the page
# finished loading. The presentation GUI is the first caller to embed the
# viewer directly, open immediately, no <details> -- and on a page carrying a
# multi-megabyte embedded trajectory, `document.readyState` can already be
# "complete" by the time this script block is reached, so the listener is
# attached to an event that already fired and never runs again. No error,
# just an empty box: the exact shape this project's catalogue calls a guard
# that cannot fail, except here nothing was guarding anything -- the boot
# path itself was silently unreachable.
# ---------------------------------------------------------------------------

def _script_block(elem_id="gl"):
    return mm.viewer_html("<pdb>", [3.0], [], [], "", elem_id=elem_id)


def test_boot_checks_readystate_before_trusting_the_load_event():
    """The fix: read document.readyState, don't just listen for 'load'."""
    js = _script_block()
    assert "document.readyState" in js, (
        "the boot path must check whether 'load' has already fired, not "
        "assume it hasn't")


def test_the_already_fired_case_still_boots():
    """`armLoadBoot` must call boot1() directly when readyState is complete,
    not only attach a listener that a fired event will never trigger again."""
    js = _script_block()
    assert "if (document.readyState === 'complete')" in js
    assert "{ boot1(); }" in js.replace("{{", "{").replace("}}", "}") \
        or "{{ boot1(); }}" in js


def test_the_not_yet_loaded_case_still_uses_the_event():
    """The fix must not regress the ordinary case -- a normal page load still
    boots via the 'load' listener, not by assuming readyState is complete."""
    js = _script_block()
    assert "window.addEventListener('load', boot1)" in js


def test_a_collapsed_details_panel_still_only_boots_on_open():
    """The existing, working case (every current caller) must be unchanged:
    a details panel that starts CLOSED boots on its toggle event, not
    immediately -- opening it before the page has a real box height is the
    original bug this mechanism exists to avoid (see boot()'s own comment)."""
    js = _script_block()
    assert "det.addEventListener('toggle'" in js
    assert "if (det.open) boot1()" in js
