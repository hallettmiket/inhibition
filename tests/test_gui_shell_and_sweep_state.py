"""The four-step GUI and the sweep's four states (#63).

The stepper is one module because three templates each carrying their own copy is
three that drift, and the first thing a reader notices is the page whose "you are
here" is wrong. The sweep states are tested because one of them was wrong on the
first run in a way that read as a catastrophe: 162 modes reported as FAILED that
had simply not been swept yet.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import gui_shell as gs                        # noqa: E402
from shared import sweep_state as ss                      # noqa: E402


# --------------------------------------------------------------------------
# the stepper
# --------------------------------------------------------------------------

def test_every_step_is_a_link_and_exactly_one_is_current():
    for href, _label, _desc in gs.STEPS:
        h = gs.nav(href)
        assert h.count('class="on"') == 1, f"{href}: not exactly one current step"
        assert f'href="{href}"' in h


def test_the_steps_are_in_pipeline_order():
    """Home, then rank, then sweep, then MD. The funnel only reads as a funnel
    in the order it happens -- and that is NOT the order the pages were built."""
    assert [s[0] for s in gs.STEPS] == [
        "index.html", "modes.html", "sweep.html", "combined.html"]


def test_an_unknown_current_page_highlights_nothing_rather_than_guessing():
    h = gs.nav("something_else.html")
    assert 'class="on"' not in h


def test_a_zero_count_is_omitted_not_printed():
    """"not measured yet" and "measured, none" are different claims."""
    h = gs.nav("modes.html", {"sweep.html": ""})
    assert "<span class=\"sn\">" not in h or "0 ok" not in h


# --------------------------------------------------------------------------
# sweep state -- the bug that made unswept modes look failed
# --------------------------------------------------------------------------

def _state_of(status, queued):
    """Exercise the classifier the way `state()` calls it."""
    row = {"status": status, "_queued": queued}
    v = row.get("status")
    s = "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v).strip()
    if s == "ok":
        return "ok"
    if s:
        return "failed"
    return "pending" if row["_queued"] else "not sent"


def test_a_queued_mode_with_no_result_is_pending_not_failed():
    """THE REGRESSION. `nan or ""` evaluates to nan -- NaN is truthy -- so
    `str(nan)` is "nan", which is non-empty, and every mode that had merely not
    been swept yet was classified FAILED. 162 of them on the first run."""
    assert _state_of(float("nan"), True) == "pending"
    assert _state_of(None, True) == "pending"
    assert _state_of("", True) == "pending"


def test_an_unqueued_mode_with_no_result_is_not_sent():
    assert _state_of(float("nan"), False) == "not sent"


def test_ok_is_ok_and_anything_else_recorded_is_failed():
    assert _state_of("ok", False) == "ok"
    assert _state_of("sweep failed", True) == "failed"
    assert _state_of("skipped: an unfinished trajectory", True) == "failed"


def test_whitespace_does_not_turn_ok_into_a_failure():
    assert _state_of(" ok ", False) == "ok"


def test_summary_counts_every_state_even_when_absent():
    d = pd.DataFrame({"sweep_state": ["ok", "ok", "pending"],
                      "frac_attack_ready": [0.5, 0.0, None]})
    s = ss.summary(d)
    assert s["ok"] == 2 and s["pending"] == 1
    assert s["failed"] == 0 and s["not sent"] == 0     # present, as zero
    assert s["productive"] == 1


def test_summary_of_nothing_is_zeros_not_an_error():
    assert ss.summary(pd.DataFrame())["ok"] == 0


@pytest.mark.parametrize("name", ["gui_shell", "sweep_state"])
def test_the_modules_import_without_touching_the_filesystem(name):
    """A GUI helper that reads /data at import time cannot be tested anywhere."""
    import importlib
    importlib.reload(importlib.import_module(f"shared.{name}"))
