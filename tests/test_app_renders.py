"""
Purpose: Every panel must actually RENDER, with and without a curation filter.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: integration/app/app.py, driven through streamlit's own AppTest harness
Output: pass/fail

WHY THIS EXISTS. `tests/test_app_names.py` catches a name a panel never binds,
which is one failure mode. It cannot catch a panel that binds every name and
then raises on the first row of real data, because a Streamlit panel's body only
executes when someone clicks its tab. Both of this project's GUI bugs so far
reached a user through that gap: the pose lookup that returned nothing, and the
curation filter that reached exactly one panel.

`AppTest` runs the script headlessly — no browser, no port — and surfaces any
exception the script raised. That is the difference between "it compiles" and
"it renders".

WHICH ENVIRONMENT THIS RUNS IN. streamlit lives in `dwi_gui`, not in the
`dwi_cheminf` env the suite normally runs under, so these skip by default and
execute when the suite is pointed at the GUI env. They are worth having in the
file either way: skipped-because-absent is a different signal from not-written.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tests"))
import gui_harness  # noqa: E402
APP = REPO / "integration" / "app" / "app.py"
sys.path.insert(0, str(REPO / "integration" / "app"))

#: Constraints that bite on the real shortlists — a spec that removed nothing
#: would let a filter that does not carry over pass every assertion below.
SPEC = "no chlorine\nmw < 450"
BAD_SPEC = "no chlorien"


def _app_test(panel: str, spec: str | None = None):
    pytest.importorskip("streamlit", reason="streamlit lives in the dwi_gui env")
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(str(APP), default_timeout=300)
    at.run()
    assert not at.exception, f"app failed on first run: {at.exception}"
    at.sidebar.radio[0].set_value(panel).run()
    assert not at.exception, f"{panel} raised: {at.exception}"
    if spec is not None:
        try:
            gui_harness.set_spec(at, spec)
        except gui_harness.HarnessLimitation as exc:
            pytest.skip(str(exc))
        assert not at.exception, f"{panel} raised under curation: {at.exception}"
    return at


def _messages(at) -> str:
    out: list[str] = []
    for kind in ("info", "success", "warning", "error", "markdown", "caption"):
        out += [str(e.value) for e in getattr(at, kind)]
    return "\n".join(out)


def _panels(filtered: bool | None = None) -> list[str]:
    """Panel names from the single declaration in `curate.PANEL_SCOPE`.

    Read at collection time so the parametrised cases below track the scope
    table rather than a second hand-maintained list that could disagree with it.
    """
    import curate

    return [s.panel for s in curate.PANEL_SCOPE
            if filtered is None or s.filtered is filtered]


@pytest.mark.parametrize("panel", _panels())
def test_panel_renders_without_a_filter(panel):
    _app_test(panel)


@pytest.mark.parametrize("panel", _panels())
def test_panel_renders_under_a_curation_filter(panel):
    _app_test(panel, SPEC)


@pytest.mark.parametrize("panel", _panels())
def test_panel_renders_when_the_constraint_is_unparseable(panel):
    """A typo must produce a refusal, not a traceback and not silent data."""
    at = _app_test(panel, BAD_SPEC)
    assert "not understood" in _messages(at), (
        f"{panel} neither filtered nor refused — a mis-parsed constraint that "
        "quietly filters nothing is indistinguishable from a working one")


@pytest.mark.parametrize("panel", _panels(filtered=True))
def test_a_filtered_panel_shows_the_persistent_indicator(panel):
    """THE ISSUE #3.2 REGRESSION TEST, at the rendered-output level."""
    at = _app_test(panel, SPEC)
    assert "Curation active" in _messages(at), (
        f"{panel} is declared filtered but rendered no curation banner — the "
        "reader cannot tell how many rows were removed")


@pytest.mark.parametrize("panel", _panels(filtered=False))
def test_an_unfiltered_panel_says_so_out_loud(panel):
    """Silence is exactly how the original bug read to the user."""
    at = _app_test(panel, SPEC)
    assert "not applied here" in _messages(at)


def test_the_dossier_stops_offering_excluded_candidates():
    """The panel the issue was reported against."""
    plain = _app_test("Candidate dossier")
    filtered = _app_test("Candidate dossier", SPEC)
    n_all = len(plain.selectbox[1].options)
    n_curated = len(filtered.selectbox[1].options)
    assert 0 < n_curated < n_all, (
        f"the dossier offered {n_curated} of {n_all} candidates under a filter "
        "that removes rows from the shortlist — it is not carrying over")


def test_the_dossier_escape_hatch_flags_what_it_brings_back():
    """Inspecting a rejected molecule must never look like endorsing it."""
    at = _app_test("Candidate dossier", SPEC)
    n_curated = len(at.selectbox[1].options)
    at.checkbox[0].set_value(True).run()
    assert not at.exception
    options = [str(o) for o in at.selectbox[1].options]
    flagged = [o for o in options if o.startswith("✗")]
    assert len(options) > n_curated and flagged, (
        "excluded candidates came back unlabelled — a reader cannot tell them "
        "from candidates still in contention")


def test_the_pose_viewer_offers_every_mode_and_an_export():
    """Issue #3.1(b) and (c), at the rendered-output level."""
    at = _app_test("Candidate dossier")
    modes = [r for r in at.radio if r.label == "which modes"]
    assert modes, "no pose-mode selector rendered; only one pose is reachable"
    assert "overlay all" in list(modes[0].options)
    assert list(at.get("download_button")), "no external hand-off offered"
