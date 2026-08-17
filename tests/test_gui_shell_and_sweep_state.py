"""The four-step GUI and the sweep's four states (#63).

The stepper is one module because three templates each carrying their own copy is
three that drift, and the first thing a reader notices is the page whose "you are
here" is wrong. The sweep states are tested because one of them was wrong on the
first run in a way that read as a catastrophe: 162 modes reported as FAILED that
had simply not been swept yet.
"""

from __future__ import annotations

import re
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


# --------------------------------------------------------------------------
# the sweep page's pose viewer -- the four ways a 3Dmol panel comes up blank
# --------------------------------------------------------------------------

def test_the_pose_viewer_calls_render():
    """3Dmol draws nothing without it. Cost: one round of "I can't see anything"."""
    from shared import pose_viewer as pv
    assert "render()" in pv.mount_js(113, "[]")


def test_the_pose_container_has_an_explicit_height():
    """A zero-height container yields a 0x0 canvas and no error."""
    from shared import pose_viewer as pv
    assert re.search(r"\.pvbox\{[^}]*height:\d+px", pv.CSS)


def test_the_canvas_is_sized_after_layout_settles():
    """Sizing before the grid settles gives an off-centre, unrotatable view."""
    from shared import pose_viewer as pv
    js = pv.mount_js(113, "[]")
    assert js.count("requestAnimationFrame") >= 2, "needs the double rAF"


def test_the_page_puts_the_library_and_data_before_the_code_that_uses_them():
    """Both inversions produce a silently blank viewer."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "build_gui", REPO / "scripts" / "build_gui.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    page = m._page("t", "sweep.html", {}, "<div id='b'></div>",
                   extra_js="USECODE", head_js="<script>LIBRARY</script>",
                   tail_data="<pre id='recpdb'>DATA</pre>")
    assert page.index("LIBRARY") < page.index("USECODE")
    assert page.index("DATA") < page.index("USECODE")


def test_a_missing_pose_asset_is_reported_not_left_blank():
    from shared import pose_viewer as pv
    assert "pvempty" in pv.CSS


def test_poses_are_separate_models_not_frames():
    """`addModelsAsFrames` builds ONE model with n frames, so `{model: i+1}`
    addresses nothing and 3Dmol dereferences undefined -- "Cannot read
    properties of undefined (reading 'setStyle')". Frames animate; these are
    alternatives to be styled independently."""
    from shared import pose_viewer as pv
    # Comments are stripped: the prose above the fix names the broken call in
    # order to explain it, and matching that would make this test unfailable.
    js = "\n".join(ln for ln in pv.mount_js(113, "[]").splitlines()
                   if not ln.lstrip().startswith("//"))
    assert "addModelsAsFrames" not in js
    assert "addModel(" in js
    assert "split('ENDMDL')" in js, "poses are not split per MODEL block"


def test_the_mode_is_read_from_the_model_record_not_counted():
    """Counting file position is #53, one layer down."""
    from shared import pose_viewer as pv
    assert "MODEL" in pv.mount_js(113, "[]")
    assert "exec(b)" in pv.mount_js(113, "[]")


def test_sweep_assets_are_named_per_mode_not_per_molecule():
    """Two modes of one molecule are DIFFERENT trajectories. Naming an asset
    after the parent would let one mode's movie stand for another's, which is
    the collision `shared/mode_key` exists to prevent."""
    src = (REPO / "scripts" / "sweep_assets.py").read_text()
    assert 'OUT / f"{ident}.pdb"' in src and 'OUT / f"{ident}.png"' in src
    assert 'OUT / f"{parent}' not in src


def test_the_distance_panel_reads_the_window_from_the_criterion():
    """A hand-typed 2.8-4.2 here could drift from the number that scores."""
    src = (REPO / "scripts" / "sweep_assets.py").read_text()
    assert "nac.NAC_DIST_MIN" in src and "nac.NAC_DIST_MAX" in src


def test_only_a_finished_sweep_supplies_assets():
    """A partial trajectory rendered as a whole one is the #53 neighbourhood."""
    src = (REPO / "scripts" / "sweep_assets.py").read_text()
    assert "Finished mdrun" in src


def test_every_page_that_renders_the_nav_also_ships_its_css():
    """THE REGRESSION. `mdprio_combine` interpolated the stepper markup but its
    CSS line was deleted when the shell stylesheet was moved into
    `shared/results_shell` -- the extraction replaced everything between <style>
    and </style>, and `{_stepcss}` was inside that range. The page still built,
    still had the nav, and rendered it as a run-on line of plain hyperlinks.
    @tt8804: "the nav bar on the gui is just hyperlinks".

    Markup without styles is the failure mode, so the test is markup implies
    styles -- checked at the BUILDER, not on built artefacts, so it runs
    anywhere."""
    builders = {
        "scripts/mdprio_combine.py": ("_stepnav", "_stepcss"),
        "shared/mode_ranking.py": ("__STEPNAV__", "__STEPCSS__"),
        "scripts/build_gui.py": ("gs.nav(", "gs.CSS"),
        "scripts/sweep_combine.py": ("gs.nav(", "gs.CSS"),
    }
    for path, (nav, css) in builders.items():
        src = (REPO / path).read_text()
        if nav not in src:
            continue                      # this builder does not render the nav
        assert css in src, (
            f"{path} renders the step nav ({nav}) but never interpolates its "
            f"CSS ({css}) -- it will render as plain hyperlinks")
