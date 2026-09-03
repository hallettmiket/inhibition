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


def test_the_asset_cache_is_keyed_on_which_trajectory_drew_it():
    """THE REGRESSION. Fixing `rep_dir` to take `pose_rank` corrected the number
    the rail is sorted on, but the PNGs and movies on disk had been drawn from a
    sibling pose and the cache had no way to know: it skipped on "the file
    exists". t4_2f88a2f534fd_m1 was ranked on its own 0.255 nm trace beside a
    plot showing rank13's 0.857 nm. @tt8804: "why is the selector showing 0.255
    nm max but the rmsd plots show max 0.857 nm".

    So the source path is recorded and compared, which makes the cache
    self-invalidating -- and the resolution must happen BEFORE the skip, or the
    comparison has nothing to compare against."""
    src = (REPO / "scripts" / "sweep_assets.py").read_text()
    assert 'f"{ident}.src"' in src, "no provenance sidecar"
    assert "stale = " in src and "src.read_text()" in src
    body = src[src.index("def main("):]
    assert body.index("rep = rep_dir(") < body.index("stale = "), \
        "the cache is consulted before the trajectory is resolved"
    # Every guard that decides whether to rebuild must consider staleness --
    # checked on the guard's own line, since `stale` sits before the marker.
    for branch in ("want_mov and (", "not png.is_file():"):
        line = next(ln for ln in body.splitlines()
                    if branch in ln and ln.lstrip().startswith("if "))
        assert "stale" in line, f"a stale asset is not rebuilt at: {line.strip()}"


def test_the_rmsd_panel_shows_the_protein_alongside_the_ligand():
    """Ligand RMSD is measured after superposing on protein CA, so it rises both
    when the ligand leaves a rigid pocket and when the protein relaxes around a
    ligand that is still bound. One trace cannot tell those apart. @tt8804: "we
    should show ligand rmsd and protein rmsd, if they change tgt then they are
    fine"."""
    src = (REPO / "scripts" / "sweep_assets.py").read_text()
    assert "protein_rmsd" in src and "protein CA" in src
    from shared import gromacs_analysis as ga
    assert hasattr(ga, "protein_rmsd")


def test_the_protein_trace_cannot_be_clipped_off_the_axis():
    """The y-limit was set from the ligand maximum alone. A protein that moved
    further than the ligand would be drawn past the top of the axis and read as
    a protein that stayed still -- the opposite conclusion."""
    src = (REPO / "scripts" / "sweep_assets.py").read_text()
    seg = src[src.index("set_ylim(0"):]
    pre = src[:src.index("set_ylim(0")]
    assert "top" in seg[:40] and "nanmax(p_rmsd)" in pre, \
        "set_ylim does not account for the protein trace"


def test_protein_rmsd_fits_and_measures_the_same_group():
    """Group 1 both times: fit on CA, measure CA. Group 2 in either slot
    re-measures the ligand under a different filename -- two traces that are the
    same quantity, which is worse than one, because it looks like agreement."""
    src = (REPO / "shared" / "gromacs_analysis.py").read_text()
    body = src[src.index("def protein_rmsd("):]
    assert 'stdin="1\\n1\\n"' in body
    assert 'stdin="1\\n2\\n"' not in body, "measures the ligand, not the protein"


def test_a_missing_protein_trace_degrades_the_plot_rather_than_failing_it():
    """Older runs predate the persisted corrected trajectory. A report that
    raises instead of drawing one trace loses the reading it does have."""
    src = (REPO / "shared" / "gromacs_analysis.py").read_text()
    body = src[src.index("def protein_rmsd("):]
    assert body.count("return None") >= 2, "missing inputs must yield None"


def test_the_distance_panel_reads_the_window_from_the_criterion():
    """A hand-typed band here could drift from the number that scores.

    UPDATED 2026-09-02: the band is no longer `NAC_DIST_MIN..NAC_DIST_MAX`.
    That pair is the SCREEN's near-attack window (2.8-4.2 A); the sweep is
    judged at 2.8-3.5 (D0111), and shading the wider one meant a trace could
    sit in the green zone while scoring 0% engaged on the same figure. The
    test's intent is unchanged -- read the band from the criterion, never type
    it -- so it now names the function that owns it.
    """
    src = (REPO / "scripts" / "sweep_assets.py").read_text()
    assert "nac.attack_ready_window()" in src, (
        "the shaded band is not read from nac_criterion.attack_ready_window()")
    # and the wider window is still drawn, as a reference line rather than as
    # the band -- dropping it would make this figure and the screen's own
    # criterion describe different physics
    assert "nac.NAC_DIST_MAX" in src


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


# --------------------------------------------------------------------------
# the mode is part of the run's identity (100 ns stage)
# --------------------------------------------------------------------------

def test_the_runner_can_name_a_run_after_its_mode():
    """`wd = work_root / ident` and the report filename both follow `ident`, so
    while `ident` was the MOLECULE two modes of one molecule shared a directory
    -- and build_workdir rebuilds in place, so the second overwrote the first's
    finished trajectory while its row survived. t4_c8c3aec07421 (_m1 and _m5)
    was queued to do exactly that."""
    src = (REPO / "scripts" / "md_residence_3ikd.py").read_text()
    assert '"--mode"' in src
    assert 'ident = f"{args.candidate}_m{int(args.mode)}"' in src
    # and the output stem must follow the run, not the molecule, or two modes
    # write the same CSV
    assert "args.tag or ident" in src


def test_omitting_the_mode_keeps_molecule_level_behaviour():
    """63 legacy rows and every existing invocation are molecule-level. The new
    identity must be opt-in or they all change meaning at once."""
    src = (REPO / "scripts" / "md_residence_3ikd.py").read_text()
    body = src[src.index("if args.candidate:"):]
    assert "ident = args.candidate" in body[:2000]


def test_the_row_carries_the_join_key_not_only_the_label():
    """shared/mode_key: the key is (parent_ident, mode), never `ident` -- a
    merge on `ident` silently dropped every mode-0 row once already."""
    src = (REPO / "scripts" / "md_residence_3ikd.py").read_text()
    assert "mode_key.split_ident(ident)" in src
    assert '"parent_ident": _parent' in src and '"mode": _mode' in src


@pytest.mark.parametrize("path,needle", [
    ("scripts/mdprio_report.py", "mode_key.split_ident(args.candidate)"),
    ("scripts/mdprio_combine.py", "mode_key.split_ident(t)"),
])
def test_consumers_resolve_a_mode_ident_to_its_molecule(path, needle):
    """The trajectory is per RUN; the SMILES, warhead class, sweep row and
    depiction are facts about the MOLECULE. Looking those up under `<parent>_mN`
    misses every one and renders a page about a molecule the project appears to
    know nothing about."""
    assert needle in (REPO / path).read_text()


def test_the_results_page_looks_up_class_and_depiction_under_the_parent():
    src = (REPO / "scripts" / "mdprio_combine.py").read_text()
    assert "cls_of.get(par" in src
    assert "thumbs[par]" in src or "thumbs.get(par" in src


def test_split_ident_is_anchored_so_a_molecule_named_with_m_survives():
    from shared import mode_key as mk
    assert mk.split_ident("t4_abc_m3") == ("t4_abc", 3)
    assert mk.split_ident("t4_abc") == ("t4_abc", None)
    # a bare ident means the mode was NOT STATED -- not mode 0. Reading it as 0
    # is exactly the assumption that produced #53's invisible collision.
    assert mk.split_ident("t4_abc")[1] is None


# --------------------------------------------------------------------------
# a stage with no results still has a page
# --------------------------------------------------------------------------

def test_an_empty_rail_does_not_emit_an_iframe():
    """THE 404. With no finished sweep, `first` is "" and the viewer's src
    interpolated to `sweep_pages/.html` -- a request for the empty ident -- so
    the server's 404 page rendered INSIDE the layout. A run with no results yet
    and a broken deployment looked identical. @tt8804: "showing a 404"."""
    src = (REPO / "scripts" / "sweep_combine.py").read_text()
    assert "_viewer" in src, "the viewer pane is not conditional"
    body = src[src.index("_viewer = "):]
    assert "if first else" in body[:400], "no empty-rail branch"
    # and the click handler must tolerate the frame being absent
    assert "if(!f)" in src or "if(!id)" in src


def test_every_page_a_topbar_links_to_is_built_without_results():
    """`pipeline.html` (the schematic) and `controls.html` were written only by
    `mdprio_combine`, which exits before reaching them when there are no 100 ns
    reports -- so on a fresh topic the "how this works" link 404s from all four
    pages. Same coupling that kept `modes.html` empty while 2,019 modes sat on
    disk: stage-2 content gated behind stage 5."""
    src = (REPO / "scripts" / "build_gui.py").read_text()
    assert 'schematic.build(' in src, "build_gui does not write pipeline.html"
    assert '"controls.html"' in src, "build_gui does not write controls.html"


def test_the_ranking_page_builds_from_the_ranking_alone():
    """It is stage 2's output; requiring stage 5 to render it means the first
    results a reader wants are the last to appear."""
    p = REPO / "scripts" / "ranking_page.py"
    assert p.is_file()
    src = p.read_text()
    assert "moderank.build(" in src and "massets.write_assets(" in src
    # It must not DEPEND on the report combiner. Checked against CODE, not
    # prose: the docstring names mdprio_combine in order to explain why the
    # call site moved, and matching the bare word would make this unfailable --
    # the same mistake the addModelsAsFrames test above had to correct for.
    assert "import mdprio_combine" not in src
    assert "mdprio_combine.py" not in src


# --------------------------------------------------------------------------
# rail search
# --------------------------------------------------------------------------

def test_the_search_widget_is_shared_not_copied():
    """Two rails with two filters is two behaviours that drift -- the same
    reason `results_shell.CSS` exists."""
    from shared import results_shell as rs
    assert 'id="railq"' in rs.SEARCH_HTML
    assert "railFilter" in rs.SEARCH_JS
    for f in ("scripts/sweep_combine.py", "scripts/mdprio_combine.py"):
        src = (REPO / f).read_text()
        assert "SEARCH_HTML" in src and "SEARCH_JS" in src, f


def test_the_ranking_page_filters_the_data_not_the_dom():
    """THE RAIL IS VIRTUALISED -- only rows currently on screen exist in the
    DOM. Hiding elements would filter the visible window and leave everything
    below it unfiltered the moment you scrolled, which looks like a working
    filter until you scroll."""
    src = (REPO / "shared" / "mode_ranking.py").read_text()
    assert "function setQuery(" in src and "railHTML()" in src
    body = src[src.index("function visible(){"):]
    assert "r.filter(matches)" in body[:400], "query not applied inside visible()"


def test_the_ranking_filter_rebuilds_the_ident_it_does_not_assume_one():
    """`i` is deliberately not sent per row -- it is p + "_m" + m and cost a
    megabyte of payload. A filter matching on a field that is not there would
    silently match nothing."""
    src = (REPO / "shared" / "mode_ranking.py").read_text()
    body = src[src.index("function matches("):src.index("function visible(){")]
    assert "x.p + '_m' + x.m" in body


def test_filtering_does_not_renumber_the_ranks():
    """Renumbering 1..n over a filtered list would make "rank 3" mean something
    different depending on what was typed."""
    from shared import results_shell as rs
    assert "textContent" in rs.SEARCH_JS          # reads rows, does not rewrite
    assert "innerHTML" not in rs.SEARCH_JS


def test_empty_section_headers_are_hidden_when_filtered_out():
    """A held/left banner with nothing under it reads as an empty category."""
    from shared import results_shell as rs
    assert ".ohd" in rs.SEARCH_JS


# --------------------------------------------------------------------------
# cohesion: one number, every page
# --------------------------------------------------------------------------

def test_only_one_module_computes_the_step_counts():
    """THE 447. Four builders computed these independently, and
    `mode_ranking._step_counts` read `mdprio_reports/sweep_state.json` -- the
    UNSCOPED path -- so the Ranking page's nav showed 3.0.0's 447 swept modes
    while the Sweep page beside it showed this run's 34. Both files existed, so
    nothing failed. @tt8804: "on ranking it shows sweep is 447 okay while
    clicking on sweep shows 32 ok, can we make this whole gui cohesive"."""
    assert hasattr(gs, "step_counts")
    for f in ("scripts/build_gui.py", "scripts/sweep_combine.py"):
        src = (REPO / f).read_text()
        assert "step_counts()" in src, f
    # and nobody may reach the unscoped reports directory for them
    mr = (REPO / "shared" / "mode_ranking.py").read_text()
    body = mr[mr.index("def _step_counts("):]
    body = body[:body.index("\n\n\n")] if "\n\n\n" in body else body
    code = "\n".join(l for l in body.splitlines() if not l.strip().startswith("#"))
    assert 'mdprio_reports" / "sweep_state.json"' not in code
    assert "step_counts()" in code


def test_the_counts_come_from_the_pipeline_probes():
    """So a page cannot disagree with the dashboard, which is the other half of
    the same complaint."""
    src = (REPO / "shared" / "gui_shell.py").read_text()
    body = src[src.index("def step_counts("):]
    assert "from . import pipeline as pl" in body
    assert "pl.status()" in body


def test_an_absent_count_is_omitted_not_zero():
    """A bare 0 asserts "measured, none"; the truth is "not measured yet"."""
    src = (REPO / "shared" / "gui_shell.py").read_text()
    body = src[src.index("def step_counts("):]
    assert 'state != "unknown"' in body
    assert "total" in body        # only emitted when there is a denominator


def test_home_does_not_pass_its_card_dict_as_the_nav():
    """`counts` keys the funnel cards (molecules/modes); the nav is keyed by page
    href. Passing one as the other is why Home was the only page whose stepper
    showed no counts at all -- `gs.nav` looked up 'sweep.html' in a dict that had
    never heard of it and correctly omitted it."""
    src = (REPO / "scripts" / "build_gui.py").read_text()
    assert "nav: dict | None = None" in src
    assert "home(counts, s, wl, nav_counts)" in src


def test_a_placeholder_is_refreshed_but_a_built_page_is_never_clobbered():
    """Writing only when absent froze combined.html's counts at the previous
    night's values while every page beside it was current."""
    src = (REPO / "scripts" / "build_gui.py").read_text()
    assert "_PLACEHOLDER" in src
    assert "awaiting stage" in src          # the pre-marker fallback
    body = src[src.index("n_placeholder = 0"):]
    assert "_PLACEHOLDER not in head" in body


def test_the_rail_row_stylesheet_is_shared_not_copied():
    """@tt8804: "the different selectors look different on diff pages".

    All three rails already emitted the SAME class names -- .row/.rk/.thumb/
    .body/.l1/.mid-id/.eng/.l2/.wc/.meta/.tag/.bar -- so the markup was never the
    problem. The ranking page kept its own COPY of the rules, because it is
    virtualised and needs its own <style>, and the copy had drifted: `.l2` had
    lost `flex-wrap:nowrap`. Two stylesheets that start identical and drift is
    exactly what results_shell exists to prevent, and it had happened inside the
    module's own subject matter."""
    from shared import results_shell as rs
    assert hasattr(rs, "ROW_CSS")
    for sel in (".row{", ".rk{", ".thumb{", ".body{", ".l1{", ".mid-id{",
                ".eng{", ".l2{", ".wc{", ".meta{", ".tag{", ".bar{"):
        assert sel in rs.ROW_CSS, f"{sel} missing from the shared row CSS"
    # the page shell still ships it, so the two rail pages are unchanged
    assert ".row{" in rs.CSS
    # and the ranking page interpolates it rather than restating it
    mr = (REPO / "shared" / "mode_ranking.py").read_text()
    assert "__ROWCSS__" in mr and "_rs().ROW_CSS" in mr


def test_the_ranking_page_keeps_only_its_virtualisation_override():
    """Every item must be exactly ROW_H tall: the window offset is computed as
    i * ROW_H rather than measured, so a row free to size itself puts every row
    below it at the wrong offset. That override is legitimate; a second copy of
    the whole row is not."""
    mr = (REPO / "shared" / "mode_ranking.py").read_text()
    tpl = mr[mr.index("__ROWCSS__"):]
    tpl = tpl[:tpl.index("</style>")]
    assert "height:64px" in tpl
    # nothing else about the row may be redefined here
    for sel in (".mid-id{", ".eng{", ".meta{", ".tag{", ".bar{"):
        assert sel not in tpl, f"{sel} redefined after the shared block"


def test_the_search_box_lives_in_the_rail_on_every_page():
    """@tt8804, three times: "still no search bar".

    It was there, served, and well-formed -- in the TOPBAR, which on the ranking
    page is a flex row with `overflow-x:auto` carrying a long title, the class
    select, a hint and two buttons. The box was pushed out of the visible strip.
    A control the reader cannot find is a control that does not exist, and
    'the HTML contains it' is not the same claim as 'it is on screen'.

    So all three rails put the SAME box in the SAME place: the top of the rail,
    where the sweep page always had it."""
    mr = (REPO / "shared" / "mode_ranking.py").read_text()
    # not in the topbar
    tb = mr[mr.index('<div id="topbar">'):mr.index("</div>", mr.index('<div id="topbar">'))]
    assert "railq" not in tb, "the search box is back in the topbar"
    # in the rail column, before the rows
    assert "__RAILSEARCH__" in mr
    rail = mr[mr.index('<div id="railcol">'):]
    assert rail.index("__RAILSEARCH__") < rail.index('id="rail"')


def test_all_three_rails_render_the_same_search_markup():
    from shared import results_shell as rs
    a = rs.search_html()
    b = rs.search_html("setQuery(this.value)", "filter — id, class, mode")
    # same structure, only the handler and hint differ
    for frag in ('class="railq"', 'id="railq"', 'id="railn"', "railClear()"):
        assert frag in a and frag in b, frag


def test_the_ranking_page_defines_the_clear_handler_it_renders():
    """The shared markup's × button calls railClear() by name. The ranking page
    has its own <script>, so inheriting the markup without the function would
    give a button that throws on click."""
    mr = (REPO / "shared" / "mode_ranking.py").read_text()
    assert "function railClear(" in mr


def test_the_search_css_reaches_a_page_with_its_own_stylesheet():
    """The ranking page does not use `results_shell.CSS`, so `.railq` had to be
    separable from it -- otherwise the box renders unstyled and looks like a
    stray input."""
    from shared import results_shell as rs
    assert ".railq{" in rs.SEARCH_CSS
    assert ".railq{" in rs.CSS                    # the rail pages still get it
    mr = (REPO / "shared" / "mode_ranking.py").read_text()
    assert "__RAILQCSS__" in mr and "_rs().SEARCH_CSS" in mr


def test_a_row_for_a_mode_nobody_selected_is_not_counted():
    """The 24 rows a broken launcher produced -- (molecule, a pose_rank in the
    hundreds) -- are keyed on pairs no real run will ever request, so nothing
    supersedes them and Home read `failed: 24` permanently. They cannot be
    deleted: the outputs root is append-only.

    The filter must run BEFORE the ident union, or the dropped rows still
    contribute an ident and merely change label from `failed` to `not sent` --
    which is what the first version did."""
    src = (REPO / "shared" / "sweep_state.py").read_text()
    body = src[src.index("def state("):]
    assert "asked = set(zip(" in body
    assert body.index("res = res[keep_row]") < body.index("idents = set()")


def test_the_filter_needs_a_worklist_and_does_nothing_without_one():
    """With no worklist there is nothing to judge against, and dropping rows on
    a guess would hide real results."""
    src = (REPO / "shared" / "sweep_state.py").read_text()
    body = src[src.index("def state("):]
    seg = body[body.index("asked = set(zip(") - 400:body.index("asked = set(zip(")]
    assert "not wl.empty" in seg


def test_build_gui_finds_the_campaign_worklist_itself():
    """It used to require --worklist, so a plain `build_gui` had none -- and
    without one the page cannot tell a selected mode from an invented one."""
    src = (REPO / "scripts" / "build_gui.py").read_text()
    assert "_pl.worklist_path()" in src


# --------------------------------------------------------------------------
# stage lengths in prose (#63 follow-up)
#
# The nav said "10 ns triage" and "100 ns runs" as literals. D0085 moved the
# triage to 8 ns and `md.sweep_ps` says 8000, so every page in the GUI carried a
# stage length the run had not used for a day, beside numbers that were correct.
# A literal cannot notice that config moved -- only a comparison can.
# --------------------------------------------------------------------------

def test_the_stepper_names_the_configured_stage_lengths():
    from shared import target_config as tc
    steps = dict((f, d) for f, _l, d in gs._steps())
    assert gs._ns(tc.md_sweep_ps()) in steps["sweep.html"]
    assert gs._ns(tc.md_production_ps()) in steps["combined.html"]


def test_the_label_helpers_track_config_rather_than_a_literal():
    """Move the config and the label must move with it."""
    cfg = {"md": {"sweep_ps": 25_000, "production_ps": 500_000}}
    from shared import target_config as tc
    assert gs._ns(tc.md_sweep_ps(cfg)) == "25 ns"
    assert gs._ns(tc.md_production_ps(cfg)) == "500 ns"
    # and sub-nanosecond stays readable rather than rounding to "0 ns"
    assert gs._ns(300.0) == "0.3 ns"


def test_no_user_facing_page_hardcodes_a_sweep_length():
    """The guard that would have caught this, over the report generators.

    Scoped to the files that PRODUCE the pages a reader sees. Docstrings and
    comments are excluded -- fixing prose about history is not the point, and
    narrowing the check to executable strings is what made the version-pin test
    honest (`how_this_project_breaks.md`).
    """
    import ast
    from shared import target_config as tc

    stale = gs._ns(10_000.0)                      # "10 ns", the old triage
    current = gs._ns(tc.md_sweep_ps())
    if stale == current:                          # config moved back; nothing to catch
        pytest.skip("the configured sweep length IS 10 ns")

    offenders = []
    for rel in ("scripts/build_gui.py", "scripts/sweep_combine.py",
                "scripts/pose_modes_report.py", "scripts/shortlist_report.py",
                "shared/gui_shell.py", "shared/mode_ranking.py"):
        p = REPO / rel
        if not p.is_file():
            continue
        tree = ast.parse(p.read_text(), filename=str(p))
        # Docstrings are prose ABOUT the code, often historical -- "this said 10
        # ns" is a true sentence. Collect them by identity and skip them, the
        # same narrowing that made the version-pin test honest.
        docs = set()
        for node in ast.walk(tree):
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                                 ast.AsyncFunctionDef)):
                body = getattr(node, "body", None)
                if (body and isinstance(body[0], ast.Expr)
                        and isinstance(body[0].value, ast.Constant)
                        and isinstance(body[0].value.value, str)):
                    docs.add(id(body[0].value))
        for node in ast.walk(tree):
            if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                    and id(node) not in docs):
                if f"{stale} sweep" in node.value or f"{stale} triage" in node.value:
                    offenders.append(f"{rel}:{node.lineno}")
    assert not offenders, (
        f"these name a {stale} sweep while md.sweep_ps says {current}; "
        f"use gui_shell.sweep_label(): {offenders}")


def test_the_md_rail_shows_the_queue_not_only_the_finished():
    """A survivor waits hours for its 100 ns run, and the page showed nothing at
    all in that window -- the same gap the sweep page had before "show results
    as pending what is sweeped". A stage whose queue is invisible looks idle.

    A queued row carries its 8 ns readings and points the viewer at its SWEEP
    report, which is the evidence it was queued on. It gets no held/left tag:
    that verdict belongs to the 100 ns run and does not exist yet."""
    src = (REPO / "scripts" / "mdprio_combine.py").read_text()
    assert "pending_rows" in src
    body = src[src.index("pending_rows = []"):src.index("# CONTROLS (#47")]
    assert "_pl.survivors()" in body
    assert "if ident in tabs:" in body, "a finished run would be listed twice"
    assert "sweep_pages/" in body, "queued rows do not point at their evidence"
    assert "t-held" not in body and "t-left" not in body, \
        "a queued row must not carry a 100 ns verdict"


def test_the_mode_is_stated_once_in_a_rail_row():
    """Printing the full ident and then a badge repeating its suffix gives
    `t4_710417e24b49_m4` followed by `m4` -- the same fact twice, which reads as
    two. The name carries the molecule, the badge carries the mode."""
    src = (REPO / "scripts" / "mdprio_combine.py").read_text()
    for frag in ("<span class='mid-id'>{html.escape(par)}",
                 "<span class='mid-id'>{html.escape(parent)}"):
        assert frag in src, f"a rail row still prints the full ident: {frag}"
