"""The ranking rail must cost the viewport, not the library.

WHY THIS IS A TEST AND NOT A NOTE. The rail rendered every visible mode as a
<button> holding an <img> and eight <span>s, in one innerHTML assignment, and
`pick()` called it -- so every row click and every mode checkbox rebuilt the
whole list. At 8,097 modes that was merely slow. Sub-splitting (#61) took the
library to 34,076, roughly 350,000 nodes and 34,000 lazy images, and @tt8804
reported the page as "extremely laggy".

The defect is structural, not a constant factor: the view scaled with the
library instead of with the screen. Nothing about the page LOOKS wrong when it
regresses -- it renders correctly and simply becomes unusable as the library
grows -- so the invariants are asserted here rather than left to be noticed
again at the next size increase.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import mode_ranking as mr                     # noqa: E402

TPL = mr._TPL


def test_the_rail_renders_a_window_not_the_whole_list():
    """A bounded slice, chosen from scroll position."""
    assert "function renderRail" in TPL
    body = TPL[TPL.index("function renderRail"):]
    body = body[:body.index("\nfunction ")]
    assert "scrollTop" in body, "the window is not derived from scroll position"
    assert re.search(r"for \(let i = s; i < e; i\+\+\)", body), \
        "renderRail no longer iterates a bounded [s, e) slice"


def test_selection_does_not_rebuild_the_rail():
    """`pick()` moves a class. It must not re-render 34,076 rows."""
    pick = TPL[TPL.index("async function pick("):]
    pick = pick[:pick.index("\nfunction ")]
    assert "markSel()" in pick
    assert "railHTML()" not in pick, \
        "pick() rebuilds the whole rail again -- this is the reported lag"


def test_the_scroll_handler_is_frame_throttled():
    """Scroll fires far more often than the screen refreshes."""
    assert "requestAnimationFrame" in TPL
    m = re.search(r"addEventListener\('scroll'.*?\}, \{passive: true\}\)", TPL, re.S)
    assert m, "no passive throttled scroll listener"
    assert "requestAnimationFrame" in m.group(0), "scroll renders per event"


def test_the_scrollbar_reflects_the_whole_library():
    """A windowed list still has to be honest about how much there is."""
    assert "railPad" in TPL and "TOTAL" in TPL
    build = TPL[TPL.index("function buildItems"):]
    build = build[:build.index("\nfunction ")]
    assert "railPad" in build and "TOTAL + 'px'" in build


def test_row_geometry_is_fixed_so_the_window_is_computable():
    """Variable row heights drift the list out of register with the scrollbar.

    Every item is a row of one height -- there are no inline class headers, so
    the offset of item i is exactly i * ROW_H. The CSS and the JS must agree on
    that number or the list slides away from the scrollbar, which looks like a
    rendering bug rather than a units mismatch.
    """
    css_row = int(re.search(r"\.row\{height:(\d+)px", TPL).group(1))
    js_row = int(re.search(r"const ROW_H = (\d+)", TPL).group(1))
    assert css_row == js_row, \
        "CSS and JS disagree about row height -- the window will drift"


def test_the_class_label_appears_once():
    """One label, not a banner plus an inline header saying the same word.

    @tt8804: "says acrylamide twice". The sticky header could not survive a
    transformed window so a banner replaced it -- and the inline header was left
    in, so both rendered, one above the other.
    """
    render = TPL[TPL.index("function renderRail"):]
    render = render[:render.index("\nfunction ")]
    assert "railBanner" in render, "nothing names the class being viewed"
    assert "chd" not in render, "inline class headers are back alongside the banner"
    build = TPL[TPL.index("function buildItems"):]
    build = build[:build.index("\nfunction ")]
    assert "{h:" not in build and "chd" not in build, \
        "buildItems emits header items again"


def test_a_mode_with_no_pose_in_the_asset_is_reported():
    """An empty viewer is indistinguishable from a broken one."""
    assert "vmiss" in TPL
    draw = TPL[TPL.index("function draw("):]
    assert "if (sel < 0)" in draw, "a missing model is not reported"
    assert 'id="vmiss"' in TPL, "the element the message goes into does not exist"


def test_the_pose_error_reports_the_actual_reason():
    """It used to blame a missing file for every failure, including its own."""
    assert "No pose drawn" in TPL
    # Comment lines are stripped: the prose above the fix names the old string in
    # order to explain it, and matching that would make this test unfailable.
    code = "\n".join(ln for ln in TPL.splitlines()
                     if not ln.lstrip().startswith("//"))
    assert "no pose file for " not in code, \
        "the catch still reports every error as a missing file"


def test_the_sort_is_not_done_per_frame():
    """`visible()` slices and sorts 34,076 rows; that belongs to scope changes."""
    render = TPL[TPL.index("function renderRail"):]
    render = render[:render.index("\nfunction ")]
    assert "visible()" not in render, "renderRail sorts the library every frame"


def test_main_has_exactly_as_many_children_as_it_has_grid_columns():
    """A stray child silently relocates the viewer.

    `main` is a two-column grid. Adding the class banner as a THIRD child put it
    in column 1 and pushed the viewer onto row 2 -- @tt8804: "why is the viewer
    on the left". Grid does not error on overflow, it wraps, so this is invisible
    until someone looks at the page.
    """
    body = TPL[TPL.index("<main>"):TPL.index("</main>")]
    depth, children = 0, 0
    for m in re.finditer(r"<(/?)(\w+)[^>]*?(/?)>", body):
        closing, tag, selfclose = m.group(1), m.group(2), m.group(3)
        if tag == "main":
            continue
        if closing:
            depth -= 1
        elif not selfclose and tag not in ("br", "img", "input", "hr"):
            if depth == 0:
                children += 1
            depth += 1
        elif depth == 0:
            children += 1
    cols = re.search(r"main\{[^}]*grid-template-columns:([^;}]+)", TPL).group(1)
    n_cols = len(cols.split())
    assert children == n_cols, (
        f"<main> has {children} direct children for {n_cols} grid columns; "
        "the extra one will wrap onto a new row and move the viewer")


def test_size_containment_is_never_used_on_the_rail():
    """`contain:strict` collapses the rail to a strip."""
    m = re.search(r"#rail\{[^}]*\}", TPL).group(0)
    assert "contain:strict" not in m, \
        "size containment makes the rail contribute zero height"
    assert "contain:layout paint" in m


def test_the_banner_lives_inside_the_rail_column():
    rc = TPL[TPL.index('<div id="railcol">'):]
    rc = rc[:rc.index("</div>\n </div>") + 14] if "</div>\n </div>" in rc else rc[:800]
    assert 'id="railBanner"' in rc and 'id="rail"' in rc


def test_the_ident_is_not_sent_on_the_wire_but_exists_in_js():
    """Derivable fields are rebuilt on load, not repeated 34,076 times."""
    assert "x.i = x.p + '_m' + x.m" in TPL
    src = (REPO / "shared" / "mode_ranking.py").read_text()
    body = src[src.index("def _rows_json"):src.index("_TPL = r")]
    assert '"i": str(x.ident)' not in body, "ident is back on the wire"
