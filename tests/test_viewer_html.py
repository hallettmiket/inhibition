"""A blank 3Dmol viewer has failed four distinct ways in this project in one day.

Every one produced the SAME symptom — an empty box, no console error a reader
would see, a page that otherwise rendered perfectly — and every one was found by
a person looking at a screenshot rather than by a test:

1. `render()` never called. 3Dmol draws nothing without it; labels are DOM
   overlays and appear anyway, so the panel looked like a working viewer with the
   molecule missing.
2. Zero-height container. `report_theme` sizes `.glbox` but has no rule for its
   child, and a viewer built into a 0-height div draws nothing. Made worse by
   panels that start closed: a closed `<details>` has no height at all.
3. Library after use. 3Dmol emitted at the end of `<body>`, so a module-level
   `window.$3Dmol` lookup captured `undefined`.
4. Data element after use. `<script type="text/plain" id="rec">` emitted below the
   script that reads it, so `getElementById` returned null and the block threw on
   its first line.

These are all orderings and omissions, which means they are all statically
checkable. This file checks them on the HTML each builder actually produces, so
the fifth one fails a test instead of a screenshot.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

#: Builders that emit a 3Dmol viewer, and the marker that says they did.
BUILDERS = [
    "scripts/crystal_recovery_report.py",
    "scripts/pose_modes_report.py",
    "scripts/shortlist_report.py",
    "shared/mode_ranking.py",
    "shared/md_movie.py",
]


def _src(rel: str) -> str:
    p = REPO / rel
    if not p.is_file():
        pytest.skip(f"{rel} not present")
    return p.read_text()


@pytest.mark.parametrize("rel", BUILDERS)
def test_every_viewer_calls_render(rel):
    """3Dmol draws nothing until render() is called."""
    s = _src(rel)
    if "createViewer" not in s:
        pytest.skip("no viewer in this builder")
    assert re.search(r"\.render\(\)", s), (
        f"{rel} creates a 3Dmol viewer and never calls render(). The box will be "
        f"empty and nothing will say so.")


@pytest.mark.parametrize("rel", BUILDERS)
def test_the_viewer_container_is_given_a_size(rel):
    """A viewer built into a zero-height element draws nothing."""
    s = _src(rel)
    if "createViewer" not in s:
        pytest.skip("no viewer in this builder")
    sized = ("inset:0" in s or "inset: 0" in s
             or re.search(r"\.glbox\s*>\s*div", s)
             or re.search(r"height:\s*\d+px", s))
    assert sized, (
        f"{rel} builds a viewer without giving its container an explicit size. "
        f"report_theme sizes .glbox but not its child.")


@pytest.mark.parametrize("rel", BUILDERS)
def test_the_library_is_resolved_late_or_declared_first(rel):
    """`window.$3Dmol` captured at module scope is undefined if the library is
    emitted later in the document. Either resolve it inside the boot function,
    or emit the library in the <head>."""
    s = _src(rel)
    if "createViewer" not in s:
        pytest.skip("no viewer in this builder")
    late_bound = re.search(r"(function|=>)[^{]*\{[^}]*window\.\$3Dmol", s, re.S)
    in_head = re.search(r"<script>\{?three\}?</script>\s*(<style|</head)", s)
    head_slot = "__THREE__</script>" in s.replace("{three}", "__THREE__") and \
                s.find("</head>") > s.find("__THREE__") if "__THREE__" in s else False
    assert late_bound or in_head or head_slot or "vendored in the HEAD" in s, (
        f"{rel} resolves $3Dmol at module scope and may emit the library after "
        f"it. Resolve it inside boot(), or vendor it in the <head>.")


def test_data_elements_are_declared_before_they_are_read():
    """`getElementById(x)` returns null if <… id=x> is emitted below the script.

    Checked on a BUILT page rather than the source, because this is a property of
    the emitted document.
    """
    import glob
    pages = sorted(glob.glob(
        "/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
        "crystal_recovery/crystal_recovery_*.html"))
    if not pages:
        pytest.skip("no built crystal_recovery page to check")
    h = Path(pages[-1]).read_text()
    for eid in re.findall(r"getElementById\(['\"]([\w-]+)['\"]\)", h):
        decl = h.find(f'id="{eid}"')
        if decl < 0:
            decl = h.find(f"id='{eid}'")
        if decl < 0:
            continue                      # created at runtime, not in the markup
        use = h.find(f"getElementById('{eid}')")
        if use < 0:
            use = h.find(f'getElementById("{eid}")')
        # A DOMContentLoaded wrapper makes order irrelevant; without one it is not.
        if "DOMContentLoaded" in h[:use] or "addEventListener('load'" in h[:use]:
            continue
        assert decl < use, (
            f"element #{eid} is declared at {decl} and read at {use}: the script "
            f"runs before the element exists and will throw on that line.")
