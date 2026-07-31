"""
Purpose: The curation filter must reach every candidate panel — and only those.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: integration/app/{app.py, curate.py}
Output: pass/fail

WHY THIS EXISTS. Issue #3.2: "the curation feature does not carry over through
the dossier and rest of gui". It was real and it was worse than a missing
feature. The constraint box lived inside the Shortlists panel, wrote the spec to
`st.session_state["_curate_spec"]`, and nothing else ever read it — so a chemist
who excluded chlorines was TOLD the compounds were gone, then met them again in
the dossier, in the convergence pairs and in the axis medians. The panel that
applied the filter is the one that made the other panels misleading.

Two halves, and the second matters as much as the first:

  1. Every panel that shows candidates must honour the filter. Enforced by
     walking app.py's AST and requiring each such panel to route its frames
     through `curated()` — an assertion about the code, because a panel's body
     only runs when its tab is clicked and nothing else catches an omission.

  2. Nothing that reports on the FULL generated population may be filtered.
     `rank` denominators, the synthesizability rebuild's counts, the
     cross-approach convergence lookup and the gate statistics are facts about
     what the pipeline produced. Curating those would not yield a curated fact;
     it would yield a false one — "ranked 3rd of 1,204 docked" would silently
     become "3rd of 46".
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pandas as pd
import pytest

APP_DIR = Path(__file__).resolve().parent.parent / "integration" / "app"
sys.path.insert(0, str(APP_DIR))
import curate  # noqa: E402

APP_SRC = (APP_DIR / "app.py").read_text()
APP_TREE = ast.parse(APP_SRC)


def _panels_declared_in_app() -> dict[str, str]:
    """{panel title: function name} from app.py's PANELS dict, without importing.

    app.py imports streamlit, which is not installed in the test environment
    (dwi_cheminf) by design — the GUI env is separate. The AST is the honest way
    to ask what panels exist, and it is the approach tests/test_app_names.py
    already takes.
    """
    for node in ast.walk(APP_TREE):
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == "PANELS"
                        for t in node.targets)
                and isinstance(node.value, ast.Dict)):
            return {k.value: v.id for k, v in zip(node.value.keys, node.value.values)
                    if isinstance(k, ast.Constant) and isinstance(v, ast.Name)}
    raise AssertionError("app.py no longer defines a PANELS dict")


def _fn(name: str) -> ast.FunctionDef:
    for node in ast.walk(APP_TREE):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"app.py has no function {name}()")


def _calls(fn: ast.FunctionDef) -> set[str]:
    return {n.func.id for n in ast.walk(fn)
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}


# --- 1. every panel is declared, and the filtered ones actually filter ------

def test_every_panel_has_a_declared_scope():
    """A panel nobody decided about is how the original bug survived."""
    missing = set(_panels_declared_in_app()) - {s.panel for s in curate.PANEL_SCOPE}
    assert not missing, (
        f"panels with no curation scope declared: {sorted(missing)}. Add them "
        "to curate.PANEL_SCOPE with a reason — a candidate view must not "
        "default to unfiltered by omission.")


def test_no_scope_is_declared_for_a_panel_that_does_not_exist():
    stale = {s.panel for s in curate.PANEL_SCOPE} - set(_panels_declared_in_app())
    assert not stale, f"PANEL_SCOPE names panels that no longer exist: {sorted(stale)}"


def test_scope_for_refuses_an_undeclared_panel():
    with pytest.raises(KeyError):
        curate.scope_for("Some Panel Nobody Declared")


@pytest.mark.parametrize(
    "panel", [s.panel for s in curate.PANEL_SCOPE if s.filtered])
def test_a_filtered_panel_actually_calls_the_filter(panel):
    """THE REGRESSION GUARD for issue #3.2.

    Declaring a panel filtered is not the same as filtering it — that gap is
    exactly what shipped. This asserts the panel's own body routes its frames
    through `curated()`.
    """
    fn = _fn(_panels_declared_in_app()[panel])
    assert "curated" in _calls(fn), (
        f"{fn.name}() is declared filtered but never calls curated(); the "
        "curation filter would not carry over to it — the exact defect in "
        "issue #3.2")


@pytest.mark.parametrize("panel", [s.panel for s in curate.PANEL_SCOPE])
def test_every_panel_states_its_relationship_to_the_filter(panel):
    """Silence is what the bug looked like, so no panel may be silent."""
    fn = _fn(_panels_declared_in_app()[panel])
    assert "curation_header" in _calls(fn), (
        f"{fn.name}() never calls curation_header(); a reader cannot tell "
        "whether their filter is in force here")


@pytest.mark.parametrize(
    "panel", [s.panel for s in curate.PANEL_SCOPE if not s.filtered])
def test_an_unfiltered_panel_never_calls_the_filter(panel):
    fn = _fn(_panels_declared_in_app()[panel])
    assert "curated" not in _calls(fn), (
        f"{fn.name}() is declared UNFILTERED but calls curated()")


def test_the_constraint_box_is_not_inside_a_panel_any_more():
    """It lived in panel_candidates, which is why it could not carry over.

    The sidebar is the only thing rendered on every page.
    """
    assert "text_area" not in APP_SRC.split("def panel_candidates")[1].split(
        "\ndef ")[0], (
        "the constraint box is back inside panel_candidates; a filter that is "
        "invisible from the panels it is not affecting cannot be noticed to be "
        "off")
    assert "st.sidebar.text_area" in APP_SRC


def test_every_unfiltered_panel_explains_why():
    for scope in curate.PANEL_SCOPE:
        if not scope.filtered:
            assert len(scope.why) > 40, (
                f"{scope.panel} is excluded from the filter with no real "
                "reason recorded")


def test_the_facts_the_filter_must_never_touch_are_on_the_record():
    names = {n for n, _ in curate.UNFILTERED_FACTS}
    assert {"rank denominators", "enrichment gate verdicts"} <= names
    for _, why in curate.UNFILTERED_FACTS:
        assert len(why) > 40


# --- 2. the filter's own reporting -----------------------------------------

@pytest.fixture
def frame() -> pd.DataFrame:
    return pd.DataFrame({
        "candidate_id": ["a", "b", "c", "d"],
        "canonical_smiles": ["Clc1ccccc1", "c1ccccc1", "CCCCCCCCCCCCCCCC",
                             "O=S(=O)(N)c1ccccc1"],
        "rank": [1, 2, 3, 4],
    })


def test_banner_states_kept_and_removed_both(frame):
    _, rules = curate.apply(frame, "no chlorine")
    text = curate.banner(rules, 4, 3)
    assert "3 of 4" in text and "1 hidden" in text and "no chlorine" in text


def test_banner_with_no_rules_says_the_filter_is_off():
    assert "No curation filter active" in curate.banner([], 25, 25)


def test_describe_gives_every_rule_its_own_count(frame):
    _, rules = curate.apply(frame, "no chlorine\nheavy_atoms <= 7")
    text = curate.describe(rules)
    assert "`no chlorine` −1" in text and "`heavy_atoms <= 7` −2" in text


def test_bounded_axes_maps_a_constraint_to_the_axis_it_truncates(frame):
    _, rules = curate.apply(frame, "mw < 200\nheavy_atoms <= 8\nno chlorine")
    bounded = curate.bounded_axes(rules)
    assert bounded == {"MW": "mw < 200", "HAC": "heavy_atoms <= 8"}, (
        "a numeric bound truncates the shared axis it names, and the axes "
        "panel must refuse to present that median as a property of the approach")


def test_a_substructure_rule_bounds_no_axis(frame):
    _, rules = curate.apply(frame, "no chlorine")
    assert curate.bounded_axes(rules) == {}


def test_every_bounded_axis_is_a_real_shared_axis():
    """The mapping is only meaningful if the axis actually exists in the GUI."""
    sys.path.insert(0, str(APP_DIR.parent.parent))
    from integration.app import data as D

    for prop, axis in curate.PROPERTY_AXIS.items():
        assert prop in curate.PROPERTIES, f"{prop} is not a curatable property"
        assert axis in D.SHARED_AXES, f"{axis} is not a shared physicochemical axis"
