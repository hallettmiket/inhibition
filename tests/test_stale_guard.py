"""
Purpose: the stale-module guard must run BEFORE any helper attribute is touched.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: integration/app/app.py source
Output: pass/fail

WHY THIS EXISTS. Streamlit re-runs app.py on every interaction but does not
re-import local helper modules. After curate.py gained PANEL_SCOPE, a running
process executed the NEW app.py against the OLD curate module and died with

    AttributeError: module 'curate' has no attribute 'PANEL_SCOPE'

The app already had a stale-module detector whose entire purpose is to say
"restart me" in exactly this situation. It sat at the BOTTOM of the file, ~30
lines AFTER the attribute access that crashed. A diagnostic that only runs once
the crash it diagnoses has already happened is not a diagnostic.

The guard is easy to move and easy to defeat by adding a module-level helper
access above it, so the ordering is pinned here rather than trusted.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

APP = Path(__file__).resolve().parent.parent / "integration" / "app" / "app.py"
HELPERS = {"D", "depict", "curate", "p3d", "syn"}


@pytest.fixture(scope="module")
def source() -> str:
    return APP.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def tree(source: str) -> ast.Module:
    return ast.parse(source)


def _guard_stop_line(tree: ast.Module) -> int:
    """Line of the `st.stop()` that ends the stale guard."""
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "stop"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "st"):
            return node.lineno
    pytest.fail("no `st.stop()` found — the stale guard must halt the page, "
                "not merely warn, or the AttributeError still reaches the user")


def test_the_guard_exists(tree):
    assert _guard_stop_line(tree) > 0


def test_stale_modules_helper_is_defined(tree):
    names = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    assert "stale_modules" in names, (
        "stale_modules() should be a named function, not inlined — it is "
        "exercised directly by tests and read by humans debugging a restart")


def _module_level_helper_attrs(tree: ast.Module) -> list[tuple[int, str]]:
    """(line, 'mod.attr') for helper attribute reads at MODULE level only.

    Accesses inside a function body do not run at import, so they cannot be the
    thing that crashes before the guard. Only top-level statements matter.

    BOTH forms are counted. `curate.PANEL_SCOPE` raises on a stale module;
    `getattr(curate, "PANEL_SCOPE", ())` does not, but it still READS the
    module and is still the thing whose position matters. Counting only the
    dotted form made this check vacuous the moment the crashing line was made
    defensive -- zero hits, guaranteed pass, no coverage.
    """
    hits = []
    for stmt in tree.body:
        if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            continue
        for node in ast.walk(stmt):
            if (isinstance(node, ast.Attribute)
                    and isinstance(node.value, ast.Name)
                    and node.value.id in HELPERS):
                hits.append((node.lineno, f"{node.value.id}.{node.attr}"))
            elif (isinstance(node, ast.Call)
                  and isinstance(node.func, ast.Name)
                  and node.func.id == "getattr"
                  and node.args
                  and isinstance(node.args[0], ast.Name)
                  and node.args[0].id in HELPERS):
                attr = (node.args[1].value
                        if len(node.args) > 1
                        and isinstance(node.args[1], ast.Constant)
                        else "?")
                hits.append((node.lineno,
                             f"getattr({node.args[0].id}, {attr!r})"))
    return sorted(hits)


def test_the_ordering_check_is_not_vacuous(tree):
    """A test that cannot fail is not a test.

    If nothing at module level reads a helper, `test_no_helper_attribute_is_
    read_before_the_guard` passes for free and would keep passing after
    someone reintroduced the defect in a form it does not recognise.
    """
    hits = _module_level_helper_attrs(tree)
    assert hits, (
        "no module-level helper access found at all — either app.py changed "
        "shape or this detector no longer recognises the access form it is "
        "meant to police. Fix the detector, do not delete the test.")


def test_no_helper_attribute_is_read_before_the_guard(tree):
    """The ordering defect itself, pinned.

    Anything reading `curate.X` / `D.X` / `p3d.X` at module level above the
    guard can raise AttributeError on a stale process before the guard has a
    chance to explain why.
    """
    stop = _guard_stop_line(tree)
    early = [(ln, name) for ln, name in _module_level_helper_attrs(tree)
             if ln < stop]
    assert not early, (
        f"helper attributes read at module level BEFORE the stale guard "
        f"(line {stop}): {early}. On a stale process these raise "
        "AttributeError and the user sees a traceback instead of 'restart "
        "me'. Move them below the guard.")


def test_guard_runs_before_the_panel_scope_check(source: str):
    """The specific access that crashed, kept behind the guard."""
    stop_m = re.search(r"^\s*st\.stop\(\)", source, re.M)
    scope_m = re.search(r"getattr\(curate,\s*[\"']PANEL_SCOPE[\"']|"
                        r"curate\.PANEL_SCOPE", source)
    assert stop_m and scope_m
    assert stop_m.start() < scope_m.start(), (
        "the PANEL_SCOPE check must sit below the stale guard — it is the "
        "exact access that produced the reported AttributeError")


def test_panel_scope_access_is_defensive(source: str):
    """Belt and braces: even reached on a stale module, it must not raise."""
    assert re.search(r"getattr\(curate,\s*[\"']PANEL_SCOPE[\"']", source), (
        "read PANEL_SCOPE via getattr with a default — a bare "
        "curate.PANEL_SCOPE reintroduces the crash for any process that "
        "somehow gets past the guard")


def test_stale_modules_reports_a_changed_file(tmp_path, monkeypatch):
    """Behavioural check on the helper itself, not just its presence."""
    import importlib.util
    import types

    mod_path = tmp_path / "fake_helper.py"
    mod_path.write_text("LOADED_MTIME = 0.0\n", encoding="utf-8")
    spec = importlib.util.spec_from_file_location("fake_helper", mod_path)
    fake = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fake)

    # LOADED_MTIME of 0 is far in the past, so the file is "changed".
    assert isinstance(fake, types.ModuleType)
    assert Path(fake.__file__).stat().st_mtime > fake.LOADED_MTIME + 1

    # And a module whose recorded mtime is current is NOT stale.
    fake.LOADED_MTIME = Path(fake.__file__).stat().st_mtime
    assert not Path(fake.__file__).stat().st_mtime > fake.LOADED_MTIME + 1
