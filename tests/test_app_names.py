"""
Purpose: No panel in the GUI may reference a name it never binds.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: integration/app/*.py
Output: pass/fail

WHY THIS EXISTS. `panel_candidates` iterates `for col, (key, cfg) in ...` and a
block inside it referred to `a`, a name from a different function. Python only
raises on the line that runs, and a Streamlit panel runs only when someone
clicks its tab -- so the app imported cleanly, passed every test, and threw
NameError in front of a user.

Import alone cannot catch it and neither can a syntax check. This walks each
function's own scope and reports any load of a name that is neither bound
locally nor available globally. It is not a type checker; it catches exactly the
mistake that reached a user.
"""

from __future__ import annotations

import ast
import builtins
from pathlib import Path

import pytest

APP_DIR = Path(__file__).resolve().parent.parent / "integration" / "app"
APP_FILES = sorted(APP_DIR.glob("*.py"))


#: Always present at module scope, bound by the import machinery rather than by
#: any statement the AST can see.
MODULE_DUNDERS = {"__file__", "__name__", "__doc__", "__package__",
                  "__spec__", "__loader__", "__builtins__", "__path__"}


def _module_globals(tree: ast.Module) -> set[str]:
    names: set[str] = set(dir(builtins)) | MODULE_DUNDERS
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for t in node.targets:
                names |= {n.id for n in ast.walk(t) if isinstance(n, ast.Name)}
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            names |= {a.asname or a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    return names


def _bound_in(fn: ast.FunctionDef) -> set[str]:
    bound = {a.arg for a in fn.args.args} | {a.arg for a in fn.args.kwonlyargs}
    if fn.args.vararg:
        bound.add(fn.args.vararg.arg)
    if fn.args.kwarg:
        bound.add(fn.args.kwarg.arg)
    for n in ast.walk(fn):
        if isinstance(n, ast.Assign):
            for t in n.targets:
                bound |= {x.id for x in ast.walk(t) if isinstance(x, ast.Name)}
        elif isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name):
            bound.add(n.target.id)
        elif isinstance(n, (ast.For, ast.AsyncFor, ast.comprehension)):
            bound |= {x.id for x in ast.walk(n.target) if isinstance(x, ast.Name)}
        elif isinstance(n, ast.withitem) and n.optional_vars:
            bound |= {x.id for x in ast.walk(n.optional_vars)
                      if isinstance(x, ast.Name)}
        elif isinstance(n, ast.ExceptHandler) and n.name:
            bound.add(n.name)
        elif isinstance(n, (ast.Lambda, ast.FunctionDef, ast.AsyncFunctionDef)):
            if n is not fn:
                bound.add(getattr(n, "name", ""))
            bound |= {a.arg for a in n.args.args}
        elif isinstance(n, ast.NamedExpr) and isinstance(n.target, ast.Name):
            bound.add(n.target.id)
        elif isinstance(n, (ast.Import, ast.ImportFrom)):
            bound |= {a.asname or a.name.split(".")[0] for a in n.names}
        elif isinstance(n, ast.Global):
            bound |= set(n.names)
    return bound


@pytest.mark.parametrize("path", APP_FILES, ids=lambda p: p.name)
def test_no_unbound_names_in_any_function(path):
    tree = ast.parse(path.read_text())
    module_names = _module_globals(tree)
    problems = []
    for fn in [n for n in ast.walk(tree)
               if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]:
        bound = _bound_in(fn)
        for n in ast.walk(fn):
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load):
                if n.id not in bound and n.id not in module_names:
                    problems.append(f"{path.name}:{n.lineno} {fn.name}() -> {n.id!r}")
    assert not problems, (
        "name(s) used but never bound — a Streamlit panel raises these only "
        "when its tab is clicked:\n  " + "\n  ".join(sorted(set(problems))))


def test_the_check_would_have_caught_the_real_bug(tmp_path):
    """A guard that cannot fail is not a guard."""
    bad = tmp_path / "bad_panel.py"
    bad.write_text(
        "def panel():\n"
        "    for key, cfg in {}.items():\n"
        "        print(f'pose_{a}')\n")
    tree = ast.parse(bad.read_text())
    module_names = _module_globals(tree)
    fn = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)][0]
    unbound = [n.id for n in ast.walk(fn)
               if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)
               and n.id not in _bound_in(fn) and n.id not in module_names]
    assert "a" in unbound
