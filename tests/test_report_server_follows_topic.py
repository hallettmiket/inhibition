"""The report server must serve the CURRENT topic, not the one it started under.

THE DEFECT. `serve_reports.py` resolved its document root once, in `main()`,
and handed it to the handler. Its own docstring says the root comes from
`run_paths` "so the server cannot end up serving a superseded topic's pages --
which it did once, for hours, from a literal path." That fixed the literal and
left the LIFETIME: resolving at startup is a pin on the topic as it was at
startup, and these processes run for weeks.

Measured 2026-08-31: the server on :8931 had been up 14 days and was returning
`mdprio_reports_nac_v5/index.html` byte-for-byte while `run.topic` was `nac_v6`
and `mdprio_reports_nac_v6/` was being rebuilt beside it every few minutes.
Nothing looked wrong -- the pages render and the numbers are a real screen's,
just the previous screen's. Catalogue #25, in the component that claims immunity.

WHY A UNIT TEST AND NOT AN INTEGRATION ONE. Standing up a server, bumping
`run.topic` and re-fetching would exercise it end to end, but `run.topic` is
global state that detached supervisors poll (CLAUDE.md), so a test must never
write it. `translate_path` is the seam the whole behaviour hangs off, so it is
tested directly with the topic faked at the `run_paths` boundary.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def _load():
    spec = importlib.util.spec_from_file_location(
        "serve_reports", REPO / "scripts" / "serve_reports.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def m():
    return _load()


class _FakeHandler:
    """The bits of SimpleHTTPRequestHandler `translate_path` touches."""

    def __init__(self, directory):
        self.directory = directory

    def translate_path(self, path):        # stands in for the base class
        return f"{self.directory}/{path.lstrip('/')}"


def test_live_handler_reresolves_the_root_on_every_request(m, monkeypatch, tmp_path):
    """The root must follow `reports_dir()`, not the value it was built with."""
    a, b = tmp_path / "reports_topic_a", tmp_path / "reports_topic_b"
    a.mkdir(); b.mkdir()

    class H(m.LiveRun, _FakeHandler):
        def __init__(self, directory):
            _FakeHandler.__init__(self, directory)

    h = H(str(a))
    monkeypatch.setattr(m.rp, "reports_dir", lambda *_a, **_k: a)
    assert h.translate_path("/index.html").startswith(str(a))

    # The topic moves on. The SAME handler instance must follow it.
    monkeypatch.setattr(m.rp, "reports_dir", lambda *_a, **_k: b)
    got = h.translate_path("/index.html")
    assert got.startswith(str(b)), (
        f"the handler is still serving {a}; the root was pinned at construction, "
        "which is the 14-day nac_v5 defect")


def test_the_guard_can_fail(m, monkeypatch, tmp_path):
    """A handler that pins its root must FAIL this test — prove it does.

    Written because this repo has shipped two tests that passed for free.
    `NoCache` is the pinned handler (correct for `--archive`, a frozen
    snapshot), so it is the honest negative control.
    """
    a, b = tmp_path / "ra", tmp_path / "rb"
    a.mkdir(); b.mkdir()

    class Pinned(m.NoCache, _FakeHandler):
        def __init__(self, directory):
            _FakeHandler.__init__(self, directory)

    h = Pinned(str(a))
    monkeypatch.setattr(m.rp, "reports_dir", lambda *_a, **_k: b)
    assert h.translate_path("/index.html").startswith(str(a)), \
        "the negative control no longer pins, so the positive test proves nothing"


def test_archive_mode_still_pins(m):
    """An archive is a FROZEN snapshot and must NOT follow the live topic.

    The two modes want opposite things, and that is deliberate: a released run
    has to stay browsable after the topic moves on.
    """
    import inspect
    src = inspect.getsource(m.main)
    arch = src.split("if args.archive:")[1].split("else:")[0]
    assert "NoCache" in arch and "LiveRun" not in arch, \
        "archive mode must keep the pinned handler"
    live = src.split("else:")[1]
    assert "LiveRun" in live, "live mode must use the re-resolving handler"


def test_live_root_is_never_a_typed_path():
    """No literal report directory may reappear in EXECUTABLE code.

    Prose must stay free to name the directories -- the docstrings above exist
    to explain which topic was wrongly served, and cannot do that without
    writing it down. An earlier version of this test stripped only `#` comments
    and flagged its own module docstring, which is the same mistake the
    version-pin test made when it reported four docstrings as offenders.
    """
    import ast
    src = (REPO / "scripts" / "serve_reports.py").read_text()
    tree = ast.parse(src)
    # Every string CONSTANT that is not a docstring, plus every attribute path.
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            d = ast.get_docstring(node, clean=False)
            if d is not None:
                docstrings.add(d)
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if node.value in docstrings:
                continue
            for lit in ("mdprio_reports_nac_v", "/data/lab_vm"):
                if lit in node.value:
                    bad.append(f"line {node.lineno}: {node.value[:60]!r}")
    assert not bad, "a typed report root is back in executable code:\n  " + \
                    "\n  ".join(bad)


# ---------------------------------------------------------------------------
# One server, every run (@tt8804: "combine them all and I just select which run
# I want to see"). Four servers on four ports, with no way to tell from a page
# which run you were looking at, is what produced D0103.
# ---------------------------------------------------------------------------


def test_run_roots_finds_live_superseded_and_archived(m, monkeypatch, tmp_path):
    """All three kinds must appear; they are not equivalent and are labelled."""
    bl = tmp_path / "blacksmith"
    (bl / "mdprio_reports_run_a").mkdir(parents=True)
    (bl / "mdprio_reports_run_b").mkdir(parents=True)
    (bl / "gui_archive" / "frozen_1").mkdir(parents=True)
    monkeypatch.setattr(m.rp, "BLACKSMITH", bl)
    got = m.run_roots()
    assert set(got) == {"run_a", "run_b", "frozen_1"}
    assert got["run_a"].name == "mdprio_reports_run_a"
    assert got["frozen_1"].parent.name == "gui_archive"


def test_the_first_path_segment_selects_the_run(m, monkeypatch, tmp_path):
    bl = tmp_path / "blacksmith"
    for r in ("mdprio_reports_run_a", "mdprio_reports_run_b"):
        (bl / r).mkdir(parents=True)
    monkeypatch.setattr(m.rp, "BLACKSMITH", bl)

    class H(m.MultiRun, _FakeHandler):
        def __init__(self, directory):
            _FakeHandler.__init__(self, directory)

    h = H(str(bl))
    assert h.translate_path("/run_a/modes.html").startswith(
        str(bl / "mdprio_reports_run_a"))
    assert h.translate_path("/run_b/modes.html").startswith(
        str(bl / "mdprio_reports_run_b"))


def test_an_unknown_run_does_not_fall_through_to_another_one(m, monkeypatch, tmp_path):
    """Serving run A's page at run B's URL is the whole failure mode."""
    bl = tmp_path / "blacksmith"
    (bl / "mdprio_reports_run_a").mkdir(parents=True)
    monkeypatch.setattr(m.rp, "BLACKSMITH", bl)

    class H(m.MultiRun, _FakeHandler):
        def __init__(self, directory):
            _FakeHandler.__init__(self, directory)

    h = H(str(bl))
    got = h.translate_path("/not_a_run/index.html")
    assert str(bl / "mdprio_reports_run_a") not in got, \
        "an unknown run resolved into a real run's tree"


def test_the_chooser_labels_the_current_run_and_says_which_are_stale(m, monkeypatch, tmp_path):
    """A superseded page renders like a current one — the label is the guard."""
    bl = tmp_path / "blacksmith"
    for r in ("mdprio_reports_cur", "mdprio_reports_old"):
        d = bl / r
        d.mkdir(parents=True)
        (d / "index.html").write_text("<html></html>")
    (bl / "gui_archive" / "frozen_1").mkdir(parents=True)
    (bl / "gui_archive" / "frozen_1" / "index.html").write_text("<html></html>")
    monkeypatch.setattr(m.rp, "BLACKSMITH", bl)
    monkeypatch.setattr(m.rp, "topic", lambda *_a, **_k: "cur")
    page = m._chooser_html().decode()
    assert "CURRENT RUN" in page
    assert "SUPERSEDED" in page
    assert "ARCHIVED" in page
    # and the current one must not ALSO be marked superseded
    cur_row = [l for l in page.split("<tr") if ">cur<" in l]
    assert cur_row and "SUPERSEDED" not in cur_row[0]


def test_the_chooser_is_generated_not_stored(m):
    """A menu built once at startup is D0103 with a nicer front end."""
    import inspect
    src = inspect.getsource(m.MultiRun.do_GET)
    assert "_chooser_html()" in src, \
        "the chooser must be produced per request, not read from a file"
