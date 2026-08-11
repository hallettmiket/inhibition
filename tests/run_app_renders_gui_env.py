"""
Purpose: run the GUI render tests in the env that actually has streamlit.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-08
Input: integration/app/app.py, driven through streamlit's own AppTest harness
Output: pass/fail per panel

WHY THIS EXISTS. `tests/test_app_renders.py` is the strongest test in the repo --
it executes every panel headlessly and catches the class of bug that only appears
when someone clicks a tab. It has almost certainly **never run**: it needs
streamlit, streamlit lives in `dwi_gui`, and `dwi_gui` has no pytest. So it
`importorskip`s under the suite's own env and reports as a skip, which reads as
"covered" in a summary line.

That gap is exactly where tonight's changes landed. The 2.2.0 panel was rewritten
twice, the viewer three times, and a new entry was added to `curate.PANEL_SCOPE`
-- a missing scope declaration had already crashed the app once this session, and
that is precisely what this harness catches.

A PLAIN SCRIPT, NOT A PYTEST FILE, so it runs under `dwi_gui/bin/python` with no
new dependency. The assertions mirror `test_app_renders.py`; when pytest exists
in the GUI env this can be deleted in favour of it.

    /data/lab_vm/envs/dwi_gui/bin/python tests/run_app_renders_gui_env.py
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tests"))
import gui_harness  # noqa: E402
APP = REPO / "integration" / "app" / "app.py"
sys.path.insert(0, str(REPO / "integration" / "app"))
sys.path.insert(0, str(REPO))

SPEC = "no chlorine\nmw < 450"
BAD_SPEC = "no chlorien"


def _messages(at) -> str:
    out: list[str] = []
    for kind in ("info", "success", "warning", "error", "markdown", "caption"):
        try:
            out += [str(e.value) for e in getattr(at, kind)]
        except Exception:                              # noqa: BLE001
            pass
    return "\n".join(out)


def _run(panel: str, spec: str | None = None):
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(str(APP), default_timeout=300)
    at.run()
    if at.exception:
        raise AssertionError(f"app failed on first run: {at.exception}")
    at.sidebar.radio[0].set_value(panel).run()
    if at.exception:
        raise AssertionError(f"{panel} raised: {at.exception}")
    if spec is not None:
        gui_harness.set_spec(at, spec)
        if at.exception:
            raise AssertionError(f"{panel} raised under curation: {at.exception}")
    return at


def main() -> int:
    import curate
    scopes = list(curate.PANEL_SCOPE)
    print(f"{len(scopes)} panels declared in curate.PANEL_SCOPE\n")

    fails = 0
    unreachable: list[tuple[str, str]] = []
    for s in scopes:
        p = s.panel
        # 1. renders at all
        try:
            _run(p)
            print(f"  PASS  render            {p}")
        except Exception as exc:                       # noqa: BLE001
            fails += 1
            print(f"  FAIL  render            {p}\n        {exc}")
            continue
        # 2. renders under a real curation filter
        try:
            at = _run(p, SPEC)
            print(f"  PASS  under filter      {p}")
        except gui_harness.HarnessLimitation as exc:
            # NOT a pass. Counted and named, so the summary cannot read as
            # covered -- which is the whole complaint in #45.
            unreachable.append((p, str(exc).split(";")[0]))
            print(f"  UNREACHABLE  under filter  {p}")
            continue
        except Exception as exc:                       # noqa: BLE001
            fails += 1
            print(f"  FAIL  under filter      {p}\n        {exc}")
            continue
        # 3. says which side of the curation line it is on. Silence is how the
        #    original bug read to the user: a panel that quietly ignored the
        #    filter looked identical to one that honoured it.
        msg = _messages(at)
        want = "Curation active" if s.filtered else "not applied here"
        if want in msg:
            print(f"  PASS  declares scope    {p}  ({'filtered' if s.filtered else 'unfiltered'})")
        else:
            fails += 1
            print(f"  FAIL  declares scope    {p}\n        expected {want!r} in the rendered output")
        # 4. a mis-typed constraint must REFUSE, never silently filter nothing
        try:
            at2 = _run(p, BAD_SPEC)
            if "not understood" in _messages(at2):
                print(f"  PASS  refuses a typo    {p}")
            else:
                fails += 1
                print(f"  FAIL  refuses a typo    {p}\n        a mis-parsed constraint filtered nothing, silently")
        except Exception as exc:                       # noqa: BLE001
            fails += 1
            print(f"  FAIL  refuses a typo    {p}\n        {exc}")

    if unreachable:
        print("\n  NOT COVERED -- the harness cannot drive these, the app is fine:")
        for panel, why in unreachable:
            print(f"    {panel}: {why}")
    print(f"\n{'all panels render' if not fails else f'{fails} FAILURES'}"
          f"{f'; {len(unreachable)} case(s) not coverable by AppTest' if unreachable else ''}")
    return 1 if fails else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:                                  # noqa: BLE001
        traceback.print_exc()
        sys.exit(2)
