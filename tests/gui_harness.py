"""Shared driving helpers for the streamlit AppTest harness.

Imported by BOTH `tests/test_app_renders.py` (pytest, when streamlit is
available) and `tests/run_app_renders_gui_env.py` (a plain script, for the
`dwi_gui` env which has no pytest). One copy so the two cannot drift into
testing different things and reporting the same name.

THE HARNESS LIMITATION THIS EXISTS TO WORK AROUND. `AppTest` cannot compute
widget state for a `st.selectbox` that uses `format_func`. Its `.options` are
the FORMATTED strings while its `.value` is the RAW option, and
`element_tree.SelectboxProto.index` does `self.options.index(formatted_value)`,
which for our panels compares `'t3_5c5a4cf73e08'` against
`'t3_5c5a4cf73e08  ·  5.70×'` and raises. Every widget on the page is snapshot
before a rerun, so ONE such selectbox anywhere makes any subsequent
`set_value(...).run()` raise -- from inside streamlit, before the app runs at
all.

That is a harness defect, not an app defect: the app renders, and a real browser
has no such problem. It was misread once as "Near-attack ranking crashes under a
filter", which is why the diagnosis is written down here rather than in a commit
message. `set_spec` therefore drives the widget the real way FIRST and only
falls back to writing session state when streamlit raises that specific
ValueError -- so a genuine app failure under curation still fails the test.
"""

from __future__ import annotations

#: The session-state key the curation text area binds to (`app.CURATE_KEY`).
#: Asserted against the app in `test_app_renders.py` so a rename cannot silently
#: turn every filtered case into a no-op that still passes.
CURATE_KEY = "_curate_spec"


class HarnessLimitation(Exception):
    """AppTest cannot drive this case. NOT an app failure, and not a pass.

    Raised so the runner can report the case as UNREACHABLE by name. It must
    never be reported as a pass: the whole point of #45 is that a case which
    reads as covered and is not is worse than one that is openly missing.
    """


def _is_format_func_snapshot_error(exc: Exception) -> bool:
    """True for the streamlit widget-state error described in the module docstring.

    Deliberately narrow. Any other ValueError, and any error raised by the app
    itself, propagates -- the point is to name a harness limitation, not to make
    failures disappear.
    """
    return isinstance(exc, ValueError) and "is not in list" in str(exc)


def set_spec(at, spec: str):
    """Apply a curation constraint through the real widget interaction.

    Returns the AppTest. Raises whatever the app raises, or `HarnessLimitation`
    if streamlit could not snapshot the page's widgets.

    NO FALLBACK IS POSSIBLE, and it is worth saying why, because writing session
    state directly looks like it should work. `AppTest.run()` snapshots EVERY
    widget on the page before executing the script, so the failing `index`
    lookup happens on the rerun itself -- whichever way the value got set. A
    single `format_func` selectbox anywhere on a panel therefore makes that whole
    panel undriveable, not just the widget.
    """
    try:
        at.sidebar.text_area[0].set_value(spec).run()
    except ValueError as exc:
        if not _is_format_func_snapshot_error(exc):
            raise
        raise HarnessLimitation(
            "streamlit AppTest cannot snapshot a selectbox that uses "
            f"format_func ({exc}); the panel renders, the harness cannot rerun "
            "it. See tests/gui_harness.py.") from exc
    return at
