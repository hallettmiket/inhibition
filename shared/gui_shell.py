"""The step navigation every page in the GUI carries (#63).

@tt8804: "integrate the ranking and results pages into one uniform gui with
arrows showing each page as a distinct step ... Home page showing input target
and parameters, ranking results, sweep results, MD results".

WHY A SHARED MODULE AND NOT A COPY PER PAGE. The pages were built at different
times by different scripts -- `modes.html` from `shared/mode_ranking`,
`combined.html` from `scripts/mdprio_combine`, `pipeline.html` from
`shared/pipeline_schematic` -- and each carries its own `<style>`. A stepper
pasted into three templates is three steppers that drift, and the first thing a
reader would notice is the one page whose "you are here" is wrong.

THE STEPS ARE THE PIPELINE'S REAL ORDER, which is not the order the pages were
written in. A molecule is docked and scored (RANK), a subset of its modes earn a
10 ns triage (SWEEP), a subset of those earn 100 ns (MD). Home states what the
run was configured to do before any of it. Showing them as arrows makes the
funnel legible: each step's count is smaller than the last, and a reader who
cannot say why has found something worth asking about.
"""

from __future__ import annotations

#: (file, label, one-line description). Order IS the pipeline order.
STEPS = [
    ("index.html", "Home", "target, receptor and the rules this run was given"),
    ("modes.html", "Ranking", "every molecule and every binding mode, scored"),
    ("sweep.html", "Sweep", "10 ns triage: which modes reach attack geometry"),
    # "MD results", not "MD": the other three steps name an OUTPUT (a ranking, a
    # sweep), and a step called "MD" names a method instead, which reads as a
    # setting rather than somewhere to go. @tt8804.
    ("combined.html", "MD results",
     "100 ns runs, interactions, and the shortlist"),
]

CSS = """
/* --- step navigation (#63), identical on every page ------------------- */
#steps{display:flex;align-items:stretch;gap:0;padding:0 14px;background:var(--raise);
 border-bottom:1px solid var(--rule);overflow-x:auto;flex:0 0 auto}
#steps a{display:flex;flex-direction:column;justify-content:center;gap:1px;
 padding:7px 20px 7px 26px;text-decoration:none;color:var(--muted);position:relative;
 white-space:nowrap;min-width:0}
#steps a:first-child{padding-left:12px}
/* The arrow is a CSS chevron rather than a character, so it scales with the row
   and cannot be selected or read out as punctuation. */
#steps a:not(:first-child)::before{content:"";position:absolute;left:6px;top:50%;
 width:7px;height:7px;border-top:2px solid var(--rule);border-right:2px solid var(--rule);
 transform:translateY(-50%) rotate(45deg)}
#steps a .sl{font:600 12.5px var(--sans);letter-spacing:.01em}
#steps a .sd{font:400 10px var(--sans);opacity:.85;max-width:30ch;overflow:hidden;
 text-overflow:ellipsis}
#steps a .sn{font-family:var(--mono);font-size:10px;opacity:.8}
#steps a:hover{background:var(--blue-pale);color:var(--ink)}
#steps a.on{color:var(--navy);background:var(--paper);
 box-shadow:inset 0 -2px 0 var(--blue)}
#steps a.on .sl{font-weight:700}
@media(max-width:760px){#steps a .sd{display:none}}
"""


def nav(current: str, counts: dict | None = None) -> str:
    """The stepper, with `current` highlighted.

    `counts` maps a step's file to a short count string shown under its label --
    the funnel in numbers. A step with no count simply omits it rather than
    printing a zero, because "not measured yet" and "measured, none" are
    different and a bare 0 says the second.
    """
    counts = counts or {}
    out = ['<nav id="steps">']
    for href, label, desc in STEPS:
        cls = " class=\"on\"" if href == current else ""
        n = counts.get(href)
        out.append(
            f'<a href="{href}"{cls} title="{desc}">'
            f'<span class="sl">{label}</span>'
            + (f'<span class="sn">{n}</span>' if n else
               f'<span class="sd">{desc}</span>')
            + '</a>')
    out.append('</nav>')
    return "".join(out)


def step_counts() -> dict:
    """The counts under each step label — ONE source for every page.

    WHY THIS EXISTS. Four builders computed these independently: `build_gui`
    from `sweep_state.summary()`, `sweep_combine` from its own inline dict,
    `mdprio_combine` and the ranking page from `mode_ranking._step_counts`,
    which read `mdprio_reports/sweep_state.json` -- the UNSCOPED directory. So
    the Ranking page's nav said "447 ok" (3.0.0's sweep) while the Sweep page
    beside it said 34 (this run's), and both were reading a file that existed.
    @tt8804: "on ranking it shows sweep is 447 okay while clicking on sweep
    shows 32 ok, can we make this whole gui cohesive".

    Four sources for one row of numbers is four numbers. This is the one, and it
    reads the pipeline's own probes -- which count artefacts on disk under
    `run.topic`, are already tested, and are what the dashboard shows. A page
    cannot now disagree with the dashboard or with another page.

    Absent counts are OMITTED, never shown as 0: "not measured yet" and
    "measured, none" are different claims and a bare 0 makes the second. A count
    the pipeline reports as `unknown` is omitted for the same reason.
    """
    try:
        from . import pipeline as pl
        st = {s["name"]: s for s in pl.status()["stages"]}
    except Exception:                                      # noqa: BLE001
        return {}

    out: dict[str, str] = {}

    def n(stage):
        s = st.get(stage) or {}
        return s.get("done"), s.get("total"), s.get("state")

    # Ranking: the number of MODES ranked, which is the rail's own length --
    # not the (1/1) "did the table get written" the pipeline tracks.
    try:
        from . import mode_ranking as mr
        r = mr.gather()
        if r is not None and len(r):
            out["modes.html"] = f"{len(r):,} modes"
    except Exception:                                      # noqa: BLE001
        pass

    done, total, state = n("sweep")
    if state != "unknown" and done is not None and total:
        pend = max(0, total - done)
        out["sweep.html"] = f"{done} ok · {pend} pending"

    done, total, state = n("production")
    if state != "unknown" and done is not None and total:
        out["combined.html"] = f"{done} of {total} at 100 ns"

    return out
