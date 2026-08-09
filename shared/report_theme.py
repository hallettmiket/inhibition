"""
Purpose: the house style for every report this project emits — one theme, defined once.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: none
Output: CSS and small HTML helpers importable by any report script

WHY A SHARED MODULE. Three report generators existed with three hand-written
stylesheets, so "the elevation report" and "the overnight report" looked like
documents from different organisations and a number formatted one way in one and
another way in the next. A house style is not decoration here: these reports go
to chemists who have to compare a candidate across two of them.

WHITE GROUND, BLUE ACCENT, at @tt8804's direction — the register of a pharma
technical report rather than a lab notebook. Deliberately NOT theme-switching:
a report that renders differently for two readers is a report two people cannot
discuss, and these get printed and pasted into decks.

THE RULES THE STYLE ENCODES, which are project rules rather than taste:

  - Numbers are tabular-lining everywhere, so a column of them can be scanned.
  - A "not measured" is an em-dash, never a zero. Rendering an unmeasured
    quantity as 0.000 is how "we did not run it" becomes "we ran it and it
    failed", and that has cost this project real time.
  - Verdicts carry a colour AND a word. Colour alone fails for the ~8% of men
    with red-green deficiency, and these tables are pass/fail.
  - Caveats get a box, not a footnote. Every serious error in 2.0.0 survived
    because its qualification was somewhere the reader was not looking.
"""

from __future__ import annotations

#: Deep navy for structure, mid blue for accents, a cool grey ramp for the rest.
#: Semantic colours are separate from the accent so "good" never reads as
#: "branded".
PALETTE = {
    "ink": "#10233f",
    "navy": "#003087",
    "blue": "#0072ce",
    "blue_pale": "#e8f1fb",
    "rule": "#d6dee8",
    "muted": "#5b6b80",
    "paper": "#ffffff",
    "raise": "#f5f8fc",
    "good": "#0f7a54",
    "warn": "#8a5a00",
    "bad": "#b3261e",
}

CSS = """
:root{
  --ink:#10233f; --navy:#003087; --blue:#0072ce; --blue-pale:#e8f1fb;
  --rule:#ccd6e2; --muted:#5b6b80; --paper:#eef2f7; --raise:#e3eaf2;
  --card:#ffffff;
  --good:#0f7a54; --warn:#8a5a00; --bad:#b3261e;
  --sans:"Helvetica Neue",Helvetica,Arial,system-ui,sans-serif;
  --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
}
:root[data-theme="dark"]{
  --ink:#dfe7f0; --navy:#8ab4e8; --blue:#6aa9e0; --blue-pale:#16283a;
  --rule:#25333f; --muted:#93a3b4; --paper:#0e151c; --raise:#16202a;
  --card:#131c25; --good:#4fc4a0; --warn:#e0b66a; --bad:#e08a70;
}
:root[data-theme="dark"] img[src^="data:image/png"]{filter:invert(1) hue-rotate(180deg)}
:root[data-theme="dark"] .structbox{background:#fff;border-radius:4px;padding:.3rem}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
  font-size:15.5px;line-height:1.55;-webkit-font-smoothing:antialiased}
.wrap{max-width:900px;margin:0 auto;padding:0 28px 70px}
/* NO .wrap ELEMENT IS EMITTED -- the report writes its blocks straight into
   <body>, so the measure has to be applied to those blocks or it does nothing.
   That is why the content sat hard against the frame edge in the catalogue
   viewer instead of centring in it. */
*,*::before,*::after{box-sizing:border-box}
body>*{max-width:940px;margin-left:auto;margin-right:auto;
  padding-left:26px;padding-right:26px}
/* The mast, the structure card and each panel are all direct children, so they
   must resolve to the SAME left edge. Without border-box the bordered blocks
   grew by their own padding and sat further left than the masthead text -- which
   is what read as "the rest of the viewer is to the left". */
body>header.mast{padding-top:26px}
body>.structrow,body>details.panel{margin-left:auto;margin-right:auto}
/* PANELS. In the catalogue viewer these reports are read inside a frame beside
   a selector, so everything is collapsible and the heavy things (movie, plots)
   start closed -- the reader opens what they want instead of scrolling past it. */
details.panel{border:1px solid var(--rule);border-radius:5px;margin:0 0 .9rem;
  background:var(--card);overflow:hidden}
details.panel>summary{cursor:pointer;list-style:none;padding:.6rem .9rem;
  font-weight:600;color:var(--navy);display:flex;align-items:center;gap:.5rem;
  background:var(--raise);border-radius:3px}
details.panel>summary::-webkit-details-marker{display:none}
details.panel>summary::before{content:"";width:0;height:0;flex:none;
  border-left:6px solid var(--blue);border-top:4.5px solid transparent;
  border-bottom:4.5px solid transparent;transition:transform .12s ease}
details.panel[open]>summary::before{transform:rotate(90deg)}
details.panel>summary:focus-visible{outline:2px solid var(--blue);outline-offset:2px}
details.panel>summary .hint{font-weight:400;color:var(--muted);font-size:.8rem;
  margin-left:auto}
details.panel .pbody{padding:.9rem;background:var(--card)}
/* the molecule itself, directly under the masthead -- what the reader wants to
   see before any number */
.structrow{display:flex;gap:1.1rem;align-items:center;background:var(--card);
  border:1px solid var(--rule);border-radius:5px;padding:.8rem 1rem;margin:0 0 1rem}
.structbox{flex:none;width:300px;max-width:52%}
.structbox svg{width:100%;height:auto;display:block}
.structnote{color:var(--muted);font-size:.9rem}
.structnote b{color:var(--navy);text-transform:none}
@media(max-width:640px){.structrow{flex-direction:column;align-items:flex-start}
  .structbox{width:100%;max-width:none}}
details.panel .pbody>img{max-width:100%;height:auto;display:block;margin:0 auto}
@media(prefers-reduced-motion:reduce){details.panel>summary::before{transition:none}}

/* masthead ---------------------------------------------------------------- */
.mast{border-top:5px solid var(--navy);padding:24px 0 18px;margin-bottom:22px;
  border-bottom:1px solid var(--rule)}
.eyebrow{font-family:var(--mono);font-size:.68rem;letter-spacing:.16em;
  text-transform:uppercase;color:var(--blue);font-weight:600;margin-bottom:.9rem}
h1{font-size:clamp(1.7rem,3.4vw,2.4rem);margin:0 0 .5rem;font-weight:600;
  letter-spacing:-.015em;color:var(--navy);text-wrap:balance}
.standfirst{font-size:1.05rem;color:var(--muted);max-width:62ch;margin:0 0 1.4rem}
.facts{display:flex;flex-wrap:wrap;gap:0 2.2rem;border-top:1px solid var(--rule);
  padding-top:.2rem}
.facts div{padding-top:.75rem;font-family:var(--mono);font-size:.72rem;
  letter-spacing:.06em;text-transform:uppercase;color:var(--muted)}
.facts b{display:block;font-family:var(--sans);font-size:1.15rem;font-weight:600;
  letter-spacing:-.01em;color:var(--ink);text-transform:none}

/* sections ---------------------------------------------------------------- */
section{margin:0 0 2.8rem}
.shead{display:flex;gap:.9rem;align-items:baseline;border-bottom:2px solid var(--navy);
  padding-bottom:.5rem;margin-bottom:1.1rem}
.snum{font-family:var(--mono);font-size:.72rem;color:var(--blue);font-weight:600;
  letter-spacing:.08em;white-space:nowrap;padding-top:.3rem}
h2{font-size:1.35rem;margin:0;font-weight:600;color:var(--navy);letter-spacing:-.01em}
h3{font-size:1rem;margin:1.6rem 0 .4rem;font-weight:600;color:var(--ink)}
.sub{color:var(--muted);font-size:.9rem;margin:.15rem 0 0}
p{margin:0 0 .85rem;max-width:78ch}
a{color:var(--blue)}
code{font-family:var(--mono);font-size:.87em;background:var(--raise);
  padding:.08em .34em;border-radius:2px;border:1px solid var(--rule)}

/* tables ------------------------------------------------------------------ */
.scroll{overflow-x:auto;margin:1rem 0;border:1px solid var(--rule)}
table{border-collapse:collapse;width:100%;font-size:.85rem}
caption{caption-side:top;text-align:left;font-size:.78rem;color:var(--muted);
  padding:.5rem .8rem;background:var(--raise);border-bottom:1px solid var(--rule)}
th,td{padding:.5rem .75rem;text-align:left;border-bottom:1px solid var(--rule);
  white-space:nowrap}
thead th{background:var(--navy);color:#fff;font-family:var(--mono);font-size:.66rem;
  letter-spacing:.07em;text-transform:uppercase;font-weight:600;border-bottom:none}
tbody tr:nth-child(even){background:var(--raise)}
tbody tr:last-child td{border-bottom:none}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums lining-nums;
  font-family:var(--mono)}
tr.ref td{background:var(--blue-pale)}
tr.ref td:first-child{box-shadow:inset 3px 0 0 var(--blue)}
.na{color:var(--muted)}

/* verdicts: colour AND word, never colour alone ---------------------------- */
.pill{display:inline-block;font-family:var(--mono);font-size:.66rem;font-weight:600;
  letter-spacing:.06em;padding:.14rem .5rem;border-radius:2px;border:1px solid currentColor;
  text-transform:uppercase;white-space:nowrap}
.pill.pass{color:var(--good)} .pill.fail{color:var(--bad)}
.pill.pend{color:var(--muted)} .pill.warn{color:var(--warn)}

/* callouts ---------------------------------------------------------------- */
.callout{border:1px solid var(--rule);border-left:4px solid var(--blue);
  background:var(--raise);padding:.9rem 1.1rem;margin:1.2rem 0}
.callout.warn{border-left-color:var(--warn)}
.callout.bad{border-left-color:var(--bad)}
.callout.good{border-left-color:var(--good)}
.callout p:last-child{margin-bottom:0}
.ctitle{font-family:var(--mono);font-size:.66rem;letter-spacing:.1em;
  text-transform:uppercase;color:var(--muted);font-weight:600;margin-bottom:.35rem}

/* figures ----------------------------------------------------------------- */
figure{margin:1.2rem 0}
figure img{width:100%;height:auto;display:block;border:1px solid var(--rule)}
figcaption{font-size:.79rem;color:var(--muted);margin-top:.5rem;max-width:80ch}

/* molecule cards ---------------------------------------------------------- */
.card{border:1px solid var(--rule);border-top:3px solid var(--blue);
  background:var(--paper);padding:1rem 1.2rem;margin:1.1rem 0}
.card h3{margin-top:0;display:flex;align-items:center;gap:.6rem;flex-wrap:wrap}
.row{display:flex;gap:1.3rem;flex-wrap:wrap;align-items:flex-start}
.row img{max-width:100%;border:1px solid var(--rule);background:#fff}

/* trajectory viewer ------------------------------------------------------- */
.viewer{border:1px solid var(--rule);background:var(--raise);margin:1rem 0}
/* 3Dmol absolutely-positions its canvas inside the element it is given, so that
   element MUST be positioned and sized or the canvas escapes onto the page. */
.glbox{position:relative;width:100%;height:430px;overflow:hidden;background:#fff}
.glbox canvas{position:absolute;top:0;left:0}
.vctl{display:flex;align-items:center;gap:.85rem;padding:.6rem .9rem;
  border-top:1px solid var(--rule);flex-wrap:wrap;background:var(--paper)}
button.play{font-family:var(--mono);font-size:.75rem;font-weight:600;
  background:var(--navy);color:#fff;border:none;padding:.4rem .9rem;
  border-radius:2px;cursor:pointer;min-width:76px}
button.play:hover{background:var(--blue)}
button.play:focus-visible,input:focus-visible{outline:2px solid var(--blue);
  outline-offset:2px}
.vctl input[type=range]{flex:1;min-width:140px;accent-color:var(--blue)}
.readout{font-family:var(--mono);font-size:.74rem;color:var(--muted);white-space:nowrap}
.readout b{color:var(--ink);font-variant-numeric:tabular-nums}
.legend{display:flex;gap:1.1rem;flex-wrap:wrap;font-size:.72rem;color:var(--muted);
  padding:.5rem .9rem;border-top:1px solid var(--rule);font-family:var(--mono);
  background:var(--paper)}
.sw{display:inline-block;width:9px;height:9px;margin-right:.35rem;
  vertical-align:middle;border:1px solid rgba(0,0,0,.15)}

ul,ol{padding-left:1.15rem;max-width:78ch}
li{margin-bottom:.4rem}
.foot{border-top:2px solid var(--navy);padding-top:1rem;margin-top:2.4rem;
  font-size:.76rem;color:var(--muted);font-family:var(--mono);line-height:1.7}

@media print{
  body{font-size:10.5pt} .wrap{max-width:none;padding:0}
  .viewer,.vctl,.legend{display:none}   /* an animation does not print */
  section{break-inside:avoid} .card{break-inside:avoid}
  thead th{background:#003087 !important;-webkit-print-color-adjust:exact;
    print-color-adjust:exact}
}
@media (max-width:680px){.wrap{padding:0 16px 60px}.glbox{height:320px}}
"""

#: matplotlib rcParams matching the CSS, so a figure looks like it belongs to the
#: page rather than merely sitting on it.
MPL = {
    "figure.facecolor": "#ffffff", "axes.facecolor": "#ffffff",
    "savefig.facecolor": "#ffffff",
    "text.color": "#10233f", "axes.labelcolor": "#10233f",
    "xtick.color": "#5b6b80", "ytick.color": "#5b6b80",
    "axes.edgecolor": "#d6dee8", "grid.color": "#d6dee8",
    "font.family": "DejaVu Sans", "font.size": 9.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 130, "savefig.bbox": "tight",
}

#: Series colours. `accent` is the data, `ref` marks a reference or anchor,
#: `alert` marks the failure region. Semantic, not decorative.
SERIES = {"accent": "#0072ce", "ref": "#0f7a54", "alert": "#b3261e",
          "muted": "#5b6b80", "grid": "#d6dee8"}


def num(v, fmt: str = "{:.3f}", dash: str = "—") -> str:
    """Format a number, or an em-dash when it is absent.

    NEVER returns 0 for a missing value. Rendering "not measured" as 0.000 is how
    a run that did not happen becomes a run that failed, and this project has
    paid for that more than once.
    """
    if v is None:
        return f'<span class="na">{dash}</span>'
    try:
        if v != v:                                   # NaN
            return f'<span class="na">{dash}</span>'
        return fmt.format(v)
    except (TypeError, ValueError):
        return f'<span class="na">{dash}</span>'


def pill(verdict: str) -> str:
    """A verdict as colour AND word. Colour alone is not an accessible signal."""
    cls = {"pass": "pass", "fail": "fail", "warn": "warn"}.get(
        str(verdict).lower(), "pend")
    return f'<span class="pill {cls}">{verdict}</span>'


def callout(title: str, body: str, kind: str = "") -> str:
    k = f" {kind}" if kind else ""
    return (f'<div class="callout{k}"><div class="ctitle">{title}</div>{body}</div>')


def masthead(title: str, standfirst: str, eyebrow: str, facts: list[tuple]) -> str:
    f = "".join(f"<div><b>{v}</b>{k}</div>" for v, k in facts)
    return (f'<header class="mast"><div class="eyebrow">{eyebrow}</div>'
            f'<h1>{title}</h1><p class="standfirst">{standfirst}</p>'
            f'<div class="facts">{f}</div></header>')


def section(num_: str, title: str, sub: str = "") -> str:
    s = f'<p class="sub">{sub}</p>' if sub else ""
    return (f'<div class="shead"><div class="snum">{num_}</div>'
            f'<div><h2>{title}</h2>{s}</div></div>')
