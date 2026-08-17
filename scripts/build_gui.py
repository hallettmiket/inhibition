#!/usr/bin/env python3
"""
Purpose: the Home and Sweep pages, and the JSON the Sweep page polls (#63).
Author: Timothy Wu (with Claude Code)
Date: 2026-08-12
Input: --worklist <sweep_gaps_N.csv> (the campaign the sweep is executing)
Output: mdprio_reports/{index.html, sweep_state.json}
        (sweep.html belongs to scripts/sweep_combine.py -- one file, one owner)

@tt8804 (#63): "integrate the ranking and results pages into one uniform gui with
arrows showing each page as a distinct step ... Home page showing input target
and parameters, ranking results, sweep results, MD results", and: "show results
as pending what is sweeped and show sweep results as it goes".

TWO NEW PAGES, AND THE OTHER TWO GET THE SAME NAV. Home states what the run was
configured to do -- read from `config/target.yaml`, so the page cannot claim a
parameter the pipeline is not using. Sweep is the one that did not exist: 170
modes over two days, and until now nothing showed what was queued, what was
running, or what had come back.

IT POLLS JSON RATHER THAN REBUILDING. A full GUI rebuild is minutes and rewrites
5,700 pose assets; the sweep table is a few hundred rows. So this writes a small
`sweep_state.json` beside the page and the page re-renders from it every 30 s.
Re-running this script with --json-only refreshes that file in about a second,
which is what a watch loop should call.
"""

from __future__ import annotations

import argparse
import html
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import gui_shell as gs                        # noqa: E402
from shared import sweep_state as ss                      # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("build-gui")
OUT = rp.reports_dir()

#: The palette is the ranking view's, verbatim, because "uniform" is the whole
#: point of #63 -- a second palette that merely looks similar is what makes two
#: pages read as two tools.
BASE_CSS = """
:root{--ink:#10233f;--navy:#003087;--blue:#0072ce;--blue-pale:#e8f1fb;
 --rule:#d6dee8;--muted:#5b6b80;--paper:#fff;--raise:#f5f8fc;--rail:#fafcfe;
 --good:#0f7a54;--bad:#b3261e;--warn:#8a6d1f;
 --sans:"Helvetica Neue",Helvetica,Arial,system-ui,sans-serif;
 --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace}
:root[data-theme="dark"]{--ink:#dfe7f0;--navy:#8ab4e8;--blue:#6aa9e0;
 --blue-pale:#16283a;--rule:#25333f;--muted:#93a3b4;--paper:#0e151c;
 --raise:#16202a;--rail:#121b24;--good:#4fc4a0;--bad:#e08a70;--warn:#d0ae5a}
*{box-sizing:border-box}
html,body{height:100%}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
 font-size:14px;line-height:1.5;font-variant-numeric:tabular-nums;
 display:flex;flex-direction:column}
#topbar{display:flex;align-items:center;gap:10px;padding:6px 14px;min-height:38px;
 border-bottom:1px solid var(--rule);background:var(--raise);flex:0 0 auto}
#topbar h1{font:600 13px var(--sans);margin:0;color:var(--navy)}
.msep{flex:1}
.mbtn{font:600 11px var(--sans);padding:.25rem .6rem;border:1px solid var(--rule);
 border-radius:3px;background:var(--paper);color:var(--ink);cursor:pointer;
 text-decoration:none}
.mbtn:hover{background:var(--blue-pale)}
main{flex:1;min-height:0;overflow-y:auto;padding:18px 22px 40px}
h2{font:600 15px var(--sans);color:var(--navy);margin:22px 0 8px}
h2:first-child{margin-top:0}
p.note{font-size:12px;color:var(--muted);margin:.3rem 0 .8rem;max-width:80ch}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:.34rem .6rem;border-bottom:1px solid var(--rule);text-align:right}
th{font:600 10px var(--sans);color:var(--muted);text-transform:uppercase;
 letter-spacing:.04em;position:sticky;top:0;background:var(--raise);z-index:1}
th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left}
td{font-family:var(--mono)}
tr:hover td{background:var(--blue-pale)}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
 gap:10px 16px;margin:6px 0 14px}
.card{border:1px solid var(--rule);border-radius:4px;padding:9px 12px;background:var(--rail)}
.card b{display:block;font-family:var(--mono);font-size:1.35rem;line-height:1.1}
.card span{font-size:11px;color:var(--muted)}
.kv{display:grid;grid-template-columns:auto 1fr;gap:2px 14px;font-size:12.5px;
 max-width:70ch}
.kv dt{font-family:var(--mono);color:var(--muted)}
.kv dd{margin:0}
.tag{display:inline-block;padding:0 .4rem;border-radius:3px;font:600 10px var(--sans);
 text-transform:uppercase;letter-spacing:.03em}
.t-ok{background:#dff2e9;color:var(--good)}
.t-pending{background:#fdf3d8;color:var(--warn)}
.t-failed{background:#fbe3e0;color:var(--bad)}
.t-notsent{background:var(--raise);color:var(--muted)}
:root[data-theme="dark"] .t-ok{background:#12312a}
:root[data-theme="dark"] .t-pending{background:#2e2716}
:root[data-theme="dark"] .t-failed{background:#331c19}
.bar{display:block;height:5px;border-radius:3px;background:var(--rule);overflow:hidden;
 margin-top:8px}
.bar i{display:block;height:100%;background:var(--blue)}
.bar i.g{background:var(--good)}
.prog{display:flex;height:9px;border-radius:4px;overflow:hidden;border:1px solid var(--rule)}
.prog i{display:block;height:100%}
/* A NARROW SELECTOR AND A LARGE VIEWER -- the proportions of the MD results and
   ranking pages, which is what "look just like MD results" means in practice.
   The first attempt gave the table 1fr and the inspector 420px, so the thing
   being selected FROM dominated the thing being looked AT, and the plots ended
   up in a 420px column below the fold. Reversed: the rail is fixed and the
   viewer takes everything else. */
#two{display:grid;grid-template-columns:380px 1fr;gap:18px;align-items:start;
 min-height:0}
@media(max-width:1000px){#two{grid-template-columns:1fr}}
#tbl{max-height:calc(100vh - 250px);overflow-y:auto;border:1px solid var(--rule);
 border-radius:4px}
#tbl table{font-size:11.5px}
#tbl th{font-size:9px;padding:.3rem .4rem}
#tbl td{padding:.28rem .4rem}
/* The rail is 380px, so only the columns that decide a click survive: which
   mode, how it did, and how strongly. The rest is on the viewer side, where
   there is room to read it. */
#tbl td.dim,#tbl th.dim{display:none}
#side{min-width:0}
#side details{margin:.7rem 0}
#side summary{font:600 12.5px var(--sans);cursor:pointer;padding:.3rem 0}
#side .hint{font-weight:400;color:var(--muted);margin-left:.4rem;font-size:11px}
#side table{font-size:12.5px}
#side table td:first-child{color:var(--muted)}
#sgrid{display:grid;grid-template-columns:260px 1fr;gap:16px;align-items:start}
@media(max-width:1200px){#sgrid{grid-template-columns:1fr}}
tr.pick{cursor:pointer}
tr.pick.cur td{background:var(--blue-pale);font-weight:700}
/* Same callout the other pages use, so a caveat looks like a caveat everywhere. */
.warnbox{border-left:3px solid var(--warn);background:var(--rail);padding:.6rem .9rem;
 font-size:12px;margin:.7rem 0 .2rem;border-radius:0 3px 3px 0;max-width:95ch}
"""


def _page(title: str, current: str, counts: dict, body: str, extra_js: str = "",
          extra_css: str = "", head_js: str = "", tail_data: str = "") -> str:
    return f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>{BASE_CSS}{gs.CSS}{extra_css}</style>
{head_js}</head><body>
<div id="topbar"><h1>{html.escape(title)}</h1><span class="msep"></span>
<button class="mbtn" onclick="tt()">dark</button></div>
{gs.nav(current, counts)}
<main>{body}</main>
{tail_data}
<script>
function tt(){{const d=document.documentElement.getAttribute('data-theme')==='dark';
document.documentElement.setAttribute('data-theme',d?'light':'dark');}}
{extra_js}
</script></body></html>"""


def _fmt(v, nd=3, dash="—"):
    try:
        if v is None or v != v:
            return dash
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return dash


def home(counts: dict, s: dict, worklist: Path | None) -> str:
    from shared import target_config as tc
    c = tc.load()
    t, dk, sp = c.get("target", {}), c.get("docking", {}), c.get("splitting", {})
    sr, md, ch = c.get("sweep_rule", {}), c.get("md", {}), c.get("chemistry", {})
    run = c.get("run", {})
    st2 = sp.get("stage2", {})
    fam = (sr.get("scope", {}) or {}).get("families", {}) or {}

    def kv(pairs):
        return ('<dl class="kv">'
                + "".join(f"<dt>{html.escape(str(k))}</dt><dd>{v}</dd>"
                          for k, v in pairs) + "</dl>")

    floor = sr.get("floor")
    body = [
        '<p class="note">Every number on this page is read from '
        '<code>config/target.yaml</code>, so it cannot claim a setting the '
        'pipeline is not using.</p>',
        "<h2>Target</h2>",
        kv([("protein", f"{t.get('name')} — {t.get('domain')} domain"),
            ("receptor", f"{t.get('pdb')} (chemist-prepared, D0059)"),
            ("anchor", f"{t.get('anchor')} {t.get('anchor_atom')}"),
            ("run / topic", f"{run.get('topic')}"),
            ("tiers ranked", ", ".join(run.get("tiers", []) or ["all"]))]),
        "<h2>The funnel</h2>",
        '<div class="cards">'
        + "".join(f'<div class="card"><b>{v}</b><span>{k}</span></div>'
                  for k, v in [
                      # "RANKED", NOT "SCREENED", and the difference is 4,000
                      # molecules. The screen ran 5,697; the ranking carries the
                      # 1,634 that are T_4 (D0081). A card labelled "screened"
                      # showing the ranked count understates the screen by 4x and
                      # would be read as the funnel's first step.
                      ("molecules ranked (T_4)", counts.get("molecules", "—")),
                      ("binding modes ranked", counts.get("modes", "—")),
                      ("modes on this campaign", s.get("ok", 0) + s.get("pending", 0)
                                                 + s.get("failed", 0)),
                      ("swept ok", s.get("ok", 0)),
                      ("reach attack geometry", s.get("productive", 0)),
                  ]) + "</div>",
        '<p class="note">Each step is smaller than the last. A reader who cannot '
        'say why has found something worth asking about.</p>',
        "<h2>Pose generation and splitting</h2>",
        kv([("docking runs / molecule", dk.get("n_runs")),
            ("every pose persisted", dk.get("persist_all_poses")),
            ("stage 1", sp.get("stage1")),
            ("stage 2 (sub-split)", f"{st2.get('enabled')}, cut "
                                    f"{st2.get('cut_diameter_a')} Å, "
                                    f"max {st2.get('max_sub')} per mode")]),
        # THE CASCADE, WITH THE GATE AT EVERY STEP. Each stage is cheaper than
        # the next and asks a different question, and a reader who cannot see
        # the gates cannot tell why the counts fall the way they do.
        "<h2>The cascade</h2>",
        '<p class="note">Each stage is cheaper than the one after it and asks a '
        'different question. Docking ranks the <em>best case</em> — how good the '
        'pose would be if it held. Everything after that tests whether it '
        'does.</p>',
        '<table><thead><tr><th>stage</th><th>asks</th><th>gate</th>'
        '<th>cost each</th></tr></thead><tbody>'
        + "".join(f"<tr><td>{a}</td><td>{b}</td><td>{c}</td><td>{e}</td></tr>"
                  for a, b, c, e in [
            ("docking + NAC rank",
             "how good is this pose, if it held?",
             f"{dk.get('n_runs')} runs, split into modes; "
             f"enrichment &ge; {sr.get('budget_floor')}, "
             f"&ge; {sr.get('min_mode_poses')} poses/mode",
             "~25 s"),
            (f"triage sweep ({md.get('sweep_ps', 0)/1000:.0f} ns)",
             "is the pose stable at all?",
             f"max ligand RMSD &lt; {md.get('sweep_survivor_rmsd_nm')} nm",
             f"~{md.get('sweep_ps', 0)/1000*4.1:.0f} min"),
            (f"production MD ({md.get('production_ps', 0)/1000:.0f} ns)",
             "does it stay for a long time?",
             f"max ligand RMSD &lt; {md.get('sweep_survivor_rmsd_nm')} nm "
             "over the full run",
             "~4.5 h"),
            ("BPMD",
             "how hard is it to push out?",
             "promoted automatically; 3 replicates",
             "~1 h"),
        ]) + "</tbody></table>",
        "<h2>What earns a simulation</h2>",
        kv([("parameter", sr.get("parameter")),
            ("capture-validated floor",
             f"<b>{floor}</b>" if floor is not None
             else '<span style="color:var(--warn)">UNSET — the pilot has not run; '
                  'the value below is a spending rule, not a measured threshold</span>'),
            ("budget floor", sr.get("budget_floor")),
            ("min poses per mode", sr.get("min_mode_poses")),
            ("max depth / family", sr.get("max_depth")),
            ("scope", ", ".join(fam) or "UNSET"),
            ("worklist", html.escape(worklist.name) if worklist else "—")]),
        "<h2>Simulation and chemistry</h2>",
        kv([("triage sweep", f"{md.get('sweep_ps')} ps"),
            ("production", f"{md.get('production_ps')} ps"),
            ("salt", f"{md.get('salt_molar')} M"),
            ("docked species", ch.get("docked_species"))]),
    ]
    return _page("DWI covalent screen — Home", "index.html", counts, "".join(body))


def sweep_json(d, s, worklist: Path | None) -> str:
    cols = ["ident", "warhead_class", "mode_label", "class_rank", "enrichment",
            "sweep_state", "frac_attack_ready", "n_visits", "frac_in_window",
            "median_dist_a", "status"]
    rows = []
    for r in d.to_dict("records"):
        rows.append({k: (None if r.get(k) is None or r.get(k) != r.get(k)
                         else r.get(k)) for k in cols})
    return json.dumps({"generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                       "worklist": worklist.name if worklist else None,
                       "summary": s, "predicts": ss.predicts(d),
                       "rows": rows}, separators=(",", ":"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--worklist", default=None,
                    help="the campaign's worklist; named rather than guessed, "
                         "because two worklists can exist and disagree")
    ap.add_argument("--json-only", action="store_true",
                    help="refresh sweep_state.json only — what a watch loop calls")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    wl = Path(args.worklist) if args.worklist else None
    d = ss.state(wl)
    s = ss.summary(d)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "sweep_state.json").write_text(sweep_json(d, s, wl))
    if args.json_only:
        print(f"  {s} -> sweep_state.json")
        return

    from shared import mode_ranking as mr
    rk = mr.gather()
    counts = {"molecules": f"{rk.parent_ident.nunique():,}" if not rk.empty else "—",
              "modes": f"{len(rk):,}" if not rk.empty else "—"}
    nav_counts = {
        "modes.html": f"{len(rk):,} modes" if not rk.empty else "",
        "sweep.html": f"{s['ok']} ok · {s['pending']} pending",
    }
    (OUT / "index.html").write_text(home(counts, s, wl))
    # SWEEP.HTML IS NOT WRITTEN HERE ANY MORE. `scripts/sweep_combine.py` owns
    # it and builds it on the MD results shell, which is what @tt8804 asked for.
    # Both scripts used to write the file, so whichever ran last won -- running
    # build_gui after sweep_combine silently replaced the MD-shell page with the
    # earlier table layout, and the page "looked wrong again" for no reason
    # visible in either script. One file, one owner.
    #
    # A STAGE WITH NO RESULTS STILL HAS A PAGE. On a fresh topic only index.html
    # exists, so every other nav link and every bookmark 404s -- the server
    # answers "File not found", which reads as a broken deployment rather than
    # as a run that has not reached that stage. @tt8804 hit exactly that on
    # /sweep.html minutes after the topic was bumped.
    #
    # These are placeholders, and they are OVERWRITTEN by the real builders the
    # moment those have anything to show: `sweep_combine` owns sweep.html and
    # `mdprio_combine` owns modes.html and combined.html. Written only when
    # absent, so a rebuild during a live run never clobbers real results.
    n_placeholder = 0
    for href, label, why in (
        ("modes.html", "Ranking",
         "No modes ranked yet. Ranking runs after the docking + NAC screen "
         "finishes, and reads that screen's aggregates."),
        ("sweep.html", "Sweep",
         "No triage sweeps yet. The 8&nbsp;ns sweep runs on the modes the "
         "ranking selects, so it waits on the two stages before it."),
        ("combined.html", "MD results",
         "No 100&nbsp;ns runs yet. Only modes that hold under "
         "0.35&nbsp;nm through the triage sweep earn one."),
    ):
        p = OUT / href
        if p.is_file():
            continue
        p.write_text(_page(
            f"{label} — awaiting stage", href, nav_counts,
            f"<section class='card'><h2>{label}</h2><p class='muted'>{why}</p>"
            f"<p class='muted'>This page fills in on its own as the run "
            f"reaches this stage — the GUI rebuilds every few minutes.</p>"
            f"</section>"))
        n_placeholder += 1
    extra = f" (+{n_placeholder} awaiting-stage placeholders)" if n_placeholder else ""
    print(f"  index.html + sweep_state.json -> {OUT}{extra}")
    print(f"  {s}")


if __name__ == "__main__":
    main()
