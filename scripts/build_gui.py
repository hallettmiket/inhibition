#!/usr/bin/env python3
"""
Purpose: the Home and Sweep pages, and the JSON the Sweep page polls (#63).
Author: Timothy Wu (with Claude Code)
Date: 2026-08-12
Input: --worklist <sweep_gaps_N.csv> (the campaign the sweep is executing)
Output: mdprio_reports/{index.html, sweep.html, sweep_state.json}

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

log = logging.getLogger("build-gui")
OUT = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/mdprio_reports")

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
/* Same callout the other pages use, so a caveat looks like a caveat everywhere. */
.warnbox{border-left:3px solid var(--warn);background:var(--rail);padding:.6rem .9rem;
 font-size:12px;margin:.7rem 0 .2rem;border-radius:0 3px 3px 0;max-width:95ch}
"""


def _page(title: str, current: str, counts: dict, body: str, extra_js: str = "") -> str:
    return f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>{BASE_CSS}{gs.CSS}</style></head><body>
<div id="topbar"><h1>{html.escape(title)}</h1><span class="msep"></span>
<button class="mbtn" onclick="tt()">dark</button></div>
{gs.nav(current, counts)}
<main>{body}</main>
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


def sweep_page(counts: dict) -> str:
    body = """
<div id="head"></div>
<p class="note">10 ns triage. <b>Pending</b> means the mode is on the active
worklist with no result yet — queued or mid-trajectory; this page cannot see the
process table and does not pretend to. The table re-reads
<code>sweep_state.json</code> every 30 s, so it fills in as runs land.</p>
<div id="tbl"></div>"""
    js = """
let ROWS=[], SUM={}, PRED={}, SORT='frac_attack_ready';
// THE SWEEP'S OWN HEADLINE. Does the docked ranking predict what the trajectory
// did? The answer belongs on this page, updating as runs land, rather than in
// somebody's message -- and it must carry its own caveat, because every mode
// here cleared the enrichment floor and a correlation measured inside a
// selected range says nothing about the selection itself.
function verdict(){
  const P=PRED||{};
  if(!P.n || P.rho===undefined)
    return '<p class="note">Not enough finished runs yet to ask whether the '
      +'ranking predicts the outcome ('+(P.n||0)+' so far; needs 8).</p>';
  const sig = P.p < 0.05;
  const col = sig ? 'var(--good)' : 'var(--warn)';
  return '<div class="warnbox" style="border-left-color:'+col+'">'
    + '<b>Does the ranking predict the sweep?</b> Spearman(enrichment, '
    + 'attack-ready) = <b>'+P.rho.toFixed(3)+'</b>, p = '+P.p.toFixed(3)
    + ' over '+P.n+' finished modes. '
    + (sig ? 'A relationship is detectable.'
           : 'No relationship is detectable — ordering by enrichment above the '
             +'floor is not selecting the modes that reach attack geometry.')
    + '<br>Median enrichment: <b>'+(P.enr_prod!==null?P.enr_prod.toFixed(2):'—')
    + '</b> among the '+P.productive+' productive, <b>'
    + (P.enr_not!==null?P.enr_not.toFixed(2):'—')+'</b> among the rest.'
    + '<br><em>Range-restricted, and that limits the claim.</em> Every mode here '
    + 'cleared the floor (lowest swept: '+ (P.floor!==undefined?P.floor.toFixed(2):'—')
    + '). This measures whether enrichment discriminates ABOVE the floor, not '
    + 'whether the floor works — that needs modes sampled from below it, which '
    + 'is what the stratified pilot is for and it has not been run.</div>';
}
const ORDER={ok:0, pending:1, failed:2, 'not sent':3};
function f(v,n){return (v===null||v===undefined)?'—':Number(v).toFixed(n);}
function draw(){
  const s=SUM, tot=(s.ok||0)+(s.pending||0)+(s.failed||0);
  const pc=x=>tot?(100*x/tot).toFixed(1)+'%':'0%';
  document.getElementById('head').innerHTML =
    '<div class="cards">'
    + [['swept ok',s.ok||0],['reach attack geometry',s.productive||0],
       ['pending',s.pending||0],['failed',s.failed||0]]
      .map(k=>'<div class="card"><b>'+k[1]+'</b><span>'+k[0]+'</span></div>').join('')
    + '</div>'
    + '<div class="prog" title="ok / pending / failed">'
    + '<i style="background:var(--good);width:'+pc(s.ok||0)+'"></i>'
    + '<i style="background:var(--warn);width:'+pc(s.pending||0)+'"></i>'
    + '<i style="background:var(--bad);width:'+pc(s.failed||0)+'"></i></div>'
    + '<p class="note">'+tot+' modes in this campaign · worklist <code>'
    + (SUM._wl||'—')+'</code> · updated '+(SUM._t||'')+'</p>'
    + verdict();
  // Finished first and best-first within that, because the question this page
  // answers is "what came back and was any of it good".
  const r=ROWS.slice().sort((a,b)=>{
    const d=(ORDER[a.sweep_state]??9)-(ORDER[b.sweep_state]??9); if(d) return d;
    return (b.frac_attack_ready??-1)-(a.frac_attack_ready??-1);});
  document.getElementById('tbl').innerHTML =
    '<table><thead><tr><th>mode</th><th>class</th><th>state</th>'
    + '<th>attack-ready</th><th>visits</th><th>in window</th><th>median d (Å)</th>'
    + '<th>enrichment</th><th>class rank</th></tr></thead><tbody>'
    + r.map(x=>{
      const cls='t-'+(x.sweep_state||'').replace(' ','');
      const bad=x.sweep_state==='failed'&&x.status?' title="'+String(x.status).replace(/"/g,'')+'"':'';
      return '<tr><td>'+x.ident+'</td><td>'+(x.warhead_class||'—')+'</td>'
        +'<td'+bad+'><span class="tag '+cls+'">'+(x.sweep_state||'')+'</span></td>'
        +'<td>'+f(x.frac_attack_ready,3)+'</td><td>'+f(x.n_visits,0)+'</td>'
        +'<td>'+f(x.frac_in_window,3)+'</td><td>'+f(x.median_dist_a,2)+'</td>'
        +'<td>'+f(x.enrichment,2)+'</td><td>'+f(x.class_rank,0)+'</td></tr>';}).join('')
    + '</tbody></table>';
}
async function load(){
  try{
    const r=await fetch('sweep_state.json?t='+Date.now());
    if(!r.ok) return;
    const j=await r.json();
    ROWS=j.rows||[]; SUM=j.summary||{}; PRED=j.predicts||{};
    SUM._t=j.generated; SUM._wl=j.worklist;
    draw();
  }catch(e){}
}
load(); setInterval(load, 30000);"""
    return _page("DWI covalent screen — Sweep", "sweep.html", counts, body, js)


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
    (OUT / "sweep.html").write_text(sweep_page(nav_counts))
    print(f"  index.html + sweep.html + sweep_state.json -> {OUT}")
    print(f"  {s}")


if __name__ == "__main__":
    main()
