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


def sweep_page(counts: dict) -> str:
    """The table, plus one pose in the pocket -- same assets as the ranking view.

    The 3Dmol library and the receptor are vendored INTO the page and placed
    BEFORE the script that uses them. Both orderings have already produced a
    silently blank viewer in this project: a library loaded at the end of the
    body is undefined when the viewer looks for it, and a data element parsed
    after the code that reads it yields an empty string with no error.
    """
    body = """
<div id="head"></div>
<p class="note">10 ns triage. <b>Pending</b> means the mode is on the active
worklist with no result yet — queued or mid-trajectory; this page cannot see the
process table and does not pretend to. The table re-reads
<code>sweep_state.json</code> every 30 s, so it fills in as runs land.
<b>Click a row</b> to see the structure and the pose that was simulated.</p>
<div id="two">
  <div id="tbl"></div>
  <aside id="side">
    <div id="sname" class="note" style="margin:0 0 8px">select a mode</div>
    <div id="sgrid">
      <div>
        <img id="sstruct" class="pvstruct" alt="" style="display:none">
        <div class="pvbox" style="height:250px"><div id="pv"></div></div>
        <div class="pvctl">
          <label><input type="checkbox" id="pv-surf" checked onchange="redraw()">
            pocket surface</label>
        </div>
        <div id="sfacts"></div>
      </div>
      <div id="sright">
        <!-- The trajectory, at a size worth looking at. This is the panel the
             narrow-rail layout exists to make room for. -->
        <img id="splot" style="display:none;width:100%;border:1px solid var(--rule);
             border-radius:4px;background:#fff"
             alt="10 ns RMSD and warhead-sulfur distance">
        <div id="splotmiss" class="note" style="display:none"></div>
      </div>
    </div>
    <!-- The same three things the MD page shows, for the 10 ns run: the pose,
         the movie, and the trajectory plots. Details, so a 9 MB movie is
         fetched only when asked for. -->
    <details id="smovwrap" style="display:none"><summary>10 ns movie
      <span class="hint">ligand in yellow, CA-fitted — 126 frames</span></summary>
      <div class="pvbox" style="height:420px"><div id="smov"></div></div>
      <div class="pvctl"><button class="mbtn" onclick="playPause()"
        id="playbtn">play</button></div>
    </details>
  </aside>
</div>"""
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
    '<table><thead><tr><th>mode</th><th>state</th>'
    + '<th>ready</th><th>visits</th><th class="dim">in window</th>'
    + '<th class="dim">median d (Å)</th><th class="dim">enrichment</th>'
    + '<th class="dim">class rank</th></tr></thead><tbody>'
    + r.map(x=>{
      const cls='t-'+(x.sweep_state||'').replace(' ','');
      const bad=x.sweep_state==='failed'&&x.status?' title="'+String(x.status).replace(/"/g,'')+'"':'';
      return '<tr onclick="pick(\\''+x.ident+'\\')" class="pick'
        +(SEL===x.ident?' cur':'')+'"><td>'+x.ident.replace(/^t4_/,'')+'</td>'
        +'<td'+bad+'><span class="tag '+cls+'">'+(x.sweep_state||'')+'</span></td>'
        +'<td>'+f(x.frac_attack_ready,3)+'</td><td>'+f(x.n_visits,0)+'</td>'
        +'<td class="dim">'+f(x.frac_in_window,3)+'</td>'
        +'<td class="dim">'+f(x.median_dist_a,2)+'</td>'
        +'<td class="dim">'+f(x.enrichment,2)+'</td>'
        +'<td class="dim">'+f(x.class_rank,0)+'</td></tr>';}).join('')
    + '</tbody></table>';
}
// --- structure and pose, the same assets the ranking view draws -------------
let SEL=null, PDBC={};
function pick(id){
  SEL=id; draw();
  const x=ROWS.find(r=>r.ident===id); if(!x) return;
  const par=id.replace(/_m\\d+$/,''), mode=parseInt((id.match(/_m(\\d+)$/)||[0,'-1'])[1],10);
  document.getElementById('sname').innerHTML='<b>'+id+'</b>'
    +(x.mode_label&&x.mode_label!==String(mode)?' &middot; mode '+x.mode_label:'');
  const im=document.getElementById('sstruct');
  im.src='mode_thumbs/'+par+'.svg'; im.style.display='';
  document.getElementById('sfacts').innerHTML=
    '<table><tbody>'
    +[['state',x.sweep_state],['attack-ready',f(x.frac_attack_ready,3)],
      ['visits',f(x.n_visits,0)],['in window',f(x.frac_in_window,3)],
      ['median distance',f(x.median_dist_a,2)+' Å'],
      ['enrichment',f(x.enrichment,2)],['class rank',f(x.class_rank,0)]]
     .map(k=>'<tr><td>'+k[0]+'</td><td>'+k[1]+'</td></tr>').join('')
    +'</tbody></table>';
  load3d(par, mode);
  // The plot is an <img> and costs one request; the movie is ~9 MB and is
  // fetched only if the reader opens it. Both are named by MODE, not by
  // molecule -- two modes of one molecule are different trajectories.
  const pl=document.getElementById('splot'), ms=document.getElementById('splotmiss');
  // An absent figure SAYS it is absent. Hiding it silently is why "I don't see
  // the rmsd plots at all" was ambiguous between "not built yet" and "broken".
  pl.onerror=function(){ pl.style.display='none'; ms.style.display='';
    ms.innerHTML='<b>No trajectory figure for this mode yet.</b><br>'
      +'Assets are built per mode by <code>scripts/sweep_assets.py</code>; '
      +'this one has not been generated, or its 10 ns run did not finish.'; };
  pl.onload =function(){ pl.style.display=''; ms.style.display='none'; };
  pl.src='sweep_assets/'+id+'.png';
  const mw=document.getElementById('smovwrap');
  mw.style.display=''; mw.open=false; MOVID=id; MOVLOADED=false;
}
// --- the 10 ns movie, loaded on demand -------------------------------------
let MOVID=null, MOVLOADED=false, MV=null, PLAYING=false, MTIMER=null;
async function loadMovie(){
  if(MOVLOADED||!MOVID) return; MOVLOADED=true;
  const box=document.getElementById('smov');
  box.innerHTML='<div class="pvempty">loading 10 ns movie…</div>';
  try{
    const r=await fetch('sweep_assets/'+MOVID+'.pdb');
    if(!r.ok) throw new Error('no movie for this mode ('+r.status+')');
    const txt=await r.text(); const M=pvLib();
    box.innerHTML='';
    MV=M.createViewer(box,{backgroundColor:'#eef1f6'});
    MV.addModelsAsFrames(txt,'pdb');       // frames ARE animation here
    MV.setStyle({},{cartoon:{color:'#c3ccd8',opacity:0.45}});
    MV.setStyle({resn:'MOL'},{stick:{radius:0.20,colorscheme:'yellowCarbon'}});
    MV.setStyle({resi:[PV_CYS]},{stick:{radius:0.26,colorscheme:'default'}});
    MV.zoomTo({resn:'MOL'}); MV.zoom(0.55);
    requestAnimationFrame(function(){requestAnimationFrame(function(){
      MV.resize(); MV.render(); });});
    MV.render();
  }catch(e){
    box.innerHTML='<div class="pvempty"><b>No movie.</b><br>'
      +String(e.message||e)+'</div>';
  }
}
function playPause(){
  if(!MV) return;
  PLAYING=!PLAYING;
  document.getElementById('playbtn').textContent=PLAYING?'pause':'play';
  if(PLAYING){ MV.animate({loop:'forward',reps:0}); }
  else{ MV.stopAnimate(); }
}
async function load3d(par, mode){
  try{
    if(!PDBC[par]){
      const r=await fetch('mode_poses/'+par+'.pdb');
      if(!r.ok) throw new Error('no pose asset ('+r.status+')');
      PDBC[par]=await r.text();
    }
    mountPose('pv', PDBC[par], mode, document.getElementById('recpdb').textContent);
  }catch(e){
    // The reason, not a guess at it -- an empty box reads as a broken page.
    document.getElementById('pv').innerHTML=
      '<div class="pvempty"><b>No pose drawn.</b><br>'+String(e.message||e)+'</div>';
  }
}
document.addEventListener('toggle', function(e){
  if(e.target && e.target.id==='smovwrap' && e.target.open) loadMovie();
}, true);
function redraw(){ if(SEL){ const p=SEL.replace(/_m\\d+$/,'');
  const m=parseInt((SEL.match(/_m(\\d+)$/)||[0,'-1'])[1],10); load3d(p,m); } }

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
    # ORDER IS THE WHOLE POINT. The 3Dmol library goes in <head>, before any
    # code that names it; the receptor goes in the body BEFORE the script that
    # reads it. Both inversions have produced a silently blank viewer here.
    from shared import mode_ranking as mr
    from shared import pose_viewer as pv
    lib = (REPO / "scripts" / ".cache_3dmol-min.js")
    three = f"<script>{lib.read_text()}</script>" if lib.is_file() else ""
    rec = ("\n".join(l for l in mr.RECEPTOR.read_text().splitlines()
                     if l.startswith(("ATOM", "HETATM")))
           if mr.RECEPTOR.is_file() else "")
    tail = f'<pre id="recpdb" style="display:none">{html.escape(rec)}</pre>'
    vjs = pv.mount_js(mr.CYS_RESI, json.dumps(mr.pocket_residues()))
    return _page("DWI covalent screen — Sweep", "sweep.html", counts, body,
                 vjs + js, extra_css=pv.CSS, head_js=three, tail_data=tail)


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
    log.info("sweep page: viewer assets from %s", mr.B.name)
    print(f"  index.html + sweep.html + sweep_state.json -> {OUT}")
    print(f"  {s}")


if __name__ == "__main__":
    main()
