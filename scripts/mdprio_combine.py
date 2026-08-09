"""
Purpose: combine several per-molecule MD reports into one browsable page.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-08
Input: --candidates <ident...> (their reports under mdprio_reports/)
Output: 00_outputs/blacksmith/mdprio_reports/combined_<N>.html

@tt8804: *"combine them into one. I want to see the movies and plots."*

WHY NOT ONE CONCATENATED FILE. Each report is ~9.5 MB because the movie frames
and every plot are embedded as base64 -- self-contained by design, so a report can
be copied anywhere and still work. Four of them inlined into a single document is
~38 MB of HTML, which the browser must parse before it shows anything and which no
artefact host will accept.

So this builds a FRAME-BASED index instead: one page, a molecule picker across the
top with each molecule's headline numbers, and the selected report rendered whole
in an iframe beneath. Every movie and every plot is the real one from the original
report -- nothing is regenerated or downsampled -- and only the report being looked
at is loaded. The originals stay individually openable.

THE COMPARISON TABLE IS THE POINT. Flipping between four reports to remember which
molecule was 50% attack-ready is the work this is meant to remove, so the numbers
that decide the shortlist sit above the viewer where they can be read together:
the 10 ns sweep readings, the 100 ns engagement, and which binding MODE was
actually elevated -- that last one because a molecule promoted on its minority
mode is a different claim from one promoted on its dominant mode.
"""

from __future__ import annotations

import argparse
import glob
import html
import json
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                     # noqa: E402

log = logging.getLogger("mdprio-combine")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
REPORTS = B / "mdprio_reports"


def _sweep() -> pd.DataFrame:
    fs = sorted(glob.glob(str(B / "attack_sweep/attack_sweep_*.csv")),
                key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
    if not fs:
        return pd.DataFrame()
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d[(d.get("sweep_ps", 0) > 1000) & (d.status == "ok")]
    # the mode that was ELEVATED is the best-scoring surviving mode, which is how
    # the worker chose it -- not necessarily mode 0
    return d.sort_values("frac_attack_ready", ascending=False) \
            .drop_duplicates("parent_ident")


def _thumbs(idents) -> dict:
    """A small 2D depiction per molecule, for the selector.

    Base64 rather than inline SVG markup: RDKit emits an XML declaration and an
    HTML comment, and both break parsing when they land inside markup the browser
    is already mid-way through. Encoding removes the question entirely.
    """
    import base64
    import re as _re
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Draw, AllChem
    RDLogger.DisableLog("rdApp.*")
    smi = {}
    for d in ("04_t4_combinatorial/D4", "03_t3_reinvent/D3"):
        sub, stem = d.split("/")
        fs = sorted(glob.glob(f"/data/lab_vm/append_only/inhibition/{sub}/{stem}_*.parquet"),
                    key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))
        if not fs:
            continue
        fr = pd.read_parquet(fs[-1]).drop_duplicates("candidate_id")
        smi.update(dict(zip(fr.candidate_id, fr.canonical_smiles)))
    out = {}
    for i in idents:
        v = smi.get(i)
        m = Chem.MolFromSmiles(v) if isinstance(v, str) else None
        if m is None:
            continue
        AllChem.Compute2DCoords(m)
        d2 = Draw.rdMolDraw2D.MolDraw2DSVG(96, 64)
        d2.drawOptions().bondLineWidth = 1
        Draw.rdMolDraw2D.PrepareAndDrawMolecule(d2, m)
        d2.FinishDrawing()
        svg = _re.sub(r"<\?xml.*?\?>", "", d2.GetDrawingText(), flags=_re.S)
        svg = _re.sub(r"<!--.*?-->", "", svg, flags=_re.S)
        out[i] = "data:image/svg+xml;base64," + base64.b64encode(svg.encode()).decode()
    return out


def _classes() -> dict:
    """Warhead class per molecule, for the within-class ranking view."""
    out = {}
    for tier, score in (("T4", "conditional_eb"), ("T3", "enrichment_conditional")):
        fs = sorted(glob.glob(str(B / f"rank_v2/rank_v2_{tier}_{score}_*.csv")))
        if not fs:
            continue
        d = pd.read_csv(fs[-1]).drop_duplicates("parent_ident")
        out.update(dict(zip(d.parent_ident, d.warhead_class)))
    return out


def _md() -> pd.DataFrame:
    fs = glob.glob(str(B / "md_residence/*.csv"))
    if not fs:
        return pd.DataFrame()
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return d[d.get("production_ps", 0) >= 50000].drop_duplicates("ident", keep="last")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--title", default="100 ns candidates")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    sw, md, cls_of = _sweep(), _md(), _classes()
    thumbs = _thumbs(args.candidates)
    swi = sw.set_index("parent_ident") if not sw.empty else pd.DataFrame()
    mdi = md.set_index("ident") if not md.empty else pd.DataFrame()

    rows, tabs, missing = [], [], []
    for c in args.candidates:
        f = REPORTS / f"{c}.html"
        if not f.is_file():
            missing.append(c)
            log.warning("%s: no report at %s", c, f.name)
            continue
        s = swi.loc[c] if c in getattr(swi, "index", []) else None
        m = mdi.loc[c] if c in getattr(mdi, "index", []) else None

        def g(src, k, fmt="{:.3f}", dash="—"):
            if src is None or k not in src or pd.isna(src[k]):
                return dash
            try:
                return fmt.format(src[k])
            except Exception:                          # noqa: BLE001
                return str(src[k])

        mode = "—"
        if s is not None and "ident" in s and isinstance(s["ident"], str):
            mode = s["ident"].split("_m")[-1] if "_m" in s["ident"] else "0"
        rows.append(
            f"<tr><td class='id'>{html.escape(c)}</td>"
            f"<td>{mode}</td>"
            f"<td class='n'>{g(s,'frac_attack_ready')}</td>"
            f"<td class='n'>{g(s,'n_visits','{:.0f}')}</td>"
            f"<td class='n'>{g(s,'frac_in_window')}</td>"
            f"<td class='n'>{g(m,'explicit_frac_frames_engaged')}</td>"
            f"<td class='n'>{g(m,'explicit_ligand_rmsd_nm_max')}</td>"
            f"<td class='n'>{g(m,'explicit_ligand_rmsd_nm_n_independent','{:.1f}')}</td>"
            f"</tr>")
        tabs.append(c)

    if not tabs:
        raise SystemExit("no reports found for any requested candidate")

    btns = "".join(
        f"<button onclick=\"show('{html.escape(t)}')\" id='b_{html.escape(t)}'>"
        f"{html.escape(t)}</button>" for t in tabs)
    miss = (f"<p class='warn'>No report yet for: {', '.join(map(html.escape, missing))}"
            " — still running, or the trajectory is incomplete.</p>") if missing else ""

    # LEFT RAIL SELECTOR, RIGHT VIEWER (@tt8804) -- the same shape as the GUI's
    # ranking panel, because that is the layout the reading actually happens in:
    # you scan the list, click, and the pose/movie/plots replace themselves beside
    # it. A row of buttons across the top pushed the viewer below the fold and
    # made comparison a scroll.
    rows_html = []
    for k, t in enumerate(tabs):
        s_ = swi.loc[t] if t in getattr(swi, "index", []) else None
        m_ = mdi.loc[t] if t in getattr(mdi, "index", []) else None
        def g(src, key, fmt="{:.3f}"):
            if src is None or key not in src or pd.isna(src[key]):
                return "\u2014"
            try:
                return fmt.format(src[key])
            except Exception:                          # noqa: BLE001
                return str(src[key])
        ar = 0.0
        if s_ is not None and "frac_attack_ready" in s_ and not pd.isna(s_["frac_attack_ready"]):
            ar = float(s_["frac_attack_ready"])
        rmax = None
        if m_ is not None and "explicit_ligand_rmsd_nm_max" in m_ and not pd.isna(m_["explicit_ligand_rmsd_nm_max"]):
            rmax = float(m_["explicit_ligand_rmsd_nm_max"])
        held = rmax is not None and rmax < 1.2
        wcls = str(cls_of.get(t, "unclassified"))
        rows_html.append(
            f"<button class='row' data-cls=\"{html.escape(wcls)}\" data-eng='{ar:.6f}' "
            f"data-held='{1 if held else 0}' "
            f"id='b_{html.escape(t)}' onclick=\"show('{html.escape(t)}')\">"
            f"<span class='rk'>{k+1}</span>"
            + (f"<img class='thumb' alt='' src=\"{thumbs[t]}\">"
               if t in thumbs else "<span class='thumb'></span>")
            + f"<span class='body'>"
            f"<span class='l1'><span class='mid-id'>{html.escape(t)}</span>"
            f"<span class='eng'>{ar*10:.2f} ns</span></span>"
            f"<span class='l2'><span class='wc'>{html.escape(wcls)}</span>"
            f"<span class='meta'>{g(s_,'n_visits','{:.0f}')} visits &middot; "
            f"{g(m_,'explicit_ligand_rmsd_nm_max')} nm</span>"
            f"<span class='tag {'t-held' if held else 't-left'}'>"
            f"{'held' if held else 'left'}</span></span>"
            f"<span class='bar'><i style='width:{max(1.5,ar*100):.1f}%'></i></span>"
            f"</span></button>")

    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(args.title)}</title><style>
/* The per-molecule reports are light-only and use this exact palette
   (shared/report_theme.py). The shell inherits it so the frame and its
   contents read as one document rather than two. */
:root{{--ink:#10233f;--navy:#003087;--blue:#0072ce;--blue-pale:#e8f1fb;
 --rule:#d6dee8;--muted:#5b6b80;--paper:#fff;--raise:#f5f8fc;
 --good:#0f7a54;--bad:#b3261e;
 --sans:"Helvetica Neue",Helvetica,Arial,system-ui,sans-serif;
 --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace}}
*{{box-sizing:border-box}}
html,body{{height:100%}}
body{{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
 font-size:14px;line-height:1.5;font-variant-numeric:tabular-nums;
 display:flex;flex-direction:column}}
/* ONE STRIP. The title, both toggles and the hint share a single ~38px bar --
   roughly one selector row -- because everything below it is the actual work and
   a two-line masthead over a scrolling list is space the reader never gets back. */
#topbar{{display:flex;align-items:center;gap:8px;padding:6px 14px;min-height:38px;
 border-bottom:1px solid var(--rule);background:var(--raise);
 overflow-x:auto;white-space:nowrap;scrollbar-width:thin}}
h1{{margin:0;font-size:.86rem;font-weight:600;letter-spacing:-.01em;color:var(--navy);
 flex:none}}
.mbtn{{font:11.5px var(--sans);padding:3px 10px;flex:none;border:1px solid var(--rule);
 background:var(--paper);color:var(--muted);border-radius:99px;cursor:pointer}}
.mbtn.on{{background:var(--navy);border-color:var(--navy);color:#fff;font-weight:600}}
.mbtn:focus-visible{{outline:2px solid var(--blue);outline-offset:2px}}
.mhint{{font-size:11px;color:var(--muted);flex:none;margin-left:4px}}
.msep{{width:1px;height:16px;background:var(--rule);margin:0 2px;flex:none}}
.ohd{{font-family:var(--mono);font-size:.6rem;letter-spacing:.14em;text-transform:uppercase;
 font-weight:700;padding:11px 14px 7px;border-bottom:1px solid var(--rule);
 position:sticky;top:0;z-index:2}}
.o-held{{background:#e6f4ee;color:var(--good)}}
.o-left{{background:#fbeae8;color:var(--bad)}}
main{{flex:1;display:grid;grid-template-columns:376px 1fr;min-height:0}}
@media(max-width:880px){{main{{grid-template-columns:1fr;grid-template-rows:250px 1fr}}}}
#rail{{overflow-y:auto;border-right:1px solid var(--rule);background:var(--rail)}}
.chd{{font-family:var(--mono);font-size:.6rem;letter-spacing:.14em;text-transform:uppercase;
 color:var(--blue);font-weight:600;padding:10px 14px 6px;background:var(--raise);
 border-bottom:1px solid var(--rule);position:sticky;top:0;z-index:1}}
.row{{display:grid;grid-template-columns:22px 46px 1fr;gap:8px;align-items:start;
 width:100%;text-align:left;font:inherit;color:inherit;background:none;cursor:pointer;
 padding:9px 14px 8px;border:0;border-bottom:1px solid var(--rule)}}
.row:hover{{background:var(--blue-pale)}}
.row.on{{background:var(--blue-pale);box-shadow:inset 3px 0 0 var(--blue)}}
.row:focus-visible{{outline:2px solid var(--blue);outline-offset:-2px}}
.rk{{font-family:var(--mono);font-size:11px;font-weight:600;color:var(--muted);
 padding-top:3px}}
.thumb{{width:46px;height:32px;object-fit:contain;background:#fff;
 border:1px solid var(--rule);border-radius:3px;display:block}}
.body{{min-width:0;display:flex;flex-direction:column;gap:3px}}
.l1{{display:flex;align-items:baseline;justify-content:space-between;gap:8px}}
.mid-id{{font-family:var(--mono);font-size:12.5px;overflow:hidden;
 text-overflow:ellipsis;white-space:nowrap}}
.eng{{font-family:var(--mono);font-size:11.5px;font-weight:600;color:var(--navy);
 flex:none}}
.l2{{display:flex;align-items:center;gap:6px;flex-wrap:nowrap;min-width:0}}
.wc{{font-size:10.5px;color:var(--blue);white-space:nowrap;overflow:hidden;
 text-overflow:ellipsis;max-width:11ch}}
.meta{{font-size:10.5px;color:var(--muted);white-space:nowrap;overflow:hidden;
 text-overflow:ellipsis;flex:1}}
.tag{{font-size:9px;letter-spacing:.06em;text-transform:uppercase;font-weight:700;
 padding:1px 6px;border-radius:99px;flex:none}}
.t-held{{background:#e6f4ee;color:var(--good)}}
.t-left{{background:#fbeae8;color:var(--bad)}}
.bar{{height:3px;background:var(--rule);border-radius:2px;overflow:hidden;margin-top:2px}}
.bar i{{display:block;height:100%;background:var(--blue)}}
#viewer{{min-width:0;min-height:0;display:flex;flex-direction:column;background:var(--paper)}}
#vhead{{padding:9px 18px;border-bottom:1px solid var(--rule);background:var(--raise);
 display:flex;justify-content:space-between;align-items:center;gap:12px}}
#vname{{font-family:var(--mono);font-size:12.5px;font-weight:600;color:var(--navy)}}
#vhead a{{font-size:12px;color:var(--blue);text-decoration:none}}
#vhead a:hover{{text-decoration:underline}}
iframe{{flex:1;width:100%;border:0;background:var(--paper)}}
.tbtn{{margin-left:6px}}
:root[data-theme="dark"]{{--ink:#dfe7f0;--navy:#8ab4e8;--blue:#6aa9e0;
 --blue-pale:#16283a;--rule:#25333f;--muted:#93a3b4;--paper:#0e151c;
 --raise:#16202a;--rail:#121b24;--good:#4fc4a0;--bad:#e08a70}}
:root[data-theme="dark"] .thumb{{background:#fff}}
</style></head><body>
<div id="topbar">
 <h1 title="Pick a molecule on the left; its pose, movie and plots load on the right.">{html.escape(args.title)}</h1>
 <span class="msep"></span>
 <button id="m-all" class="mbtn on" onclick="setMode('all')">all classes</button>
 <button id="m-cls" class="mbtn" onclick="setMode('cls')">by warhead class</button>
 <span class="msep"></span>
 <button id="o-mix" class="mbtn on" onclick="setSplit(0)">combined</button>
 <button id="o-spl" class="mbtn" onclick="setSplit(1)">split held / left</button>
 <span class="mhint" id="mhint"></span>
 <button id="theme" class="mbtn tbtn" onclick="toggleTheme()" title="light / dark">dark</button>
</div>
<main>
 <div id="rail">{''.join(rows_html)}</div>
 <div id="viewer">
  <div id="vhead"><span id="vname">&mdash;</span>
   <a id="vopen" href="#" target="_blank" rel="noopener">open full report &#8599;</a></div>
  <iframe id="v" title="molecule report" src="{html.escape(tabs[0])}.html"></iframe>
 </div>
</main>
<script>
var RAIL=document.getElementById('rail');
var ROWS=Array.prototype.slice.call(RAIL.querySelectorAll('.row'));
var MODE='all', SPLIT=0;
function renumber(l){{l.forEach(function(b,i){{b.querySelector('.rk').textContent=i+1}});}}
function hdr(cls,txt){{var h=document.createElement('div');h.className=cls;h.textContent=txt;
  RAIL.appendChild(h);}}
function byEng(a,b){{return parseFloat(b.dataset.eng)-parseFloat(a.dataset.eng)}}
function layoutGroup(rows){{
  if(MODE==='all'){{rows.forEach(function(b){{RAIL.appendChild(b)}});renumber(rows);return;}}
  var g={{}};
  rows.forEach(function(b){{(g[b.dataset.cls]=g[b.dataset.cls]||[]).push(b)}});
  Object.keys(g).sort(function(x,y){{return byEng(g[x][0],g[y][0])}}).forEach(function(n){{
    hdr('chd',n+'  ('+g[n].length+')');
    g[n].forEach(function(b){{RAIL.appendChild(b)}}); renumber(g[n]);
  }});
}}
function relayout(){{
  RAIL.querySelectorAll('.chd,.ohd').forEach(function(h){{h.remove()}});
  var all=ROWS.slice().sort(byEng);
  if(!SPLIT){{ layoutGroup(all); }}
  else{{
    var held=all.filter(function(b){{return b.dataset.held==='1'}});
    var gone=all.filter(function(b){{return b.dataset.held!=='1'}});
    if(held.length){{hdr('ohd o-held','held the pocket  ('+held.length+')'); layoutGroup(held);}}
    if(gone.length){{hdr('ohd o-left','dissociated  ('+gone.length+')'); layoutGroup(gone);}}
  }}
  var bits=[ROWS.length+' molecules'];
  bits.push(MODE==='all'?'one ranking across all classes'
    :'ranked within warhead class \u2014 cross-class comparison is biased (#47)');
  if(SPLIT) bits.push('held and dissociated shown separately');
  document.getElementById('mhint').textContent=bits.join(' \u00b7 ');
}}
function setMode(m){{MODE=m;
  document.getElementById('m-all').classList.toggle('on',m==='all');
  document.getElementById('m-cls').classList.toggle('on',m==='cls');
  relayout();}}
function setSplit(v){{SPLIT=v;
  document.getElementById('o-mix').classList.toggle('on',!v);
  document.getElementById('o-spl').classList.toggle('on',!!v);
  relayout();}}
function applyTheme(doc){{
  try{{ doc.documentElement.setAttribute('data-theme',
        document.documentElement.getAttribute('data-theme')||'light'); }}catch(e){{}}
}}
function toggleTheme(){{
  var cur=document.documentElement.getAttribute('data-theme')==='dark'?'light':'dark';
  document.documentElement.setAttribute('data-theme',cur);
  document.getElementById('theme').textContent=cur==='dark'?'light':'dark';
  try{{localStorage.setItem('cat-theme',cur)}}catch(e){{}}
  var f=document.getElementById('v');
  if(f&&f.contentDocument) applyTheme(f.contentDocument);
}}
(function(){{
  var saved='light';
  try{{saved=localStorage.getItem('cat-theme')||'light'}}catch(e){{}}
  document.documentElement.setAttribute('data-theme',saved);
  document.addEventListener('DOMContentLoaded',function(){{
    document.getElementById('theme').textContent=saved==='dark'?'light':'dark';
  }});
}})();
function show(t){{
  var f=document.getElementById('v');
  f.onload=function(){{applyTheme(f.contentDocument)}};
  f.src=t+'.html';
  document.getElementById('vname').textContent=t;
  document.getElementById('vopen').href=t+'.html';
  document.querySelectorAll('.row').forEach(function(b){{b.classList.remove('on')}});
  var el=document.getElementById('b_'+t); if(el){{el.classList.add('on');}}
}}
relayout();
show({json.dumps(tabs[0])});
</script></body></html>"""

    dest = sout.Topic("blacksmith", "mdprio_reports").write("combined", ".html")
    dest.write_text(page)
    # Also drop a stable name beside the reports so the iframes resolve relatively.
    side = REPORTS / "combined.html"
    side.write_text(page)
    print(f"\n  {len(tabs)} reports combined -> {side}")
    print(f"  versioned copy: {dest}")


if __name__ == "__main__":
    main()
