#!/usr/bin/env python3
"""
Purpose: the per-LIGAND ranking page — one row per molecule, not per mode.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-28
Input: the run's rank_v2 engagement table (topic from config)
Output: <reports>/ligands.html

@tt8804: "can i see rank by mol in gui". The Ranking page is per MODE -- 327,167
rows for nac_v6 -- because a mode is the thing that gets simulated. But the thing
you order, buy and test is a MOLECULE, and nothing showed that.

BOTH AGGREGATIONS ARE SHOWN, SIDE BY SIDE, because they answer different
questions and picking one silently would answer the wrong one (D0098):

  best  -- can this ligand reach attack geometry AT ALL. Immune to how many
           groups it happens to have, and therefore to docking depth. But it is
           a maximum over a noisy score, and argmax selection is separately
           measured as the worst rule available at the pose level (6.7% crystal
           recovery against 33.3% for the medoid of the well-anchored quartile).
  mean  -- how well does it engage TYPICALLY. Penalises a molecule that can only
           reach attack geometry one way in twenty -- right if you believe the
           mode population, wrong if you do not, because `n_modes` grows with
           docking depth and never saturates (D0092, b = +0.69).

`n_modes` is therefore printed beside the mean on every row: a mean over a
denominator that moves with runtime is comparable only at fixed depth.

NOT A CLAIM THAT ANY OF THESE BIND. The ordering is reachability of attack
geometry, which is a precondition for covalent chemistry and not evidence of it.
`rank_validated` is False for this run as for every other.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import glob
import html
import logging
import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import engagement_rank as er           # noqa: E402
from shared import run_paths as rp                 # noqa: E402
from shared import target_config as tc             # noqa: E402

log = logging.getLogger("ligand-page")

#: THE SHELL IS THE PROJECT'S, NOT THIS PAGE'S. The first version of this page
#: carried its own palette and its own chrome, so it read as a different
#: instrument bolted onto the GUI -- no stepper, no "how this works", a different
#: dark toggle, different type. `gui_shell` exists precisely so a page cannot
#: drift from the others (#63), and the palette below is `mode_ranking._TPL`'s
#: verbatim so the two ranking views look like one product.
_PALETTE = """
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
#topbar{display:flex;align-items:center;gap:8px;padding:6px 14px;min-height:38px;
 border-bottom:1px solid var(--rule);background:var(--raise);
 overflow-x:auto;white-space:nowrap;scrollbar-width:thin}
h1{margin:0;font-size:.86rem;font-weight:600;letter-spacing:-.01em;color:var(--navy)}
.msep{width:1px;height:16px;background:var(--rule);flex:none}
.mbtn{font:inherit;font-size:12px;padding:3px 9px;border:1px solid var(--rule);
 background:var(--paper);color:var(--ink);border-radius:4px;cursor:pointer;
 text-decoration:none}
.mbtn:hover{border-color:var(--blue);color:var(--blue)}
.mhint{font-size:12px;color:var(--muted)}
input,select{font:inherit;font-size:12px;padding:3px 7px;border:1px solid var(--rule);
 background:var(--paper);color:var(--ink);border-radius:4px}
main{flex:1;display:grid;grid-template-columns:minmax(0,1fr) 460px;
 gap:1px;background:var(--rule);overflow:hidden}
#left{background:var(--paper);overflow:auto;padding:12px 14px 40px}
#right{background:var(--rail);display:flex;flex-direction:column;overflow:hidden}
#gl{flex:1;min-height:300px;position:relative}
#vhead{padding:9px 12px;border-bottom:1px solid var(--rule);background:var(--raise)}
#vhead b{color:var(--navy)}
#vmeta{font-size:12px;color:var(--muted);margin-top:3px}
#vfoot{padding:8px 12px;border-top:1px solid var(--rule);font-size:11.5px;
 color:var(--muted)}
.thumb{width:100%;max-height:150px;object-fit:contain;background:var(--paper);
 border-top:1px solid var(--rule)}
tbody tr.sel{background:var(--blue-pale);box-shadow:inset 3px 0 0 var(--blue)}
@media (max-width:1180px){main{grid-template-columns:1fr;grid-template-rows:1fr 60vh}}
.lede{color:var(--muted);font-size:12.5px;max-width:112ch;margin:0 0 12px}
.lede b{color:var(--ink)}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:5px 9px;border-bottom:1px solid var(--rule);text-align:right;white-space:nowrap}
th:nth-child(2),td:nth-child(2),th:last-child,td:last-child{text-align:left}
thead th{position:sticky;top:0;background:var(--rail);cursor:pointer;font-size:11px;
 font-weight:600;color:var(--navy);border-bottom:1px solid var(--rule)}
thead th:hover{color:var(--blue)}
tbody tr:nth-child(even){background:var(--rail)}
tbody tr:hover{background:var(--blue-pale)}
code{font-family:var(--mono);font-size:12px}
#vstats{width:100%;border-collapse:collapse;font-size:12px;margin:0}
#vstats td{padding:4px 12px;border-bottom:1px solid var(--rule)}
#vstats td:first-child{color:var(--muted);width:52%}
#vstats td:last-child{text-align:right;font-family:var(--mono);
  font-variant-numeric:tabular-nums}
#vstats tr.key td{font-weight:600;color:var(--ink);background:var(--raise)}
#vstats tr.sub td:first-child{padding-left:24px}
#vstats tr.gap td{border-bottom:2px solid var(--rule)}
.smi{max-width:320px;overflow:hidden;text-overflow:ellipsis;display:inline-block;
 vertical-align:bottom;color:var(--muted);font-family:var(--mono);font-size:11px}
.tw{overflow-x:auto;border:1px solid var(--rule);border-radius:5px}
"""

# RAW STRING, AND IT HAS TO BE. This block is JavaScript passing through a
# Python literal, and a non-raw """...""" REINTERPRETS the escapes on the way
# through: `\n` inside the regex `/MODEL[^\n]*\n/` became an actual newline,
# which split the regex literal across three lines and made the WHOLE script
# a syntax error. Nothing in the page said so -- it served 200, rendered
# perfectly, and simply had no behaviour: the family filter did nothing and
# clicking a row did nothing, because `show` was never defined. A guard is
# below (`_assert_js_parses`).
_JS = r"""
const tb=document.querySelector('tbody');
let dir={};
document.querySelectorAll('thead th').forEach((th,i)=>th.onclick=()=>{
  dir[i]=!dir[i];
  const rows=[...tb.rows];
  rows.sort((a,b)=>{
    const x=a.cells[i].dataset.v??a.cells[i].textContent, y=b.cells[i].dataset.v??b.cells[i].textContent;
    const nx=parseFloat(x), ny=parseFloat(y);
    const c=(!isNaN(nx)&&!isNaN(ny))?nx-ny:String(x).localeCompare(String(y));
    return dir[i]?c:-c;});
  rows.forEach(r=>tb.appendChild(r));});
const q=document.getElementById('q'), f=document.getElementById('f');
function filt(){const s=q.value.toLowerCase(),c=f.value;
  [...tb.rows].forEach(r=>{r.style.display=(r.textContent.toLowerCase().includes(s)
    &&(!c||r.dataset.c===c))?'':'none';});}
q.oninput=filt; f.onchange=filt;
let V=null, CACHE={}, SEL=null;
function lib(){ return window.$3Dmol || null; }
const f3 = v => (v === null || v === undefined) ? '\u2014' : Number(v).toFixed(3);
const pct = v => (v === null || v === undefined) ? '\u2014'
                 : (100 * Number(v)).toFixed(1) + '%';
/* THE PANEL NAMES WHICH ROW DID THE RANKING. A table of six summary statistics
   in which one of them is the sort key and the rest are not is unreadable
   unless it says so -- and here the distinction decides which molecule wins:
   ranking on `best` instead of `mean` changes 48 of the top 50. */
function renderStats(ident){
  const b = document.querySelector('#vstats tbody');
  const s = STATS[ident];
  if (!s){ b.innerHTML = ''; return; }
  const row = (k, v, cls) => '<tr class="' + (cls || '') + '"><td>' + k
                             + '</td><td>' + v + '</td></tr>';
  let h = row('<b>fraction in the NAC window</b> &middot; ranks this table',
              (100 * s['ge0.5'] / s.n).toFixed(1) + '%', 'key');
  h += row('mean rank_score &middot; tiebreak', f3(s.mean));
  h += row('median', f3(s.median));
  h += row('best mode', f3(s.best));
  h += row('90th percentile', f3(s.p90));
  h += row('worst mode', f3(s.worst), 'gap');
  h += row('modes total', s.n);
  h += row('<span style="color:var(--muted)">anchor_quality, the geometric '
           + 'criterion &mdash; support is excluded</span>', '', 'sub');
  CUTS.forEach(c => {
    const n = s['ge' + c];
    const nm = {0.25:'&ge; 0.25 &mdash; could be in window',
                0.5:'&ge; 0.50 &mdash; <b>inside the window</b>',
                0.7:'&ge; 0.70 &mdash; well inside'}[c] || ('&ge; ' + c);
    h += row(nm, n + '  <span style="color:var(--muted)">('
             + (100 * n / s.n).toFixed(1) + '%)</span>', 'sub');
  });
  h += row('modes H-bonding Ser114/115', pct(s.supported), 'gap');
  b.innerHTML = h;
}
async function show(ident, mode, label){
  document.querySelectorAll('tbody tr').forEach(r =>
    r.classList.toggle('sel', r.dataset.i === ident));
  SEL = ident;
  document.getElementById('vname').textContent = ident;
  document.getElementById('vmeta').textContent = label;
  document.getElementById('vthumb').src = 'mode_thumbs/' + ident + '.svg';
  renderStats(ident);
  const gl = document.getElementById('gl');
  const M = lib();
  if (!M){ gl.innerHTML = '<p style="padding:14px">3Dmol did not load. '
    + 'Asset: <code>3Dmol-min.js</code> beside this page.</p>'; return; }
  try{
    if (!CACHE[ident]){
      const r = await fetch('mode_poses/' + ident + '.pdb');
      if (!r.ok) throw new Error('HTTP ' + r.status);
      CACHE[ident] = await r.text();
    }
    if (!REC){
      const r = await fetch('receptor.pdb');
      if (!r.ok) throw new Error('receptor HTTP ' + r.status);
      REC = await r.text();
    }
    if (!V) V = M.createViewer(gl, {backgroundColor:'#eef1f6'});
    V.clear();
    V.addModel(REC, 'pdb');
    V.setStyle({}, {cartoon:{color:'#c3ccd8', opacity:0.5}});
    /* Cys113 in element colours -- it is the atom the whole screen aims at --
       with a sphere on SG, which is easy to lose at stick radius. */
    V.setStyle({resi:[113]}, {stick:{radius:0.28, colorscheme:'default'},
                              cartoon:{color:'#c3ccd8', opacity:0.5}});
    V.addStyle({resi:[113], atom:'SG'}, {sphere:{radius:0.62}});
    /* THE MODE IS READ FROM THE MODEL RECORD, never counted by position --
       counting is #53, and it is how a molecule's mode 4 came to be drawn as
       its mode 0. */
    const blocks = CACHE[ident].split('ENDMDL').filter(b => b.indexOf('MODEL') >= 0);
    let drew = false;
    blocks.forEach(b => {
      const m = /MODEL\s+(-?\d+)/.exec(b);
      if (m && parseInt(m[1], 10) === mode){
        V.addModel(b.replace(/MODEL[^\n]*\n/, ''), 'pdb');
        V.setStyle({model:-1}, {stick:{radius:0.20, colorscheme:'default'}});
        drew = true;
      }
    });
    if (!drew){
      document.getElementById('vfoot').textContent =
        'mode ' + mode + ' is not in this molecule\'s pose file — showing the receptor only.';
    } else {
      document.getElementById('vfoot').textContent =
        'best mode ' + mode + ' · grey cartoon is 3IKD · yellow sphere is Cys113 SG';
    }
    V.zoomTo({resi:[113]}); V.zoom(0.55); V.render();
  }catch(e){
    gl.innerHTML = '<p style="padding:14px"><b>No pose drawn.</b><br>'
      + String(e && e.message ? e.message : e)
      + '<br><span style="color:var(--muted)">asset: mode_poses/' + ident + '.pdb</span></p>';
  }
}
let REC = null;
function toggleTheme(){
  const d = document.documentElement.getAttribute('data-theme') === 'dark';
  document.documentElement.setAttribute('data-theme', d ? 'light' : 'dark'); }
"""





def _assert_js_parses(js: str, src: str) -> None:
    """Refuse to write a page whose script cannot parse.

    THE FAILURE THIS EXISTS FOR IS INVISIBLE. A JavaScript syntax error costs
    the page every behaviour it has -- handlers, click targets, the viewer --
    while the HTML still renders perfectly and the server still returns 200. The
    only symptom is that nothing does anything, which reads like an unfinished
    feature rather than a broken one.

    Two checks, and NEITHER CAN PASS VACUOUSLY:

    1. Escapes survived the Python literal. The escapes are read out of THIS
       FILE'S SOURCE TEXT, not out of `js` -- reading them from `js` is the
       vacuous version, because dropping the `r` prefix is exactly what removes
       them from `js`, and a check driven by an empty list passes for free
       (`docs/how_this_project_breaks.md`, disguise #4).
    2. No character class is split across lines -- the visible signature of an
       escape that became a real newline, and the shape of the original bug.
    """
    want = set(re.findall(r"\\[nrtsdwbSDWB]", src))
    missing = sorted(e for e in want if e not in js)
    if missing:
        raise SystemExit(
            f"the emitted script lost the escape(s) {missing} on the way through "
            f"the Python literal -- `_JS` needs its r-prefix. Emitting it would "
            f"produce a page that renders and does nothing.")
    if not want:
        raise SystemExit("_assert_js_parses found no escapes to check in the "
                         "source; the guard is vacuous and must be repaired.")
    for i, line in enumerate(js.splitlines(), 1):
        if line.rstrip().endswith("[^"):
            raise SystemExit(f"script line {i} ends with an open character "
                             f"class -- a regex literal is split across lines.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=None)
    # THE DEFAULT COMES FROM THE CONFIG, NOT FROM HERE. Which aggregation the
    # page uses is a project decision with a documented rationale beside it
    # (`ranking.ligand_agg`), and a second copy of it in an argparse default is
    # a pin that cannot announce it has gone stale -- disguise #3.
    ap.add_argument("--sort", default=tc.load()["ranking"]["ligand_agg"],
                    choices=("best", "mean", "median", "fraction_above"))
    ap.add_argument("--cutoff", type=float,
                    default=tc.load()["ranking"].get("ligand_cutoff"),
                    help="rank_score a mode must clear to count toward "
                         "fraction_above; REQUIRED with --sort fraction_above "
                         "and refused otherwise. There is no defensible "
                         "default: the curves are smooth and every cutoff from "
                         "0.05 to 0.60 leaves ~100%% of molecules with at least "
                         "one mode above it, so it orders molecules without "
                         "selecting them. That is why `mean` is the config "
                         "default -- it needs no threshold at all.")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if a.sort == "fraction_above" and a.cutoff is None:
        raise SystemExit("--sort fraction_above needs an explicit --cutoff; "
                         "no value is defensible from the curves, so it has to "
                         "be a stated choice rather than a default.")

    topic = a.topic or rp.topic()
    fs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_v2/rank_v2_T4_{topic}_engagement_*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit(f"no engagement ranking for topic {topic!r}")
    d = pd.read_csv(fs[-1])
    # PREFER THE SUPPORT-WEIGHTED SCORE WHEN IT EXISTS, and say which one the
    # page is showing. Two scores that both order modes are indistinguishable in
    # a table; a page that does not name its own is not comparable with one that
    # used the other.
    rs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_score_{topic}/rank_score_*.csv")),
                key=os.path.getmtime)
    score_col, score_label = "engagement", "engagement (geometry only)"
    if rs:
        d = pd.read_csv(rs[-1])
        score_col, score_label = "rank_score", "rank_score (geometry x support)"
    d["family"] = d.warhead_class.map(tc.family_of())

    # one row per MOLECULE, both aggregations
    # KEYED ON parent_ident DIRECTLY. Renaming it to `ident` collided with the
    # mode-level `ident` already in the frame, and pandas groups on the ambiguous
    # name rather than raising anything a reader would understand.
    cut_col = tc.load()["ranking"].get("ligand_cut_col")
    lig = er.rank_ligands(d, how=a.sort, ligand_key="parent_ident",
                          cutoff=(a.cutoff if a.sort == "fraction_above" else None),
                          score_col=score_col,
                          cut_col=(cut_col if a.sort == "fraction_above" else None))
    meta = (d.groupby("parent_ident")
              .agg(warhead_class=("warhead_class", "first"),
                   family=("family", "first"),
                   smiles=("smiles", "first"),
                   QED=("QED", "first"),
                   supported=("support", lambda v: float((v > 1.0).mean())
                              if v.notna().any() else float("nan")),
                   best_mode=(score_col, "idxmax"))
              .reset_index())
    meta["best_mode"] = d.loc[meta.best_mode, "mode"].values
    t = lig.merge(meta, on="parent_ident", how="left").rename(
        columns={"parent_ident": "ident"})
    t = t.sort_values("ligand_engagement", ascending=False).reset_index(drop=True)
    t.insert(0, "rank", range(1, len(t) + 1))

    # PER-MOLECULE DISTRIBUTION, for the viewer panel. The ranking column is one
    # number off a distribution of `n_modes` values, and a single number cannot
    # say whether it came from a tight population or one excellent mode dragging
    # a mediocre crowd. Both shapes are common here and they are NOT the same
    # molecule: nac_v6's top-50 by `mean` and top-50 by `best` share two members.
    # So the panel shows the distribution, and marks which row did the ranking.
    #
    # THE THRESHOLD COUNTS ARE DESCRIPTIVE, NOT A GATE. They are here because
    # "how many of its modes are actually good" is the question a reader asks
    # next, and 0.4/0.6/0.8 span the usable range. Nothing selects on them --
    # no cutoff separates (every one from 0.05 to 0.60 leaves ~100% of molecules
    # with a qualifying mode), which is exactly why `mean` is the ranking.
    STAT_CUTS = (0.25, 0.5, 0.7)   # necessary / sufficient / well inside
    gs = d.groupby("parent_ident")[score_col]
    cut_stat = cut_col if (cut_col and cut_col in d.columns) else score_col
    st = pd.DataFrame({"mean": gs.mean(), "median": gs.median(), "best": gs.max(),
                       "p90": gs.quantile(0.90), "worst": gs.min(),
                       "n": gs.size()})
    for c in STAT_CUTS:
        st[f"ge{c}"] = d[d[cut_stat] >= c].groupby("parent_ident")[cut_stat] \
                        .size().reindex(st.index, fill_value=0)
    st["supported"] = (d.assign(_s=d["support"] > 1.0)
                       .groupby("parent_ident")["_s"].mean()
                       if "support" in d.columns else float("nan"))
    stats = {i: {k: (None if pd.isna(v) else float(v)) for k, v in r.items()}
             for i, r in st.iterrows()}
    log.info("panel stats for %d molecules, cuts %s", len(stats), STAT_CUTS)
    log.info("%s: %d molecules from %d modes", topic, len(t), len(d))

    # THE TABLE CARRIES WHAT YOU SORT ON; the panel carries what you read.
    # Eight summary statistics per row is not a ranking, it is a spreadsheet --
    # and one of them (`above cut`) was empty on every row, because it is
    # `fraction_above`'s numerator and the config ranks on `mean`. A column
    # that is blank for the whole table teaches the reader to distrust the
    # others. The full distribution now renders in the viewer pane for the one
    # molecule being looked at, which is where a distribution is legible.
    cols = [("rank", "rank"), ("ident", "molecule"),
            ("ligand_engagement", "score"), ("n_modes", "modes"),
            ("warhead_class", "warhead class"), ("QED", "QED"),
            ("smiles", "SMILES")]
    for k, _ in cols:
        if k not in t.columns:
            t[k] = ""
    body = []
    for r in t.itertuples():
        fam = r.family if isinstance(r.family, str) else ""
        cells = []
        for k, _ in cols:
            v = getattr(r, k, "")
            if k == "smiles":
                cells.append(f'<td><span class="smi" title="{html.escape(str(v))}">'
                             f'{html.escape(str(v))}</span></td>')
            elif isinstance(v, float):
                cells.append(f'<td data-v="{v}">{v:.3f}</td>')
            elif k == "ident":
                cells.append(f'<td><code>{html.escape(str(v))}</code></td>')
            else:
                cells.append(f'<td data-v="{v}">{html.escape(str(v))}</td>')
        bm = "" if r.best_mode == "" else int(r.best_mode)
        lab = (f"score {getattr(r,'ligand_engagement',float('nan')):.3f} · "
               f"best mode {bm} · {getattr(r,'n_modes','?')} modes")
        body.append(
            f'<tr data-c="{html.escape(fam)}" data-i="{html.escape(str(r.ident))}" '
            f'onclick="show(\'{html.escape(str(r.ident))}\',{bm if bm != "" else -1},'
            f'\'{html.escape(lab)}\')">' + "".join(cells) + "</tr>")

    fams = sorted({x for x in t.family.dropna().unique()})
    opts = "".join(f'<option value="{html.escape(f)}">{html.escape(f)}</option>' for f in fams)
    head = "".join(f"<th>{html.escape(lbl)}</th>" for _, lbl in cols)
    from shared import gui_shell as GS
    counts = GS.step_counts()
    # THE STEPPER SAYS "Ranking", because that is what this is a view of. It is
    # not a new pipeline stage and must not draw as one -- the chevrons between
    # steps are the funnel, and a fifth arrow would claim a stage that does not
    # exist.
    nav = GS.nav("modes.html", counts)
    # `json.dumps` and not str(): a Python dict repr is not JSON (None, True),
    # and it would parse as a ReferenceError inside the script tag -- silently,
    # taking every handler on the page with it, which is the bug the guard
    # below exists for.
    stats_json = json.dumps(stats, allow_nan=False)
    cuts_json = json.dumps(list(STAT_CUTS))
    _assert_js_parses(_JS, pathlib.Path(__file__).read_text())
    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{topic} — ligands</title>
<style>{_PALETTE}
{GS.CSS}</style></head><body>
<div id="topbar">
 <h1 title="One row per molecule. Click a header to sort.">{topic} &mdash; ligands</h1>
 <span class="msep"></span>
 <input id="q" placeholder="filter &mdash; id, SMILES, warhead&hellip;" size="26">
 <select id="f"><option value="">all families</option>{opts}</select>
 <span class="mhint">{len(t):,} molecules &middot; {html.escape(score_label)}
   &middot; {html.escape(str(a.sort))}{f" @ {a.cutoff}" if a.sort == "fraction_above" else ""}</span>
 <span class="msep"></span>
 <a class="mbtn lnk" href="modes.html" title="the same ranking, one row per mode">per-mode view &#8599;</a>
 <a class="mbtn lnk" href="pipeline.html" title="how a molecule becomes a row">how this works &#8599;</a>
 <button class="mbtn" onclick="toggleTheme()">dark</button>
</div>
{nav}
<main><div id="left">
<p class="lede">One row per <b>molecule</b> &mdash; the Ranking page is per mode
({len(d):,} of them). <b>score</b> is the ordering column and <b>above cut</b> /
<b>modes</b> are its numerator and denominator. <b>best</b> is the molecule's
strongest single mode and the only depth-immune column here: every other one
moves with how long you docked, because the mode count never saturates (D0092).
<b>% supported</b> is the share of a molecule's modes H-bonding to Ser114 or
Ser115, the two residues flanking the catalytic sulfur &mdash; a bounded
multiplicative bonus (max +15%) that can never promote a pose which cannot react.
None of this is evidence that anything binds: it orders reachability of attack
geometry, and <code>rank_validated</code> is False.</p>
<div class="tw"><table><thead><tr>{head}</tr></thead><tbody>
{chr(10).join(body)}
</tbody></table></div>
</div>
<div id="right">
  <div id="vhead"><b id="vname">select a molecule</b>
    <div id="vmeta">its best-scoring mode, in the receptor it was docked into</div></div>
  <div id="gl"></div>
  <table id="vstats"><tbody></tbody></table>
  <img id="vthumb" class="thumb" alt="">
  <div id="vfoot">click a row on the left</div>
</div>
</main>
<script src="3Dmol-min.js"></script><script>const STATS={stats_json};const CUTS={cuts_json};</script>
<script>{_JS}</script></body></html>"""

    # ASSETS BESIDE THE PAGE, not inlined. modes.html embeds the whole 3Dmol
    # payload and is 57 MB for it; served from the same directory the library and
    # the receptor are cached by the browser and shared between pages, and the
    # page stays small enough to open. Same cached copy either way, so the two
    # builders cannot ship different versions.
    # THE PAGE IS WRITTEN WHERE ITS DATA CAME FROM. `--topic` selected which
    # ranking to READ and this wrote to `run.topic` regardless, so building a
    # page for any topic other than the current one silently overwrote the
    # CURRENT run's ligands.html with another run's contents -- under the
    # current run's title, on every server serving it. That is the
    # half-moved-topic defect (`how_this_project_breaks` #25) in a new place:
    # `reports_dir` has taken an optional topic all along and this caller did
    # not pass it.
    out = rp.reports_dir(topic)
    out.mkdir(parents=True, exist_ok=True)
    cache = REPO / "scripts" / ".cache_3dmol-min.js"
    if cache.is_file():
        tgt = out / "3Dmol-min.js"
        if not tgt.is_file() or tgt.stat().st_size != cache.stat().st_size:
            tgt.write_text(cache.read_text())
            log.info("wrote %s (%d KB)", tgt.name, tgt.stat().st_size // 1024)
    else:
        log.warning("no cached 3Dmol at %s — the viewer pane will say so rather "
                    "than being blank", cache.name)
    recsrc = rp.receptor_prep()
    tgt = out / "receptor.pdb"
    if not tgt.is_file() or tgt.stat().st_size != recsrc.stat().st_size:
        tgt.write_text(recsrc.read_text())
        log.info("wrote %s (%d KB)", tgt.name, tgt.stat().st_size // 1024)

    dest = rp.reports_dir(topic) / "ligands.html"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(page, encoding="utf-8")
    print(dest)
    print(f"  {len(t):,} ligands on {score_col} by {a.sort} · top: "
          f"{t.iloc[0].ident} (score {t.iloc[0].ligand_engagement:.3f}, "
          f"best {t.iloc[0].best_mode_engagement:.3f}, "
          f"{int(t.iloc[0].n_modes)} modes)")


if __name__ == "__main__":
    main()
