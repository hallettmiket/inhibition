"""The ranking view: every molecule and every mode, before anything is simulated.

TWO VIEWS, AND THIS IS THE FIRST ONE. `combined.html` shows the sweep and the
100 ns results -- the subset that was simulated. This shows what the screen
*scored*: every molecule, every mode, ranked, with its pose. In the pipeline's
real order this comes FIRST -- you read the ranked list, look at the poses, and
then choose what goes to the sweep. It was built retrospectively, after #53 found
that the sweep took mode 0 for 242 of 242 molecules while the ranking is per
mode, and that gap was invisible precisely because no view like this existed.

RANK WITHIN A WARHEAD CLASS IS THE DEFAULT; GLOBAL IS OFFERED AND FLAGGED. The
SN2 angular criterion is far stricter than the perpendicular one (#47), so a
global order compares scores computed under different bars. It is offered because
"where does this sit overall" is a real question, and refusing to answer it does
not remove the bias -- the toggle names it instead.

THE JOIN IS ON (parent_ident, mode). Never on `ident`: mode 0 is the bare ident
in the sweep table and `_m0` in the rank table, so a merge on the label silently
drops exactly the rows that were simulated (`shared/mode_key.py`).

Depictions and poses are FILES fetched on demand, not inlined. 8,096 rows of
base64 is tens of megabytes; the results GUI can inline its 59 and this cannot.
"""

from __future__ import annotations

import glob
import html
import json
from pathlib import Path

import pandas as pd

from shared import mode_key as mk

B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
RECEPTOR = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb")
CYS_RESI = 113


def _latest(pattern: str) -> Path | None:
    fs = sorted(glob.glob(str(B / pattern)))
    return Path(fs[-1]) if fs else None


def gather() -> pd.DataFrame:
    """One row per mode: rank, docking-derived scores, and what was simulated."""
    frames = []
    for tier, score in (("T4", "conditional_eb"), ("T3", "enrichment_conditional")):
        f = _latest(f"rank_v2/rank_v2_{tier}_{score}_*.csv")
        if f is None:
            continue
        d = pd.read_csv(f)
        if "mode" in d.columns:
            d = d[d["mode"].notna()]
        d["tier"] = d.get("tier", tier)
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    r = pd.concat(frames, ignore_index=True)

    sf = sorted(glob.glob(str(B / "attack_sweep/attack_sweep_*.csv")))
    sweep = (pd.concat([pd.read_csv(x) for x in sf], ignore_index=True)
             .drop_duplicates("ident", keep="last") if sf else pd.DataFrame())
    if not sweep.empty:
        keep = [c for c in ("ident", "parent_ident", "mode", "frac_attack_ready",
                            "n_visits", "status") if c in sweep.columns]
        # ATTEMPTED is not SUCCEEDED. A row exists for every mode sent;
        # frac_attack_ready is null when the run failed. Counting only the
        # successful ones as "swept" would report a mode that was tried and
        # crashed as one nobody ever chose.
        sweep = sweep.assign(_sent=True)
        # bare_is_mode_zero: these rows predate #53, when the sweep wrote the
        # bare ident for mode 0. Stated, not assumed.
        r = mk.join(r, sweep[keep + ["_sent"]].rename(
            columns={"status": "sweep_status"}),
            right_bare_is_mode_zero=True, suffixes=("", "_sw"))

    md_ids: set[str] = set()
    for f in glob.glob(str(B / "md_residence/*.csv")):
        try:
            d = pd.read_csv(f)
        except Exception:                                  # noqa: BLE001
            continue
        if "ident" not in d.columns or "production_ps" not in d.columns:
            continue
        d = d[(d.production_ps >= 50000)
              & d.status.astype(str).str.startswith("ok")]
        md_ids |= set(d.ident.astype(str))

    r["sent"] = (r["_sent"].fillna(False).astype(bool)
                 if "_sent" in r.columns else False)
    r["swept"] = (r["frac_attack_ready"].notna()
                  if "frac_attack_ready" in r.columns else False)
    # A mode gets the 100 ns badge only if its molecule ran AND this mode is the
    # one that was sent. The MD rows do not record their mode (#36), so anything
    # looser would badge a mode that never moved.
    r["ran_md"] = r.parent_ident.isin(md_ids) & r["sent"]

    # GLOBAL RANK IS COMPUTED HERE AND LABELLED BIASED WHEREVER IT IS SHOWN.
    # conditional_eb is not comparable across warhead classes (#47); this exists
    # so the page can answer "where does this sit overall" while saying so.
    if "conditional_eb" in r.columns:
        r["global_rank"] = r["conditional_eb"].rank(
            ascending=False, method="min", na_option="bottom")
    return r


def idents(r: pd.DataFrame) -> set[str]:
    """The molecules the view will ask for assets for."""
    return set(r.parent_ident.astype(str)) if not r.empty else set()


def _rows_json(r: pd.DataFrame) -> str:
    out = []
    for _, x in r.iterrows():
        if pd.isna(x.get("class_rank")):
            continue
        st = ("md" if x.ran_md else "swept" if x.swept
              else "failed" if x.sent else "none")

        def num(k, nd=None):
            v = x.get(k)
            if pd.isna(v):
                return None
            return round(float(v), nd) if nd is not None else int(v)

        out.append({
            "i": str(x.ident), "p": str(x.parent_ident), "m": int(x["mode"]),
            "c": str(x.warhead_class), "cr": int(x.class_rank),
            "gr": num("global_rank"), "n": num("n_poses_mode"),
            "np": num("n_poses"), "vf": num("viable_fraction", 4),
            "eb": num("conditional_eb", 3), "en": num("enrichment", 2),
            "sp": num("spread_a", 2), "dc": num("dir_coherence", 3),
            "fa": num("frac_attack_ready", 4), "s": st,
        })
    return json.dumps(out, separators=(",", ":"))


_TPL = r"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__ — ranking</title>
<script>__THREE__</script>
<style>
/* The same palette and the same shell as the results GUI, so the two views read
   as one instrument rather than two pages that happen to link to each other. */
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
.mbtn{font:600 11px var(--sans);padding:3px 10px;border-radius:99px;cursor:pointer;
 border:1px solid var(--rule);background:var(--paper);color:var(--ink)}
.mbtn.on{background:var(--navy);color:#fff;border-color:var(--navy)}
.mbtn.lnk{text-decoration:none;color:var(--blue)}
.msep{width:1px;height:16px;background:var(--rule);flex:none}
.mhint{font-size:11px;color:var(--muted);margin-left:4px}
main{flex:1;display:grid;grid-template-columns:376px 1fr;min-height:0}
@media(max-width:880px){main{grid-template-columns:1fr;grid-template-rows:250px 1fr}}
#rail{overflow-y:auto;border-right:1px solid var(--rule);background:var(--rail)}
.chd{font-family:var(--mono);font-size:.6rem;letter-spacing:.14em;
 text-transform:uppercase;color:var(--blue);font-weight:600;padding:10px 14px 6px;
 background:var(--raise);border-bottom:1px solid var(--rule);position:sticky;top:0;z-index:1}
.row{display:grid;grid-template-columns:30px 46px 1fr;gap:8px;align-items:start;
 width:100%;text-align:left;font:inherit;color:inherit;background:none;cursor:pointer;
 padding:9px 14px 8px;border:0;border-bottom:1px solid var(--rule)}
.row:hover{background:var(--blue-pale)}
.row.on{background:var(--blue-pale);box-shadow:inset 3px 0 0 var(--blue)}
.rk{font-family:var(--mono);font-size:11px;font-weight:600;color:var(--muted);padding-top:3px}
.thumb{width:46px;height:32px;object-fit:contain;background:#fff;
 border:1px solid var(--rule);border-radius:3px;display:block}
.body{min-width:0;display:flex;flex-direction:column;gap:3px}
.l1{display:flex;align-items:baseline;justify-content:space-between;gap:8px}
.mid-id{font-family:var(--mono);font-size:12.5px;overflow:hidden;
 text-overflow:ellipsis;white-space:nowrap}
.eng{font-family:var(--mono);font-size:11.5px;font-weight:600;color:var(--navy);flex:none}
.l2{display:flex;align-items:center;gap:6px;min-width:0}
.wc{font-size:10.5px;color:var(--blue);white-space:nowrap;overflow:hidden;
 text-overflow:ellipsis;max-width:12ch}
.meta{font-size:10.5px;color:var(--muted);white-space:nowrap;overflow:hidden;
 text-overflow:ellipsis;flex:1}
.tag{font-size:9px;letter-spacing:.06em;text-transform:uppercase;font-weight:700;
 padding:1px 6px;border-radius:99px;flex:none}
.t-md{background:#e6f4ee;color:var(--good)}
.t-swept{background:#e8f1fb;color:var(--navy)}
.t-failed{background:#faf3e0;color:var(--warn)}
.t-none{background:var(--raise);color:var(--muted)}
.bar{height:3px;background:var(--rule);border-radius:2px;overflow:hidden;margin-top:2px}
.bar i{display:block;height:100%;background:var(--blue)}
#viewer{min-width:0;min-height:0;display:flex;flex-direction:column;background:var(--paper)}
#vhead{padding:9px 18px;border-bottom:1px solid var(--rule);background:var(--raise);
 display:flex;justify-content:space-between;align-items:center;gap:12px}
#vname{font-family:var(--mono);font-size:12.5px;font-weight:600;color:var(--navy)}
#vbody{flex:1;min-height:0;overflow-y:auto;padding:16px 18px 30px}
.facts{display:grid;grid-template-columns:repeat(auto-fit,minmax(128px,1fr));
 gap:10px 18px;margin-bottom:14px}
.fact b{display:block;font-family:var(--mono);font-size:1.05rem}
.fact span{font-size:11px;color:var(--muted)}
#glbox{position:relative;width:100%;height:440px;background:#eef1f6;
 border:1px solid var(--rule);border-radius:4px;overflow:hidden}
#glbox>div{position:absolute;inset:0}
#glbox canvas{position:absolute;top:0;left:0}
.vctl{display:flex;flex-wrap:wrap;gap:.4rem 1.2rem;padding:.6rem .1rem 0;font-size:12px}
.vctl label{display:flex;align-items:center;gap:.35rem;cursor:pointer;font-weight:600}
.note{font-size:12px;color:var(--muted);margin:.9rem 0 0;max-width:78ch}
.warnbox{border-left:3px solid var(--warn);background:#fdf8ea;padding:.6rem .9rem;
 font-size:12px;margin:0 0 12px;border-radius:0 3px 3px 0}
:root[data-theme="dark"] .warnbox{background:#241f12}
:root[data-theme="dark"] .thumb{background:#fff}
a{color:var(--blue)}
</style></head><body>
<div id="topbar">
 <h1 title="Pick a mode on the left; its pose and scores load on the right.">__TITLE__ — ranking</h1>
 <span class="msep"></span>
 <button id="b-class" class="mbtn on" onclick="setMode('class')"
   title="rank within a warhead class — the only comparison the criterion supports">by warhead class</button>
 <button id="b-global" class="mbtn" onclick="setMode('global')"
   title="one order across all classes — biased, see the hint">global</button>
 <span class="msep"></span>
 <button id="b-all" class="mbtn on" onclick="setFilter('all')">all modes</button>
 <button id="b-un" class="mbtn" onclick="setFilter('unsimulated')"
   title="modes the screen scored and never simulated">never simulated</button>
 <span class="mhint" id="mhint"></span>
 <span class="msep"></span>
 <a class="mbtn lnk" href="combined.html" title="the sweep and 100 ns MD results">results &#8599;</a>
 <a class="mbtn lnk" href="pipeline.html" title="how a molecule becomes a row">how this works &#8599;</a>
 <button class="mbtn" onclick="toggleTheme()">dark</button>
</div>
<main>
 <div id="rail"></div>
 <div id="viewer">
  <div id="vhead"><span id="vname">select a mode</span>
   <span class="mhint">medoid pose, in the receptor it was docked into</span></div>
  <div id="vbody">
   <div id="vempty" class="note">This is the <strong>ranking</strong> view: every
   molecule and every mode the screen scored, simulated or not. In the pipeline's
   real order it comes first — you read the ranked list and look at the poses, then
   choose what goes to the 10&nbsp;ns sweep. It was built after the fact, which is
   how <a href="https://github.com/hallettmiket/inhibition/issues/53">#53</a> went
   unnoticed: nothing in the project had ever shown the per-mode ranking it
   computes.</div>
   <div id="vfull" style="display:none">
    <div id="gwarn" class="warnbox" style="display:none"></div>
    <div class="facts" id="vfacts"></div>
    <div id="glbox"><div id="gl"></div></div>
    <div class="vctl">
     <label><input type="checkbox" id="c-surf" checked> pocket surface</label>
     <label><input type="checkbox" id="c-other"> other modes of this molecule</label>
    </div>
    <p class="note" id="vnote"></p>
   </div>
  </div>
 </div>
</main>
<script type="text/plain" id="recpdb">__RECEPTOR__</script>
<script>
const ROWS = __ROWS__;
const MODE_COLS = ['#0072ce','#7b5ea7','#c2703d','#0f7a54','#b3261e','#8a6d1f'];
let RANKMODE = 'class', FILTER = 'all', SEL = null, V = null, SURF = null;

function lib(){ return window.$3Dmol || window['3Dmol']; }
function fmt(x, d){ return (x === null || x === undefined) ? '—' : (+x).toFixed(d); }

function visible(){
  let r = (FILTER === 'unsimulated') ? ROWS.filter(x => x.s === 'none') : ROWS.slice();
  if (RANKMODE === 'class') r.sort((a,b) => a.c.localeCompare(b.c) || a.cr - b.cr);
  else r.sort((a,b) => (a.gr === null ? 1e9 : a.gr) - (b.gr === null ? 1e9 : b.gr));
  return r;
}

function railHTML(){
  const r = visible(), out = [];
  let cls = null;
  for (const x of r){
    if (RANKMODE === 'class' && x.c !== cls){ cls = x.c;
      out.push('<div class="chd">' + cls + '</div>'); }
    const rank = (RANKMODE === 'class') ? x.cr : (x.gr === null ? '—' : x.gr);
    const pct = x.vf === null ? 0 : Math.round(x.vf * 100);
    const badge = x.s === 'md' ? '100 ns' : x.s === 'swept' ? 'swept'
                : x.s === 'failed' ? 'failed' : 'not run';
    out.push(
      '<button class="row' + (SEL === x.i ? ' on' : '') + '" onclick="pick(\'' + x.i + '\')">' +
      '<span class="rk">' + rank + '</span>' +
      '<img class="thumb" loading="lazy" alt="" src="mode_thumbs/' + x.p + '.svg">' +
      '<span class="body"><span class="l1">' +
      '<span class="mid-id">' + x.i + '</span>' +
      '<span class="eng">' + fmt(x.eb, 2) + '</span></span>' +
      '<span class="l2"><span class="wc">' + x.c + '</span>' +
      '<span class="meta">' + (x.n === null ? '—' : x.n) + ' poses · ' + pct + '% viable</span>' +
      '<span class="tag t-' + x.s + '">' + badge + '</span></span>' +
      '<span class="bar"><i style="width:' + pct + '%"></i></span></span></button>');
  }
  document.getElementById('mhint').textContent =
    r.length.toLocaleString() + ' modes' +
    (RANKMODE === 'global' ? ' · global order compares classes scored under different bars (#47)' : '');
  document.getElementById('rail').innerHTML = out.join('');
}

function setMode(m){ RANKMODE = m;
  document.getElementById('b-class').classList.toggle('on', m === 'class');
  document.getElementById('b-global').classList.toggle('on', m === 'global');
  railHTML(); }
function setFilter(f){ FILTER = f;
  document.getElementById('b-all').classList.toggle('on', f === 'all');
  document.getElementById('b-un').classList.toggle('on', f === 'unsimulated');
  railHTML(); }
function toggleTheme(){
  const d = document.documentElement.getAttribute('data-theme') === 'dark';
  document.documentElement.setAttribute('data-theme', d ? 'light' : 'dark'); }

async function pick(id){
  const x = ROWS.find(r => r.i === id); if (!x) return;
  SEL = id; railHTML();
  document.getElementById('vempty').style.display = 'none';
  document.getElementById('vfull').style.display = '';
  document.getElementById('vname').textContent = x.i;
  document.getElementById('vfacts').innerHTML = [
    ['class rank', x.cr],
    ['global rank', x.gr === null ? '—' : x.gr],
    ['poses in mode', x.n === null ? '—' : x.n + ' of ' + (x.np === null ? '?' : x.np)],
    ['viable fraction', x.vf === null ? '—' : (x.vf*100).toFixed(1) + '%'],
    ['enrichment', fmt(x.en, 2)],
    ['conditional_eb', fmt(x.eb, 3)],
    ['spread', fmt(x.sp, 2) + ' Å'],
    ['direction coherence', fmt(x.dc, 3)],
    ['10 ns sweep', x.fa === null ? (x.s === 'none' ? 'never sent' : 'no score')
                                  : (x.fa*100).toFixed(1) + '% ready'],
  ].map(kv => '<div class="fact"><b>' + kv[1] + '</b><span>' + kv[0] + '</span></div>').join('');

  const g = document.getElementById('gwarn');
  if (x.s === 'none'){
    g.style.display = '';
    g.innerHTML = '<strong>Never simulated.</strong> This mode was scored and ranked, '
      + 'and no sweep or MD was ever run from it. Every number above is docking-derived.';
  } else { g.style.display = 'none'; }

  document.getElementById('vnote').innerHTML =
    'The pose is this mode\'s <strong>medoid</strong> — the pose most central to '
    + 'the mode among its best-anchored quartile, not its lowest-energy member. The '
    + 'individual poses were not persisted '
    + '(<a href="https://github.com/hallettmiket/inhibition/issues/44">#44</a>).';

  try {
    const res = await fetch('mode_poses/' + x.p + '.pdb');
    if (!res.ok) throw new Error(res.status);
    draw(await res.text(), x);
  } catch (e) {
    document.getElementById('gl').innerHTML =
      '<p class="note" style="padding:14px">no pose file for ' + x.p + '</p>';
  }
}

function draw(pdbTxt, x){
  const M = lib(); if (!M) return;
  if (!V) V = M.createViewer(document.getElementById('gl'), {backgroundColor:'#eef1f6'});
  V.clear(); SURF = null;
  V.addModel(document.getElementById('recpdb').textContent, 'pdb');
  // Each block carries its own mode in the MODEL record, so the model->mode
  // mapping is READ rather than counted. Counting positions is #53.
  const blocks = pdbTxt.split('ENDMDL').filter(b => b.indexOf('MODEL') >= 0);
  const modes = [];
  blocks.forEach(b => {
    const m = /MODEL\s+(-?\d+)/.exec(b);
    modes.push(m ? parseInt(m[1], 10) : -1);
    V.addModel(b.replace(/MODEL[^\n]*\n/, ''), 'pdb');
  });
  V.setStyle({}, {cartoon:{color:'#c3ccd8', opacity:0.5}});
  const CC = Object.assign({}, (M.elementColors||{}).defaultColors||{}, {C:0xb3261e});
  V.setStyle({resi:[__CYS__]}, {stick:{radius:0.26, colorscheme:{prop:'elem', map:CC}},
                               cartoon:{color:'#c3ccd8', opacity:0.5}});
  const showOther = document.getElementById('c-other').checked;
  modes.forEach(function(mo, i){
    const on = (mo === x.m) || showOther;
    V.setStyle({model: i+1}, on ? {stick:{
      radius: (mo === x.m) ? 0.22 : 0.13,
      color: MODE_COLS[(mo >= 0 ? mo : 0) % MODE_COLS.length],
      opacity: (mo === x.m) ? 1 : 0.45}} : {});
  });
  if (document.getElementById('c-surf').checked){
    SURF = V.addSurface(M.SurfaceType.VDW, {opacity:0.55, color:'#b9c7db'},
                        {model:0, not:{resi:[__CYS__]}});
  }
  const sel = modes.indexOf(x.m);
  V.zoomTo(sel >= 0 ? {model: sel+1} : {resn:'MOL'});
  V.zoom(0.5); V.resize();
  V.render();                       // 3Dmol draws NOTHING without this.
}

document.getElementById('c-surf').addEventListener('change', function(){ if (SEL) pick(SEL); });
document.getElementById('c-other').addEventListener('change', function(){ if (SEL) pick(SEL); });
railHTML();
</script>
</body></html>"""


def build(title: str, date_str: str, three: str = "") -> str:
    r = gather()
    if r.empty:
        return "<!doctype html><p>no rank tables found</p>"
    rec = ("\n".join(l for l in RECEPTOR.read_text().splitlines()
                     if l.startswith(("ATOM", "HETATM")))
           if RECEPTOR.is_file() else "")
    return (_TPL
            .replace("__ROWS__", _rows_json(r))
            .replace("__RECEPTOR__", rec)
            .replace("__THREE__", three)
            .replace("__CYS__", str(CYS_RESI))
            .replace("__TITLE__", html.escape(title)))
