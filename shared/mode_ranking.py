"""The ranking view: every molecule and every mode, before anything is simulated.

TWO VIEWS, AND THIS IS THE FIRST ONE. `combined.html` shows the sweep and the
100 ns results -- the subset that was simulated. This shows what the screen
*scored*: every molecule, every mode, ranked, with its pose. In the pipeline's
real order this comes FIRST -- you read the ranked list, look at the poses, and
then choose what goes to the sweep. It was built retrospectively, after #53 found
that the sweep took mode 0 for 242 of 242 molecules while the ranking is per
mode, and that gap was invisible precisely because no view like this existed.

SCOPE IS ONE CONTROL, NOT TWO. A dropdown picks a warhead class, every class
ranked within itself, or one global order. Two separate toggles let a reader
combine "global" with a class filter, a combination with no meaning.

Within-class is the default. The SN2 angular criterion is far stricter than the
perpendicular one (#47), so a global order compares scores computed under
different bars. It is offered because "where does this sit overall" is a real
question, and refusing to answer it does not remove the bias -- the option names
it instead.

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
#: Surface shell: residues with a heavy atom within this of Cys113's SG. A
#: whole-protein VDW mesh is the expensive call in a 3Dmol viewer and the far
#: side of the protein is not the subject.
SURF_SHELL_A = 12.0


def pocket_residues() -> list[int]:
    """Residue numbers forming the pocket wall, from the receptor itself."""
    import numpy as np
    if not RECEPTOR.is_file():
        return []
    sg, res = None, {}
    for ln in RECEPTOR.read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        ri, nm = ln[22:26].strip(), ln[12:16].strip()
        try:
            xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
        res.setdefault(ri, []).append(xyz)
        if ri == str(CYS_RESI) and nm == "SG":
            sg = xyz
    if sg is None:
        return []
    return sorted(int(ri) for ri, xs in res.items()
                  if ri.lstrip("-").isdigit()
                  and min(float(np.linalg.norm(x - sg)) for x in xs) <= SURF_SHELL_A)


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
select#scope{font:600 11px var(--sans);padding:3px 26px 3px 10px;border-radius:99px;
 border:1px solid var(--rule);background:var(--paper);color:var(--ink);cursor:pointer;
 appearance:none;background-image:linear-gradient(45deg,transparent 50%,var(--muted) 50%),
 linear-gradient(135deg,var(--muted) 50%,transparent 50%);
 background-position:calc(100% - 14px) 52%,calc(100% - 9px) 52%;
 background-size:5px 5px,5px 5px;background-repeat:no-repeat}
select#scope:focus{outline:2px solid var(--blue);outline-offset:1px}
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
#vstructwrap{border:1px solid var(--rule);border-radius:4px;background:#fff;
 padding:6px;margin-bottom:10px;display:flex;justify-content:center}
#vstruct{width:100%;max-width:420px;height:auto}
#glbox{position:relative;width:100%;height:440px;background:#eef1f6;
 border:1px solid var(--rule);border-radius:4px;overflow:hidden}
#glbox>div{position:absolute;inset:0}
#glbox canvas{position:absolute;top:0;left:0}
.vctl{display:flex;flex-wrap:wrap;gap:.4rem 1.2rem;padding:.6rem .1rem 0;font-size:12px}
.vctl label{display:flex;align-items:center;gap:.35rem;cursor:pointer;font-weight:600}
.note{font-size:12px;color:var(--muted);margin:.9rem 0 0;max-width:78ch}
#sibs h3{font:600 11px var(--sans);letter-spacing:.05em;text-transform:uppercase;
 color:var(--muted);margin:1.3rem 0 .4rem}
table.sib{border-collapse:collapse;width:100%;font-size:12.5px}
table.sib th,table.sib td{padding:.34rem .6rem;border-bottom:1px solid var(--rule);
 text-align:right;white-space:nowrap}
table.sib th{font:600 10px var(--sans);color:var(--muted);text-transform:uppercase;
 letter-spacing:.05em;border-bottom:2px solid var(--rule)}
table.sib th:first-child,table.sib td:first-child,
table.sib th:last-child,table.sib td:last-child{text-align:left}
table.sib td{font-family:var(--mono)}
tr.sibrow{cursor:pointer}
tr.sibrow:hover{background:var(--blue-pale)}
tr.sibrow.cur{background:var(--blue-pale);font-weight:700}
input.mchk{margin:0 .45rem 0 0;vertical-align:-1px;cursor:pointer}
i.sw{width:11px;height:11px;border-radius:2px;display:inline-block;margin-right:.45rem;
 vertical-align:-1px}
.win{color:var(--good);font-weight:700}
.warnbox{border-left:3px solid var(--warn);background:#fdf8ea;padding:.6rem .9rem;
 font-size:12px;margin:0 0 12px;border-radius:0 3px 3px 0}
:root[data-theme="dark"] .warnbox{background:#241f12}
:root[data-theme="dark"] .thumb{background:#fff}
a{color:var(--blue)}
</style></head><body>
<div id="topbar">
 <h1 title="Pick a mode on the left; its pose and scores load on the right.">__TITLE__ — ranking</h1>
 <span class="msep"></span>
 <select id="scope" onchange="setScope(this.value)"
   title="rank within one warhead class, within every class, or across all of them"></select>
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
    <div id="vstructwrap"><img id="vstruct" alt="2D structure"></div>
    <div id="glbox"><div id="gl"></div></div>
    <div class="vctl">
     <label><input type="checkbox" id="c-surf" checked> pocket surface</label>
    </div>
    <div id="sibs"></div>
    <p class="note" id="vnote"></p>
   </div>
  </div>
 </div>
</main>
<script type="text/plain" id="recpdb">__RECEPTOR__</script>
<script>
const ROWS = __ROWS__;
const MODE_COLS = [0x0072ce, 0x7b5ea7, 0xc2703d, 0x0f7a54, 0xb3261e, 0x8a6d1f];
const MODE_CSS  = ['#0072ce','#7b5ea7','#c2703d','#0f7a54','#b3261e','#8a6d1f'];
const POCKET = __POCKET__;
// CARBONS CARRY THE MODE COLOUR; EVERY OTHER ELEMENT KEEPS ITS CONVENTIONAL ONE.
// Colouring a whole molecule by mode hides its chemistry -- the sulfur, the
// halogen and the oxygens are what a chemist reads a pose by, and they must look
// the same in every mode so the only thing that changes is the carbon skeleton.
function carbonScheme(col){
  const M = lib();
  return {prop:'elem',
          map: Object.assign({}, (M.elementColors||{}).defaultColors||{}, {C: col})};
}
// SCOPE is one control: a warhead class name, '*' for every class ranked within
// itself, or '__global__' for one order across all of them. Two orthogonal
// toggles let a reader combine "global" with a class filter, which is a
// combination with no meaning.
let SCOPE = '*', SEL = null, V = null, SURF = null;
// SHOWN holds the modes currently drawn for the selected molecule, so several
// alternatives can be compared and any of them switched off again. SEL stays the
// PRIMARY -- the one the facts panel describes -- because a panel of numbers has
// to be about one mode, and "several are visible" is a different question from
// "which one am I reading".
let SHOWN = new Set(), PDBCACHE = {};

function lib(){ return window.$3Dmol || window['3Dmol']; }
function fmt(x, d){ return (x === null || x === undefined) ? '—' : (+x).toFixed(d); }

function isGlobal(){ return SCOPE === '__global__'; }

function visible(){
  let r = ROWS.slice();
  if (SCOPE !== '*' && !isGlobal()) r = r.filter(x => x.c === SCOPE);
  if (isGlobal()) r.sort((a,b) => (a.gr === null ? 1e9 : a.gr) - (b.gr === null ? 1e9 : b.gr));
  else r.sort((a,b) => a.c.localeCompare(b.c) || a.cr - b.cr);
  return r;
}

function buildScope(){
  const n = {};
  ROWS.forEach(x => { n[x.c] = (n[x.c] || 0) + 1; });
  const opts = ['<optgroup label="ranked within its own class">',
    '<option value="*">all classes</option>'];
  Object.keys(n).sort().forEach(c =>
    opts.push('<option value="' + c + '">' + c + ' (' + n[c].toLocaleString() + ')</option>'));
  opts.push('</optgroup><optgroup label="across classes">',
    '<option value="__global__">global — biased (#47)</option></optgroup>');
  const el = document.getElementById('scope');
  el.innerHTML = opts.join('');
  el.value = SCOPE;
}

function railHTML(){
  const r = visible(), out = [];
  let cls = null;
  for (const x of r){
    if (SCOPE === '*' && x.c !== cls){ cls = x.c;
      out.push('<div class="chd">' + cls + '</div>'); }
    const rank = isGlobal() ? (x.gr === null ? '—' : x.gr) : x.cr;
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
    (isGlobal() ? ' · one order across classes scored under different bars (#47)'
                : ' · rank is within the warhead class');
  document.getElementById('rail').innerHTML = out.join('');
}

function setScope(v){ SCOPE = v; railHTML(); }

function toggleMode(m){
  // The primary cannot be hidden -- the facts panel is describing it, and a
  // panel of numbers with nothing drawn beside it reads as a rendering failure.
  const cur = SEL ? ROWS.find(r => r.i === SEL) : null;
  if (cur && m === cur.m) return;
  if (SHOWN.has(m)) SHOWN.delete(m); else SHOWN.add(m);
  if (SEL) pick(SEL);
}
function toggleTheme(){
  const d = document.documentElement.getAttribute('data-theme') === 'dark';
  document.documentElement.setAttribute('data-theme', d ? 'light' : 'dark'); }

async function pick(id){
  const x = ROWS.find(r => r.i === id); if (!x) return;
  const prev = SEL ? ROWS.find(r => r.i === SEL) : null;
  if (!prev || prev.p !== x.p) SHOWN = new Set();   // new molecule, new selection
  SHOWN.add(x.m);                                    // the primary is always drawn
  SEL = id; railHTML();
  document.getElementById('vempty').style.display = 'none';
  document.getElementById('vfull').style.display = '';
  document.getElementById('vname').textContent = x.i;
  // The same depiction the rail uses, at panel size. It is an SVG, so one file
  // serves both; drawing a second at a larger size would be a second answer to
  // "what does this molecule look like".
  document.getElementById('vstruct').src = 'mode_thumbs/' + x.p + '.svg';
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

  // EVERY MODE OF THIS MOLECULE, and how they were ranked against each other.
  // The rail orders modes across the whole library; this is the comparison the
  // pipeline claims to make -- a molecule's modes competing as separate rows --
  // and it is the one place a reader can see whether the mode that was simulated
  // is the one that scored best.
  const sibs = ROWS.filter(r => r.p === x.p).sort((a,b) => a.m - b.m);
  const best = sibs.reduce((a,b) => (b.cr < a.cr ? b : a), sibs[0]);
  document.getElementById('sibs').innerHTML =
    '<h3>modes of ' + x.p + ' — ' + sibs.length + ', ranked against each other</h3>' +
    '<p class="note" style="margin:.2rem 0 .5rem">Tick to draw a mode; click the '
    + 'row to read it. Several can be shown at once.</p>'
    + '<table class="sib"><thead><tr><th>show</th><th>class rank</th><th>poses</th>' +
    '<th>viable</th><th>enrichment</th><th>conditional_eb</th><th>spread</th>' +
    '<th>coherence</th><th>simulated</th></tr></thead><tbody>' +
    sibs.map(function(m){
      const col = MODE_CSS[(m.m >= 0 ? m.m : 0) % MODE_CSS.length];
      const badge = m.s === 'md' ? '100 ns' : m.s === 'swept' ? 'swept'
                  : m.s === 'failed' ? 'sweep failed' : 'never';
      return '<tr class="sibrow' + (m.i === x.i ? ' cur' : '') + '"'
        + ' onclick="pick(\'' + m.i + '\')">'
        + '<td onclick="event.stopPropagation();toggleMode(' + m.m + ')">'
        + '<input type="checkbox" class="mchk"' + (SHOWN.has(m.m) ? ' checked' : '')
        + ' onclick="event.stopPropagation();toggleMode(' + m.m + ')">'
        + '<i class="sw" style="background:' + col + '"></i>m' + m.m + '</td>'
        + '<td' + (m.i === best.i ? ' class="win"' : '') + '>' + m.cr + '</td>'
        + '<td>' + (m.n === null ? '—' : m.n) + '</td>'
        + '<td>' + (m.vf === null ? '—' : (m.vf*100).toFixed(1) + '%') + '</td>'
        + '<td>' + fmt(m.en, 2) + '</td><td>' + fmt(m.eb, 3) + '</td>'
        + '<td>' + fmt(m.sp, 2) + '</td><td>' + fmt(m.dc, 3) + '</td>'
        + '<td><span class="tag t-' + m.s + '">' + badge + '</span></td></tr>';
    }).join('') + '</tbody></table>' +
    (sibs.length > 1 && best.s === 'none'
      ? '<p class="note"><strong>The best-ranked mode of this molecule was never '
        + 'simulated.</strong> m' + best.m + ' ranks ' + best.cr + ' in '
        + best.c + '; the sweep took mode 0 (#53).</p>'
      : '');

  document.getElementById('vnote').innerHTML =
    'The pose is this mode\'s <strong>medoid</strong> — the pose most central to '
    + 'the mode among its best-anchored quartile, not its lowest-energy member. The '
    + 'individual poses were not persisted '
    + '(<a href="https://github.com/hallettmiket/inhibition/issues/44">#44</a>).';

  try {
    if (!PDBCACHE[x.p]){
      const res = await fetch('mode_poses/' + x.p + '.pdb');
      if (!res.ok) throw new Error(res.status);
      PDBCACHE[x.p] = await res.text();
    }
    draw(PDBCACHE[x.p], x);
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
  // CYS113 IN CONVENTIONAL ELEMENT COLOURS, so its sulfur reads as sulfur --
  // it is the atom the whole screen is aimed at. SG additionally gets a sphere,
  // because at stick radius a single S is easy to lose against the cartoon.
  V.setStyle({resi:[__CYS__]},
             {stick:{radius:0.28, colorscheme:'default'},
              cartoon:{color:'#c3ccd8', opacity:0.5}});
  V.addStyle({resi:[__CYS__], atom:'SG'}, {sphere:{radius:0.62}});
  // Draw every mode the reader has ticked. The PRIMARY is thicker and fully
  // opaque so it stays identifiable in a stack of alternatives; the rest are
  // thinner and translucent, which is the difference between comparing poses and
  // producing one unreadable object out of several.
  modes.forEach(function(mo, i){
    if (!SHOWN.has(mo)){ V.setStyle({model: i+1}, {}); return; }
    const primary = (mo === x.m);
    V.setStyle({model: i+1}, {stick:{
      radius: primary ? 0.22 : 0.14,
      opacity: primary ? 1 : 0.6,
      colorscheme: carbonScheme(MODE_COLS[(mo >= 0 ? mo : 0) % MODE_COLS.length])}});
  });
  if (document.getElementById('c-surf').checked){
    // The pocket shell only, and never over Cys113 or a ligand: a mesh drawn on
    // top of those hides the two things the panel exists to show. `and` rather
    // than a bare `not`, so the selection cannot leak onto the pose models.
    SURF = V.addSurface(M.SurfaceType.VDW, {opacity:0.62, color:'#b9c7db'},
      {and:[{model:0}, {resi:POCKET}, {not:{resi:[__CYS__]}}]});
  }
  const sel = modes.indexOf(x.m);
  V.zoomTo(sel >= 0 ? {model: sel+1} : {resn:'MOL'});
  V.zoom(0.5); V.resize();
  V.render();                       // 3Dmol draws NOTHING without this.
}

document.getElementById('c-surf').addEventListener('change', function(){ if (SEL) pick(SEL); });
buildScope();
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
            .replace("__POCKET__", json.dumps(pocket_residues()))
            .replace("__TITLE__", html.escape(title)))
