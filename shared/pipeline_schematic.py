"""
Purpose: the pipeline schematic — how one molecule becomes one row in the catalogue.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-09
Input: nothing — the geometry is illustrative and generated here
Output: a self-contained HTML page, linked from the catalogue's top bar

@tt8804: *"make a schematic showing how our pipeline works ... show 500 poses,
doesn't have to be real, then poses splitting into 3 modes ... then the criteria
next to them, showing ranking and then sweep and md."*

THE POSES ARE DRAWN, NOT MEASURED, AND THE PAGE SAYS SO IN ITS OWN MASTHEAD. This
is the one place in the project where invented numbers are on screen, so the label
is not a footnote — a reader who mistakes this for data would take exactly the kind
of plausible-and-wrong reading `how_this_project_breaks.md` catalogues. Every
PARAMETER shown (500 runs, the 2.8–4.2 Å window, the 1.2 nm residence bar, 10 ns
and 100 ns) is real and cited; every COORDINATE is fabricated.

DETERMINISTIC BY CONSTRUCTION. `random.Random(SEED)` rather than an unseeded draw,
because `outputs.py` versions every write: an unseeded schematic would produce a
new file on every rebuild and the versioned tree would fill with identical-looking
pages that differ only in dot positions.
"""

from __future__ import annotations

import collections
import json
import math
import random
from pathlib import Path

SEED = 7
N_POSES = 500

#: A real screened molecule to illustrate with: 491 poses, exactly three modes.
#: Named here rather than discovered, because the page is a worked example and it
#: should show the SAME molecule every time it is rebuilt.
EXAMPLE = "t4_0e251ffccad1"
ALLPOSES = ("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
            "nac_v3_allposes")
RECEPTOR = ("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/"
            "3IKD_noligand.pdb")
THREEDMOL = "scripts/.cache_3dmol-min.js"


def _read_sdf(path):
    """Per-pose mode, centroid and raw block — without paying for RDKit.

    The mode tag is written by the screen itself, so the grouping shown on this
    page is the real one; only the single dot standing in for each pose is a
    summary of it.
    """
    out = []
    for blk in path.read_text(errors="ignore").split("$$$$"):
        lines = blk.splitlines()
        if len(lines) < 5:
            continue
        # FIND THE COUNTS LINE BY WHAT IT IS, not by where it sits. The title line
        # is blank in these files, so any leading-newline handling shifts a fixed
        # index onto an atom line, int() fails, and every pose is silently skipped
        # -- an empty panel from a file that parses perfectly well.
        ci = next((i for i, ln in enumerate(lines) if ln.rstrip().endswith("V2000")), -1)
        if ci < 0:
            continue
        try:
            n = int(lines[ci][:3])
        except ValueError:
            continue
        xs = ys = zs = 0.0
        k = 0
        for ln in lines[ci + 1:ci + 1 + n]:
            try:
                x, y, z = float(ln[0:10]), float(ln[10:20]), float(ln[20:30])
            except ValueError:
                continue
            sym = ln[31:34].strip()
            if sym == "H":
                continue
            xs += x; ys += y; zs += z; k += 1
        if not k:
            continue
        mode = "0"
        for i, ln in enumerate(lines):
            if ln.startswith(">") and "<mode>" in ln and i + 1 < len(lines):
                mode = lines[i + 1].strip()
                break
        out.append({"mode": mode, "c": [xs / k, ys / k, zs / k],
                    "blk": "\n".join(lines) + "\n$$$$\n"})
    return out

#: The three illustrative modes: share of the 500, centre, spread, and the story
#: each one tells. Sizes are deliberately uneven — a dominant mode plus two
#: minority modes is the common real shape, and a molecule promoted on a minority
#: mode is a different claim from one promoted on its dominant mode.
MODES = [
    {"key": "A", "n": 246, "cx": 214, "cy": 150, "sx": 34, "sy": 25,
     "col": "#0072ce", "d": 3.4, "ang": 168, "ar": 0.41, "note":
     "warhead on the sulfur, roughly in line — inside the window and near-linear"},
    {"key": "B", "n": 158, "cx": 292, "cy": 186, "sx": 27, "sy": 20,
     "col": "#7b5ea7", "d": 4.9, "ang": 141, "ar": 0.06, "note":
     "sits deeper in the pocket, warhead too far out to react"},
    {"key": "C", "n": 96, "cx": 176, "cy": 208, "sx": 22, "sy": 17,
     "col": "#c2703d", "d": 3.1, "ang": 96, "ar": 0.02, "note":
     "close enough, but approaching side-on — distance alone would pass it"},
]

SG = (268, 132)          # the Cys113 sulfur, in stage-1/2 coordinates


def _points(rng):
    """One cloud, reused by every stage so the reader is following the same 500."""
    pts = []
    for m in MODES:
        for _ in range(m["n"]):
            x = rng.gauss(m["cx"], m["sx"])
            y = rng.gauss(m["cy"], m["sy"])
            z = rng.random()                      # depth: drives size + opacity
            pts.append({"x": x, "y": y, "z": z, "m": m["key"], "col": m["col"]})
    rng.shuffle(pts)                              # so no mode paints on top
    return pts


def _dots(pts, colour=None, only=None, ox=0.0, oy=0.0, scale=1.0):
    out = []
    for p in pts:
        if only and p["m"] != only:
            continue
        r = (1.1 + 1.3 * p["z"]) * scale
        op = 0.30 + 0.60 * p["z"]
        c = colour or p["col"]
        x = ox + p["x"] * scale
        y = oy + p["y"] * scale
        out.append(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r:.2f}' "
                   f"fill='{c}' opacity='{op:.2f}'/>")
    return "".join(out)


def _pocket(idp: str) -> str:
    """A cut-away pocket with depth cues, so the cloud reads as sitting *in* it."""
    return f"""
<defs>
  <radialGradient id='cav{idp}' cx='45%' cy='40%' r='70%'>
    <stop offset='0%' stop-color='var(--cav-in)'/>
    <stop offset='100%' stop-color='var(--cav-out)'/>
  </radialGradient>
  <radialGradient id='sg{idp}' cx='35%' cy='32%' r='68%'>
    <stop offset='0%' stop-color='#ffe58a'/><stop offset='100%' stop-color='#c9971a'/>
  </radialGradient>
</defs>
<path d='M40 96 C70 40 190 22 268 34 C338 45 402 84 404 150
         C406 216 344 268 262 274 C176 280 78 250 52 196 C36 164 30 124 40 96 Z'
      fill='var(--prot)' stroke='var(--rule)' stroke-width='1'/>
<ellipse cx='232' cy='176' rx='128' ry='84' fill='url(#cav{idp})'
         stroke='var(--rule)' stroke-width='1'/>
<ellipse cx='232' cy='176' rx='96' ry='60' fill='none' stroke='var(--rule)'
         stroke-width='.7' stroke-dasharray='2 4' opacity='.65'/>
<circle cx='{SG[0]}' cy='{SG[1]}' r='11' fill='url(#sg{idp})'/>
<text x='{SG[0] + 16}' y='{SG[1] + 4}' class='lbl'>Cys113 S&gamma;</text>
"""


def _stage1(pts) -> str:
    return f"""<svg viewBox="0 0 430 300" class="dia" role="img"
 aria-label="500 docked poses scattered through the Pin1 pocket around Cys113">
{_pocket('1')}
{_dots(pts, colour='var(--dot)')}
<text x="14" y="292" class="cap">500 poses &middot; one molecule</text>
</svg>"""


def _stage2(pts) -> str:
    hulls = "".join(
        f"<ellipse cx='{m['cx']}' cy='{m['cy']}' rx='{m['sx'] * 2.1:.0f}' "
        f"ry='{m['sy'] * 2.1:.0f}' fill='{m['col']}' opacity='.10' "
        f"stroke='{m['col']}' stroke-width='1' stroke-dasharray='3 3'/>"
        for m in MODES)
    tags = "".join(
        f"<text x='{m['cx']}' y='{m['cy'] - m['sy'] * 2.1 - 6:.0f}' "
        f"class='mtag' fill='{m['col']}'>mode {m['key']}</text>" for m in MODES)
    return f"""<svg viewBox="0 0 430 300" class="dia" role="img"
 aria-label="The same 500 poses, coloured into three binding modes">
{_pocket('2')}
{hulls}{_dots(pts)}{tags}
<text x="14" y="292" class="cap">clustered on the reactive atom &amp; where the warhead points</text>
</svg>"""


def _mode_panel(pts, m) -> str:
    """One mode: its own share of the cloud, plus the representative it elects."""
    ox, oy, sc = -108.0, -76.0, 0.86
    cx = ox + m["cx"] * sc
    cy = oy + m["cy"] * sc
    ang = math.radians(180 - m["ang"])
    sgx = ox + SG[0] * sc
    sgy = oy + SG[1] * sc
    wx = cx + 34 * math.cos(ang)
    wy = cy - 34 * math.sin(ang)
    return f"""<svg viewBox="0 0 200 150" class="mini" role="img"
 aria-label="Mode {m['key']}: {m['n']} poses, representative geometry">
<ellipse cx='{ox + 232 * sc:.0f}' cy='{oy + 176 * sc:.0f}'
         rx='{128 * sc:.0f}' ry='{84 * sc:.0f}' fill='var(--cav-out)'
         stroke='var(--rule)' stroke-width='.8'/>
{_dots(pts, only=m['key'], ox=ox, oy=oy, scale=sc)}
<line x1='{cx:.1f}' y1='{cy:.1f}' x2='{sgx:.1f}' y2='{sgy:.1f}'
      stroke='{m['col']}' stroke-width='1.4' stroke-dasharray='3 2'/>
<line x1='{cx:.1f}' y1='{cy:.1f}' x2='{wx:.1f}' y2='{wy:.1f}'
      stroke='{m['col']}' stroke-width='2.4' stroke-linecap='round'/>
<circle cx='{cx:.1f}' cy='{cy:.1f}' r='4.2' fill='{m['col']}'/>
<circle cx='{sgx:.1f}' cy='{sgy:.1f}' r='7' fill='url(#sg2)'/>
<text x='8' y='16' class='mtag' fill='{m['col']}'>mode {m['key']}</text>
<text x='8' y='142' class='cap'>{m['n']} poses</text>
</svg>"""


def _real_poses() -> tuple[str, dict]:
    """The same story, on real output: 491 poses of one screened molecule.

    One dot per pose at its heavy-atom centroid, coloured by the mode the screen
    assigned it, inside the real 3IKD receptor. Dots rather than 491 sets of
    sticks because the page has to load: the cloud is the point, and one
    representative per mode is drawn in full beside it.
    """
    sdf = Path(ALLPOSES) / f"{EXAMPLE}.sdf"
    rec = Path(RECEPTOR)
    js = Path(__file__).resolve().parent.parent / THREEDMOL
    if not (sdf.is_file() and rec.is_file() and js.is_file()):
        return "", {}
    poses = _read_sdf(sdf)
    if not poses:
        return "", {}
    order = [m for m, _ in
             sorted(collections.Counter(p["mode"] for p in poses).items(),
                    key=lambda kv: -kv[1])]
    cols = {m: c for m, c in zip(order, ["#0072ce", "#7b5ea7", "#c2703d",
                                         "#0f7a54", "#b3261e"])}
    # EVERY POSE AS A STRUCTURE, NOT A DOT. A centroid says where a pose sat and
    # nothing about how it sat, which is the thing the modes actually differ in.
    # Grouped into one multi-molecule SDF per mode so each group can be styled and
    # hidden on its own -- and so 3Dmol reads the explicit bond block rather than
    # inferring bonds by distance, which across 491 overlapping poses would wire
    # neighbouring molecules to each other.
    by_mode: dict[str, list[str]] = {}
    for p in poses:
        by_mode.setdefault(p["mode"], []).append(p["blk"])
    sdf_by = {m: "".join(v) for m, v in by_mode.items()}
    reps = {}
    for p in poses:
        reps.setdefault(p["mode"], p["blk"])
    counts = collections.Counter(p["mode"] for p in poses)
    btns = "".join(
        f"<button class='p3b' data-m='{m}' onclick=\"pmode('{m}')\">"
        f"<span class='dotc' style='background:{cols[m]}'></span>mode {m}"
        f" <b>{counts[m]}</b></button>" for m in order)
    block = f"""
<div class="p3wrap">
 <div class="p3ctl">
  <button class="p3b on" data-m="all" onclick="pmode('all')">all {len(poses)} poses</button>
  {btns}
  <label class="p3l"><input id="p3rep" type="checkbox" checked> show one full pose</label>
 </div>
 <div class="p3box"><div id="p3"></div>
   <div id="p3wait" class="p3wait">building {len(poses)} poses&hellip;</div></div>
 <p class="p3cap">Real output for <code>{EXAMPLE}</code> — all {len(poses)} poses,
 {len(order)} modes, in 3IKD. Every pose is drawn as a structure, coloured by the
 mode the screen assigned. Tick <em>one full pose</em> to pick a single one out in
 sticks.</p>
</div>
<script>{js.read_text(errors='ignore')}</script>
<script type="text/plain" id="p3rec">{rec.read_text(errors='ignore')}</script>
<script type="text/plain" id="p3sdf">{json.dumps(sdf_by)}</script>
<script type="text/plain" id="p3reps">{json.dumps(reps)}</script>
<script>
(function(){{
  const M = window.$3Dmol || window['3Dmol'];
  const COLS = {json.dumps(cols)};
  const SDF = JSON.parse(document.getElementById('p3sdf').textContent);
  const REPS = JSON.parse(document.getElementById('p3reps').textContent);
  let v = null, cur = 'all', groups = {{}}, repModel = null;
  function style(){{
    Object.keys(groups).forEach(function(m){{
      const on = (cur === 'all' || cur === m);
      groups[m].forEach(function(mod){{
        mod.setStyle({{}}, on ? {{line:{{colorscheme: 'default', color: COLS[m],
                                       linewidth: cur === 'all' ? 1.0 : 1.6}}}} : {{}});
      }});
    }});
  }}
  function drawRep(){{
    if (repModel) {{ try {{ v.removeModel(repModel); }} catch(e) {{}} repModel = null; }}
    if (!document.getElementById('p3rep').checked) return;
    const m = cur === 'all' ? Object.keys(REPS)[0] : cur;
    if (!REPS[m]) return;
    repModel = v.addModel(REPS[m], 'sdf');
    repModel.setStyle({{}}, {{stick:{{radius:0.17, colorscheme:'yellowCarbon'}}}});
  }}
  function render(){{ style(); drawRep(); v.render(); }}
  window.pmode = function(m){{
    cur = m;
    document.querySelectorAll('.p3b').forEach(function(b){{
      b.classList.toggle('on', b.dataset.m === m); }});
    render();
  }};
  window.addEventListener('load', function(){{
    requestAnimationFrame(function(){{ requestAnimationFrame(function(){{
      v = M.createViewer(document.getElementById('p3'), {{backgroundColor:'#eef1f6'}});
      v.addModel(document.getElementById('p3rec').textContent, 'pdb');
      v.setStyle({{}}, {{cartoon:{{color:'#9fb0c4', opacity:0.72}}}});
      // addModels returns the models it made in most builds and a single model in
      // some; normalise rather than trusting one shape.
      Object.keys(SDF).forEach(function(m){{
        const r = v.addModels(SDF[m], 'sdf');
        groups[m] = Array.isArray(r) ? r : [r];
      }});
      render();
      v.zoomTo({{model: groups[Object.keys(groups)[0]][0]}});
      v.zoom(0.8); v.resize();
      const w = document.getElementById('p3wait'); if (w) w.style.display = 'none';
      document.getElementById('p3rep').addEventListener('change', render);
    }}); }});
  }});
}})();
</script>"""
    return block, {"n": len(poses), "modes": len(order), "counts": dict(counts)}


def build() -> str:
    rng = random.Random(SEED)
    pts = _points(rng)
    real_block, real = _real_poses()
    # Built outside the f-string: a dict literal inside an f-string expression is
    # read as a set of a set, which is unhashable and fails at build time.
    _c = real.get("counts") or {}
    counts_txt = ("" if not _c else
                  " &mdash; " + " / ".join(f"mode {k}: {v}" for k, v in _c.items()))

    panels = "".join(
        f"<div class='mcard'>{_mode_panel(pts, m)}"
        f"<table class='crit'>"
        f"<tr><th>d(C&rarr;S&gamma;)</th><td class='n'>{m['d']:.1f} &Aring;</td>"
        f"<td class='v {'ok' if 2.8 <= m['d'] <= 4.2 else 'no'}'>"
        f"{'in window' if 2.8 <= m['d'] <= 4.2 else 'outside'}</td></tr>"
        f"<tr><th>attack angle</th><td class='n'>{m['ang']}&deg;</td>"
        f"<td class='v {'ok' if m['ang'] >= 150 else 'no'}'>"
        f"{'near-linear' if m['ang'] >= 150 else 'too bent'}</td></tr>"
        f"<tr><th>attack-ready</th><td class='n'>{m['ar']:.2f}</td>"
        f"<td class='v {'ok' if m['ar'] > 0.2 else 'no'}'>"
        f"{'carries' if m['ar'] > 0.2 else 'marginal'}</td></tr>"
        f"</table><p class='mnote'>{m['note']}</p></div>"
        for m in MODES)

    ranked = sorted(MODES, key=lambda m: -m["ar"])
    rank_rows = "".join(
        f"<tr><td class='n'>{i + 1}</td>"
        f"<td><span class='sw' style='background:{m['col']}'></span>mode {m['key']}</td>"
        f"<td class='n'>{m['n']}</td><td class='n'>{m['d']:.1f}</td>"
        f"<td class='n'>{m['ang']}</td><td class='n'>{m['ar']:.2f}</td></tr>"
        for i, m in enumerate(ranked))

    # 10 ns sweep: episodes in/out of attack geometry for the winning mode.
    rng2 = random.Random(SEED + 1)
    t, eps = 0.0, []
    while t < 1000:
        w = rng2.uniform(18, 62)
        if rng2.random() < 0.42:
            eps.append((t, min(w, 1000 - t)))
        t += w
    sweep = "".join(
        f"<rect x='{40 + x * 0.63:.1f}' y='26' width='{w * 0.63:.1f}' height='26' "
        f"rx='2' fill='var(--blue)' opacity='.80'/>" for x, w in eps)

    # 100 ns MD: an RMSD trace that stays under the 1.2 nm residence bar.
    rng3 = random.Random(SEED + 2)
    v, path = 0.30, []
    for i in range(121):
        v = max(0.16, min(1.05, v + rng3.gauss(0, 0.055)))
        path.append(f"{40 + i * 5.33:.1f},{126 - v * 62:.1f}")
    trace = " ".join(path)

    return f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>pipeline — how a molecule becomes a row</title><style>
:root{{--ink:#10233f;--navy:#003087;--blue:#0072ce;--blue-pale:#e8f1fb;
 --rule:#ccd6e2;--muted:#5b6b80;--paper:#fff;--raise:#f5f8fc;--card:#fff;
 --good:#0f7a54;--warn:#8a5a00;--bad:#b3261e;
 --prot:#e9eef5;--cav-in:#f7fafd;--cav-out:#dce7f3;--dot:#4a6885;
 --sans:"Helvetica Neue",Helvetica,Arial,system-ui,sans-serif;
 --mono:ui-monospace,"SF Mono",Menlo,Consolas,monospace}}
:root[data-theme="dark"]{{--ink:#dfe7f0;--navy:#8ab4e8;--blue:#6aa9e0;
 --blue-pale:#16283a;--rule:#25333f;--muted:#93a3b4;--paper:#0e151c;
 --raise:#16202a;--card:#131c25;--good:#4fc4a0;--warn:#e0b66a;--bad:#e08a70;
 --prot:#1b2530;--cav-in:#141d26;--cav-out:#1e2b38;--dot:#7f9ab5}}
*{{box-sizing:border-box}}
body{{margin:0;padding:24px 28px 60px;background:var(--paper);color:var(--ink);
 font-family:var(--sans);font-size:14px;line-height:1.55;
 font-variant-numeric:tabular-nums;max-width:1180px}}
h1{{font-size:1.2rem;color:var(--navy);margin:0 0 3px}}
.sub{{color:var(--muted);margin:0 0 4px}}
.warnbar{{border-left:3px solid var(--warn);background:var(--raise);
 padding:9px 13px;margin:14px 0 26px;border-radius:0 4px 4px 0;font-size:13px}}
.step{{display:grid;grid-template-columns:446px 1fr;gap:26px;align-items:start;
 padding:22px 0;border-top:1px solid var(--rule)}}
@media(max-width:900px){{.step{{grid-template-columns:1fr}}}}
.dia{{width:100%;height:auto;background:var(--card);border:1px solid var(--rule);
 border-radius:6px}}
.mini{{width:100%;height:auto;background:var(--card);border:1px solid var(--rule);
 border-radius:5px;display:block}}
h2{{font-size:.95rem;color:var(--navy);margin:2px 0 6px}}
.n0{{font-family:var(--mono);font-size:.6rem;letter-spacing:.14em;
 text-transform:uppercase;color:var(--blue);font-weight:700;margin:0 0 3px}}
p{{margin:.45em 0}}
.lbl{{font:600 10px var(--mono);fill:var(--muted)}}
.cap{{font:10.5px var(--sans);fill:var(--muted)}}
.mtag{{font:700 11px var(--mono);letter-spacing:.06em}}
.modes{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}}
@media(max-width:760px){{.modes{{grid-template-columns:1fr}}}}
.mcard{{border:1px solid var(--rule);border-radius:6px;padding:9px;
 background:var(--raise)}}
table{{border-collapse:collapse;width:100%;font-size:12.5px}}
.crit{{margin-top:7px}}
.crit th{{text-align:left;font-weight:500;color:var(--muted);padding:2px 0;
 font-size:11.5px}}
.crit td{{padding:2px 0}}
td.n,th.n{{text-align:right;font-family:var(--mono)}}
.v{{text-align:right;font-size:10.5px;font-weight:700;padding-left:8px}}
.v.ok{{color:var(--good)}} .v.no{{color:var(--bad)}}
.mnote{{font-size:11.5px;color:var(--muted);margin:7px 0 0;line-height:1.4}}
.rank th{{font-family:var(--mono);font-size:.58rem;letter-spacing:.1em;
 text-transform:uppercase;color:var(--muted);text-align:right;padding:5px 8px;
 border-bottom:1px solid var(--rule)}}
.rank th:nth-child(2){{text-align:left}}
.rank td{{padding:5px 8px;border-bottom:1px solid var(--rule)}}
.sw{{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:6px}}
.tl{{width:100%;height:auto;background:var(--card);border:1px solid var(--rule);
 border-radius:6px}}
code{{font-family:var(--mono);font-size:12.5px;background:var(--raise);
 padding:1px 5px;border-radius:3px}}
.step.wide{{grid-template-columns:1fr}}
.full{{min-width:0}}
.p3wrap{{margin-top:10px}}
.p3box{{position:relative;width:100%;height:460px;border:1px solid var(--rule);
 border-radius:6px;overflow:hidden;background:#eef1f6}}
.p3box > div{{position:absolute;inset:0}}
.p3ctl{{display:flex;align-items:center;gap:7px;flex-wrap:wrap;margin-bottom:8px}}
.p3b{{font:11.5px var(--sans);padding:3px 10px;border:1px solid var(--rule);
 background:var(--card);color:var(--muted);border-radius:99px;cursor:pointer;
 display:inline-flex;align-items:center;gap:6px}}
.p3b.on{{border-color:var(--navy);color:var(--navy);font-weight:600}}
.p3b b{{font-family:var(--mono);font-size:10.5px}}
.dotc{{width:9px;height:9px;border-radius:50%;display:inline-block}}
.p3l{{font-size:11.5px;color:var(--muted);margin-left:4px}}
.p3cap{{font-size:11.5px;color:var(--muted);margin-top:8px;line-height:1.45}}
.p3wait{{position:absolute;inset:0;display:flex;align-items:center;
 justify-content:center;font:12px var(--mono);color:var(--muted);
 background:#eef1f6;z-index:2}}
.foot{{margin-top:30px;padding-top:14px;border-top:1px solid var(--rule);
 font-size:12.5px;color:var(--muted)}}
</style></head><body>

<h1>How a molecule becomes a row</h1>
<p class="sub">Docking &rarr; modes &rarr; criteria &rarr; ranking &rarr; sweep &rarr; MD.</p>
<div class="warnbar"><strong>The geometry on this page is drawn, not measured.</strong>
The poses, clusters and traces are illustrative — they are there to show the shape of
the pipeline, and no number positioned on a diagram is a result. Every
<em>parameter</em> named is real: 500 runs, the 2.8&ndash;4.2&nbsp;&Aring; near-attack
window, the 150&deg; angular bar, 10&nbsp;ns and 100&nbsp;ns, and the
1.2&nbsp;nm residence cut.</div>

<div class="step">
 <div>{_stage1(pts)}</div>
 <div><p class="n0">Step 1 &middot; dock</p>
  <h2>500 poses, one molecule</h2>
  <p>We dock each molecule into Pin1 <strong>500 separate times</strong>. Every run
  searches independently, so we get a cloud of possible placements, not one answer.</p>
  <p>500 because we measured it: ~300 runs covers 95% of the poses. 500 leaves margin.</p>
  <p>The search works. The right pose is somewhere in this cloud <strong>93.3%</strong>
  of the time. It gets lost later.</p></div>
</div>

<div class="step">
 <div>{_stage2(pts)}</div>
 <div><p class="n0">Step 2 &middot; split</p>
  <h2>The cloud is several binding modes</h2>
  <p>We group the poses into <strong>modes</strong> by where the reactive atom sits
  and which way the warhead points.</p>
  <p>Not by whole-molecule shape — two poses can differ in a far-off ring and still
  be the same mode. Not by docking energy, which we know carries no signal here.</p>
  <p>Each mode then becomes its own row in the GUI.</p></div>
</div>

<div class="step wide">
 <div class="full"><p class="n0">The same two steps, on real output</p>
  <h2>{real.get('n', 0)} real poses of one screened molecule</h2>
  <p>Drag to rotate. Switch modes to see the cloud split{counts_txt}.</p>
  {real_block}</div>
</div>

<div class="step">
 <div class="modes">{panels}</div>
 <div><p class="n0">Step 3 &middot; criteria</p>
  <h2>Can this mode actually react?</h2>
  <p>Two checks, both must pass. <strong>Distance</strong>: the reactive carbon has
  to sit <code>2.8&ndash;4.2 &Aring;</code> from the sulfur. Closer means the bond
  already formed; further means no reaction.</p>
  <p><strong>Angle</strong>: it has to come in roughly head-on.</p>
  <p>Mode C shows why you need both — right distance, wrong angle. Distance alone
  would have passed it.</p></div>
</div>

<div class="step">
 <div><table class="rank">
  <tr><th class="n">#</th><th>mode</th><th class="n">poses</th>
      <th class="n">d &Aring;</th><th class="n">angle</th><th class="n">ready</th></tr>
  {rank_rows}</table>
  <p class="mnote" style="margin-top:10px">Ranked on attack-readiness, not on
  docking energy — energy correlates with reaction competence at
  &rho;&nbsp;=&nbsp;+0.009 across 115,300 poses, which is noise.</p></div>
 <div><p class="n0">Step 4 &middot; rank</p>
  <h2>The biggest mode does not automatically win</h2>
  <p>We rank on geometry, not on how many poses landed in a mode. The two can
  disagree, and when they do we carry which mode was picked on the row.</p>
  <p>Picking this way finds the right pose <strong>93.3%</strong> of the time.
  Picking by docking energy: <strong>60.0%</strong>.</p>
  <p>Every ranking is still stamped <code>rank_validated = False</code>. It is an
  ordering we produced, not proof the top molecules bind.</p></div>
</div>

<div class="step">
 <div><svg viewBox="0 0 700 150" class="tl" role="img"
  aria-label="A 10 ns sweep with attack-ready episodes marked">
  <text x="40" y="18" class="cap">10 ns sweep &middot; attack-ready episodes</text>
  <rect x="40" y="26" width="630" height="26" rx="2" fill="var(--cav-out)"/>
  {sweep}
  <line x1="40" y1="60" x2="670" y2="60" stroke="var(--rule)"/>
  <text x="40" y="74" class="cap">0</text><text x="640" y="74" class="cap">10 ns</text>
  <text x="40" y="104" class="cap">100 ns MD &middot; ligand RMSD</text>
  <polyline points="{trace}" fill="none" stroke="var(--blue)" stroke-width="1.5"/>
  <line x1="40" y1="{126 - 1.2 * 62:.1f}" x2="670" y2="{126 - 1.2 * 62:.1f}"
        stroke="var(--bad)" stroke-width="1" stroke-dasharray="4 3"/>
  <text x="676" y="{130 - 1.2 * 62:.1f}" class="cap" text-anchor="end"
        style="fill:var(--bad)">1.2 nm</text>
  <line x1="40" y1="132" x2="670" y2="132" stroke="var(--rule)"/>
 </svg></div>
 <div><p class="n0">Step 5 &middot; sweep, then MD</p>
  <h2>Does it hold up once things move?</h2>
  <p>A docked pose is frozen. The <strong>10 ns sweep</strong> lets it move and asks
  how much of the time it still looks ready to react. This is <em>triage</em> — it
  picks what is worth a long run.</p>
  <p>Those go to <strong>100 ns MD</strong>, which asks a different question: does
  the molecule stay on target at all? That is what the GUI ranks on.</p>
  <p>A molecule can be attack-ready and still leave. It can also sit there for
  100 ns facing the wrong way.</p></div>
</div>

<p class="foot">Parameters and measured results shown here come from the 2.2.0
framework and its decision records. The controls that test whether this criterion
recognises chemistry known to react are on the
<a href="controls.html">controls page</a> — read them next; they qualify everything
above.</p>
</body></html>"""
