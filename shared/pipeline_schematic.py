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

import math
import random

SEED = 7
N_POSES = 500

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


def build() -> str:
    rng = random.Random(SEED)
    pts = _points(rng)

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
  <h2>500 poses, one molecule, one pocket</h2>
  <p>Reactive docking against <strong>3IKD</strong> with a flexible Cys113,
  <code>--nrun 500</code>. Each run is an independent search, so the output is a
  cloud of candidate placements rather than one answer.</p>
  <p><strong>Why 500 and not 200.</strong> Pose coverage was measured, not guessed:
  reaching 95% of the reachable poses needs about 300 runs, so 500 is the default
  with margin. Coverage is not the same as reproducibility, and only coverage has
  been demonstrated.</p>
  <p>The search is not the bottleneck. Against 15 deposited Pin1 complexes the
  crystallographic pose is somewhere in this cloud <strong>93.3%</strong> of the
  time. What happens next is where it gets lost.</p></div>
</div>

<div class="step">
 <div>{_stage2(pts)}</div>
 <div><p class="n0">Step 2 &middot; split</p>
  <h2>The cloud is not one binding mode</h2>
  <p>Poses cluster into <strong>modes</strong> on two things: where the reactive
  atom sits, and which way the warhead points.</p>
  <p><strong>Deliberately not</strong> whole-molecule RMSD — two poses that place
  the warhead identically and differ in a distal ring are one mode, not two
  (D0062). <strong>Deliberately not</strong> docking energy, which would re-import
  the exact defect this stage exists to remove (#23/#30). And not distance-to-S&gamma;
  or the attack angle either, because those <em>are</em> the criteria: a mode defined
  partly by its own score is guaranteed to look internally consistent.</p>
  <p>Each mode becomes its own <strong>candidate row</strong>, so ranking, selection
  and this catalogue consume it unchanged.</p></div>
</div>

<div class="step">
 <div class="modes">{panels}</div>
 <div><p class="n0">Step 3 &middot; criteria</p>
  <h2>Each mode is judged on its own geometry</h2>
  <p>Two independent gates, and a mode has to clear both. <strong>Distance</strong>:
  the reactive carbon must sit in the <code>2.8&ndash;4.2 &Aring;</code> near-attack
  window — closer is a formed bond, further is no reaction.
  <strong>Angle</strong>: the approach must be near-linear for an S<sub>N</sub>2
  displacement.</p>
  <p>Mode C is why both are needed. It is comfortably inside the distance window
  and still cannot react, because it approaches side-on at 96&deg;. Distance alone
  would have passed it.</p>
  <p>The angular bar is stricter for S<sub>N</sub>2 than for a perpendicular
  addition, which is why cross-class ranking is biased and the catalogue offers a
  within-class view (#47).</p></div>
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
  <h2>The dominant mode is not automatically the winner</h2>
  <p>Mode A leads here on both population and geometry, but those can disagree, and
  when they do the population is not what decides. A molecule promoted on a minority
  mode is a different claim from one promoted on its dominant mode, so the mode that
  was elevated is carried on the row rather than left implicit.</p>
  <p>Selecting the mode this way recovers the crystal pose <strong>93.3%</strong> of
  the time against docking energy's <strong>60.0%</strong>, measured on 15 crystal
  complexes at 500 runs, each docked twice.</p>
  <p><strong>Every ranking here is stamped <code>rank_validated = False</code>.</strong>
  It is an ordering the pipeline produced, not evidence the molecules at the top
  bind.</p></div>
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
  <h2>Does the geometry survive being moved?</h2>
  <p>A docked pose is a static guess. The <strong>10 ns sweep</strong> asks how much
  of the time the mode actually holds attack geometry once the system is free to
  move — that fraction, and the number of separate <em>visits</em>, are what the
  catalogue's headline <code>ns</code> figure and visit count report.</p>
  <p>Survivors go to <strong>100 ns MD</strong>, which asks a different question:
  does the molecule stay in the pocket at all? Maximum ligand RMSD under
  <code>1.2 nm</code> is <em>held</em>; above it the molecule left.</p>
  <p><strong>These two readings are near-independent</strong> (&rho;&nbsp;=&nbsp;&minus;0.007),
  which is why the catalogue can show them combined or split. A molecule can be
  attack-ready and still leave; it can also sit in the pocket for 100 ns facing the
  wrong way.</p></div>
</div>

<p class="foot">Parameters and measured results shown here come from the 2.2.0
framework and its decision records. The controls that test whether this criterion
recognises chemistry known to react are on the
<a href="controls.html">controls page</a> — read them next; they qualify everything
above.</p>
</body></html>"""
