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

from . import run_paths as rp
from . import target_config as tc

#: Sweep length, DERIVED (@tt8804: "update the gui to say 8 ns sweep not
#: 10"). The deck said 10 ns in six places while the sweep has run at 8 ns
#: since D0085 -- an explainer that misstates the spec teaches the wrong run.
_SWEEP_NS = int(round(tc.md_sweep_ps() / 1000))

import collections
import json
import math
import random
from pathlib import Path

SEED = 7
N_POSES = 500

#: The molecule the worked example is drawn from.
#:
#: PINNED WITHIN A RUN, CHOSEN PER RUN. It was a literal -- a nac_v3 molecule --
#: so after a topic bump the page that explains the CURRENT screen had no cloud
#: to draw and fell back to schematic dots only. `example_molecule()` picks
#: deterministically from THIS topic's clouds, so the picture is the same on
#: every rebuild of a given run and still exists after the next one.
EXAMPLE_PIN: str | None = None
#: THIS RUN's pose clouds and receptor. Both were literals -- `nac_v3_allposes`
#: and `3IKD_noligand.pdb` -- so the page that explains the CURRENT screen drew
#: its worked example from a superseded one and named a receptor the config no
#: longer has to agree with.
ALLPOSES = str(rp.allposes_dir())
RECEPTOR = str(rp.receptor_prep())
THREEDMOL = "scripts/.cache_3dmol-min.js"


def example_molecule() -> str | None:
    """A molecule from this run to draw the worked example from.

    Wants a cloud that actually shows the point: several hundred poses and at
    least three modes, so "the mess is several binding modes" has three visible
    clusters rather than one. Chosen by sorted name among the candidates that
    qualify -- deterministic, so the page does not change picture on every
    rebuild, and derived, so it survives a topic bump.
    """
    import glob as _g
    if EXAMPLE_PIN:
        return EXAMPLE_PIN
    best = None
    for f in sorted(_g.glob(str(rp.allposes_dir() / "*.sdf")))[:60]:
        try:
            poses = _read_sdf(Path(f))
        except Exception:                                  # noqa: BLE001
            continue
        if len(poses) < 200:
            continue
        modes = {q["mode"] for q in poses}
        if len(modes) >= 3:
            return Path(f).stem
        if best is None:
            best = Path(f).stem
    return best


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
        try:
            nb = int(lines[ci][3:6])
        except ValueError:
            nb = 0
        # Keep every atom's coordinate AND its original index, because the bond
        # block numbers atoms including hydrogens. Filtering H first and then
        # reading the bond block would silently shift every bond by however many
        # hydrogens preceded it -- bonds drawn between the wrong atoms, which
        # still looks like a molecule.
        allpos, heavy = [], {}
        for ai, ln in enumerate(lines[ci + 1:ci + 1 + n], start=1):
            try:
                x, y, z = float(ln[0:10]), float(ln[10:20]), float(ln[20:30])
            except ValueError:
                allpos.append(None)
                continue
            allpos.append((x, y, z))
            if ln[31:34].strip() != "H":
                heavy[ai] = (x, y, z)
        if not heavy:
            continue
        bonds = []
        for ln in lines[ci + 1 + n:ci + 1 + n + nb]:
            try:
                a, b = int(ln[0:3]), int(ln[3:6])
            except ValueError:
                continue
            if a in heavy and b in heavy:
                bonds.append((a, b))
        xs = sum(p[0] for p in heavy.values())
        ys = sum(p[1] for p in heavy.values())
        zs = sum(p[2] for p in heavy.values())
        k = len(heavy)
        mode = "0"
        for i, ln in enumerate(lines):
            if ln.startswith(">") and "<mode>" in ln and i + 1 < len(lines):
                mode = lines[i + 1].strip()
                break
        out.append({"mode": mode, "c": [xs / k, ys / k, zs / k],
                    "atoms": heavy, "bonds": bonds,
                    "blk": "\n".join(lines) + "\n$$$$\n"})
    return out


def _basis(poses):
    """One shared 2D viewing frame for every panel.

    Computed over ALL poses at once so the panels are comparable: a per-panel
    projection would give each mode its own axes and the reader would be
    comparing pictures taken from different angles.
    """
    import numpy as np
    pts = np.array([p for q in poses for p in q["atoms"].values()], dtype=float)
    ctr = pts.mean(axis=0)
    u, s, vt = np.linalg.svd(pts - ctr, full_matrices=False)
    return ctr, vt[0], vt[1]


def _medoid(poses, mode=None):
    """The pose sitting closest to the middle of its group.

    A representative chosen by geometry rather than by file order — `poses[0]`
    would be whichever the docking happened to write first, which is the
    take-it-by-position habit this project keeps paying for.
    """
    import numpy as np
    sel = [p for p in poses if mode is None or p["mode"] == mode]
    if not sel:
        return None
    c = np.array([p["c"] for p in sel], dtype=float)
    return sel[int(np.argmin(((c - c.mean(axis=0)) ** 2).sum(axis=1)))]


def _sg_coord():
    """Cys113's sulfur, from the receptor, found by residue identity.

    Every distance in this project is measured to this atom, so the picture that
    shows the poses should show what they are aimed at. Matched on residue name,
    number and atom name rather than a line offset -- the same file also carries
    Cys57's SG, and taking the first SG in the file would silently anchor the
    whole diagram to the wrong cysteine.
    """
    p = Path(RECEPTOR)
    if not p.is_file():
        return None
    for ln in p.read_text(errors="ignore").splitlines():
        if (ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "SG"
                and ln[17:20].strip() == "CYS" and ln[22:26].strip() == "113"):
            return (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return None


def _pocket_bg(w, h) -> str:
    """A schematic pocket behind the poses, so they read as sitting IN something.

    Drawn, not measured — it is an outline at the extent of the cloud, not the
    molecular surface, and the caption says so.
    """
    return (f"<ellipse cx='{w/2:.0f}' cy='{h/2:.0f}' rx='{w*0.44:.0f}' "
            f"ry='{h*0.42:.0f}' fill='var(--cav-out)' stroke='var(--rule)' "
            f"stroke-width='.8'/>"
            f"<ellipse cx='{w/2:.0f}' cy='{h/2:.0f}' rx='{w*0.33:.0f}' "
            f"ry='{h*0.31:.0f}' fill='none' stroke='var(--rule)' "
            f"stroke-width='.6' stroke-dasharray='2 4' opacity='.7'/>")


def _pose_svg(poses, basis, only=None, colour="#0072ce",
              w=250, h=190, stroke=0.55, op=0.5, pocket=False, one=None,
              highlight=None, only_bg=False) -> str:
    """Every pose drawn flat: one <path> per pose, not one line per bond.

    491 poses at ~32 bonds each is ~15,000 elements as lines and 491 as paths,
    for the same picture and a fraction of the document.
    """
    import numpy as np
    ctr, e1, e2 = basis
    sel = [] if only_bg else ([one] if one is not None else
                              [p for p in poses if only is None or p["mode"] == only])
    if not sel and not only_bg:
        return ""
    allxy = np.array([[(np.array(v) - ctr) @ e1, (np.array(v) - ctr) @ e2]
                      for q in poses for v in q["atoms"].values()])
    lo, hi = allxy.min(axis=0), allxy.max(axis=0)
    pad = 8.0
    sx = (w - 2 * pad) / max(1e-6, hi[0] - lo[0])
    sy = (h - 2 * pad) / max(1e-6, hi[1] - lo[1])
    s = min(sx, sy)

    def xy(v):
        d = np.array(v) - ctr
        return (pad + (d @ e1 - lo[0]) * s, h - pad - (d @ e2 - lo[1]) * s)

    paths = []
    for p in sel:
        d = []
        for a, b in p["bonds"]:
            x1, y1 = xy(p["atoms"][a]); x2, y2 = xy(p["atoms"][b])
            d.append(f"M{x1:.1f} {y1:.1f}L{x2:.1f} {y2:.1f}")
        if d:
            paths.append(f"<path d='{''.join(d)}'/>")
    # THE ANCHOR. Cys113's sulfur, projected through the same basis as the poses,
    # so its position on the picture is measured rather than placed.
    anchor = ""
    sg = _sg_coord()
    if sg is not None:
        gx, gy = xy(sg)
        anchor = (f"<circle cx='{gx:.1f}' cy='{gy:.1f}' r='4' fill='#f0c000' "
                  f"stroke='#8a6d00' stroke-width='.7'/>"
                  f"<text x='{gx + 7:.1f}' y='{gy + 3.5:.1f}' class='sglbl'>"
                  f"Cys113 S&gamma;</text>")
    # THE MEDOID DRAWN ON TOP OF ITS OWN CLOUD. A cloud alone does not say which
    # pose represents it, and a medoid alone does not say what it was chosen from.
    # Overlaying the two answers both at once.
    hi_path = ""
    if highlight is not None:
        d = []
        for a, b in highlight["bonds"]:
            x1, y1 = xy(highlight["atoms"][a]); x2, y2 = xy(highlight["atoms"][b])
            d.append(f"M{x1:.1f} {y1:.1f}L{x2:.1f} {y2:.1f}")
        if d:
            hi_path = (f"<path d='{''.join(d)}' fill='none' stroke='{colour}' "
                       f"stroke-width='1.7' stroke-opacity='1' "
                       f"stroke-linecap='round'/>")
    return (f"<svg viewBox='0 0 {w} {h}' class='psvg' role='img' "
            f"aria-label='{len(sel)} docked poses drawn flat around Cys113'>"
            + (_pocket_bg(w, h) if pocket else "")
            + f"<g fill='none' stroke='{colour}' stroke-width='{stroke}' "
            f"stroke-opacity='{op}' stroke-linecap='round'>"
            + "".join(paths) + "</g>" + hi_path + anchor + "</svg>")

#: The three illustrative modes: share of the 500, centre, spread, and the story
#: each one tells. Sizes are deliberately uneven — a dominant mode plus two
#: minority modes is the common real shape, and a molecule promoted on a minority
#: mode is a different claim from one promoted on its dominant mode.
MODES = [
    {"key": "1", "n": 246, "cx": 214, "cy": 150, "sx": 34, "sy": 25,
     "col": "#0072ce", "d": 3.4, "ang": 168, "ar": 0.41, "note":
     "warhead on the sulfur, roughly in line — inside the window and near-linear"},
    {"key": "2", "n": 158, "cx": 292, "cy": 186, "sx": 27, "sy": 20,
     "col": "#7b5ea7", "d": 4.9, "ang": 141, "ar": 0.06, "note":
     "sits deeper in the pocket, warhead too far out to react"},
    {"key": "3", "n": 96, "cx": 176, "cy": 208, "sx": 22, "sy": 17,
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
<!-- ABOVE the sulfur, not beside it: to the right it collided with mode 2's tag -->
<text x='{SG[0]}' y='{SG[1] - 15}' class='lbl' text-anchor='middle'>Cys113 S&gamma;</text>
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


def _stage2b(pts) -> str:
    """The SECOND pass: one mode, cut again on whole-molecule shape (#61).

    Drawn on mode 1 only, and drawn as a cut through an EXISTING hull rather than
    as new blobs, because that is what it is: the warhead grouping is unchanged
    and the cut is inside it. A reader who sees three fresh clouds would think
    the first pass had been redone.
    """
    m = MODES[0]
    cx, cy, sx, sy = m["cx"], m["cy"], m["sx"] * 2.1, m["sy"] * 2.1
    col = m["col"]
    # Two lobes inside the one hull, plus the medoid each contributes.
    lobes = ""
    for i, (dx, dy, lab) in enumerate(((-0.42, -0.18, "1a"), (0.40, 0.20, "1b"))):
        lx, ly = cx + dx * sx, cy + dy * sy
        lobes += (
            f"<ellipse cx='{lx:.0f}' cy='{ly:.0f}' rx='{sx * .52:.0f}' "
            f"ry='{sy * .58:.0f}' fill='{col}' opacity='.16' stroke='{col}' "
            f"stroke-width='1.2'/>"
            f"<circle cx='{lx:.0f}' cy='{ly:.0f}' r='4.2' fill='{col}'/>"
            f"<text x='{lx:.0f}' y='{ly - sy * .58 - 5:.0f}' class='mtag' "
            f"fill='{col}' text-anchor='middle'>{lab}</text>")
    others = "".join(
        f"<ellipse cx='{o['cx']}' cy='{o['cy']}' rx='{o['sx'] * 2.1:.0f}' "
        f"ry='{o['sy'] * 2.1:.0f}' fill='{o['col']}' opacity='.06' "
        f"stroke='{o['col']}' stroke-width='1' stroke-dasharray='3 3'/>"
        f"<text x='{o['cx']}' y='{o['cy'] - o['sy'] * 2.1 - 6:.0f}' class='mtag' "
        f"fill='{o['col']}' opacity='.5'>mode {o['key']}</text>" for o in MODES[1:])
    return f"""<svg viewBox="0 0 430 300" class="dia" role="img"
 aria-label="Mode 1 cut again into sub-modes 1a and 1b on whole-molecule shape">
{_pocket('2b')}
{others}
<ellipse cx='{cx}' cy='{cy}' rx='{sx:.0f}' ry='{sy:.0f}' fill='none'
 stroke='{col}' stroke-width='1' stroke-dasharray='3 3' opacity='.55'/>
{_dots(pts, only=m['key'])}
{lobes}
<text x="14" y="292" class="cap">mode 1 cut again &mdash; anything wider than
 2 &#8491; gets its own row</text>
</svg>"""


#: Mode colours, matching MODES above.
MODE_COLS = {"1": "#0072ce", "2": "#7b5ea7", "3": "#c2703d"}
#: Molecule colours, deliberately unlike the mode colours.
MOL_COLS = {"mol A": "#1b7f79", "mol B": "#a63d7a", "mol C": "#6b7f1b"}


def _blend(a: str, b: str, w: float = 0.55) -> str:
    """Mix two hex colours, w of the first.

    A mode of a molecule is BOTH things at once, so its colour is both: the
    molecule's hue pulled toward the mode's. Every mode-of-molecule then has a
    colour nothing else has, which is what lets one ranked list be carried from
    step 4 into step 5 and still be read row by row.
    """
    def rgb(h):
        h = h.lstrip("#")
        return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))
    ca, cb = rgb(a), rgb(b)
    return "#%02x%02x%02x" % tuple(
        round(ca[i] * w + cb[i] * (1 - w)) for i in range(3))


#: ONE table behind step 4 and step 5. (molecule, mode, geometry score, 10 ns
#: attack-ready, 100 ns engagement or None, held or None).
#:
#: Kept as a single list rather than two, because step 4's ranking and step 5's
#: survival are the SAME rows read twice -- and two hand-maintained lists that
#: must agree is exactly the shape this project keeps finding broken.
#: THREE GATES, NOT ONE, and a row can stop at any of them:
#:   sweep None           -- ranked too low to be swept at all. No 10 ns figure,
#:                           because nothing was run. A sweep number on a row that
#:                           was never swept is a value with no measurement under
#:                           it, which is the defect this repo is named for.
#:   sweep set, eng None  -- swept, and the sweep did not earn it a 100 ns run.
#:   both set             -- ran 100 ns, then either held or left.
#:
#: Proportions follow the real funnel rather than flattering it: 226 swept, 58
#: elevated, so most rows here stop before 100 ns.
#: ENGAGEMENT AND RESIDENCE ARE NEAR-INDEPENDENT (rho = -0.007, #46), so all four
#: combinations appear here on purpose. A picture where every high-engagement row
#: also held would teach the reader a correlation the measurements say is not
#: there -- and the whole reason the GUI can split held from left is that a
#: molecule can do well on one and badly on the other.
ENTRIES = [
    ("mol B", "1", 0.71, 0.52, 0.91, True),    # high engagement, held
    ("mol A", "1", 0.62, 0.41, 0.86, False),   # high engagement, still LEFT
    ("mol B", "2", 0.44, 0.37, 0.29, True),    # low engagement, yet HELD
    ("mol C", "1", 0.33, 0.33, 0.34, False),   # low engagement, left
    ("mol A", "2", 0.15, 0.21, None, None),    # swept, not elevated
    ("mol B", "3", 0.12, 0.11, None, None),    # swept, not elevated
    ("mol C", "2", 0.09, 0.08, None, None),    # swept, not elevated
    ("mol A", "3", 0.07, None, None, None),    # never swept
    ("mol C", "3", 0.04, None, None, None),    # never swept
]


def _ecol(mol: str, mode: str) -> str:
    return _blend(MOL_COLS[mol], MODE_COLS[mode])


def _stage_pool() -> str:
    """Each molecule's modes on the left, all of them ranked together on the right."""
    ranked = sorted(ENTRIES, key=lambda r: -r[2])

    left, y = [], 26
    for mol in MOL_COLS:
        mine = [e for e in ENTRIES if e[0] == mol]
        if not mine:
            continue
        h = 18 + 20 * len(mine)
        left.append(f"<rect x='10' y='{y}' width='150' height='{h}' rx='5' "
                    f"fill='none' stroke='{MOL_COLS[mol]}' stroke-opacity='.5'/>"
                    f"<text x='18' y='{y + 14}' class='mtag' "
                    f"fill='{MOL_COLS[mol]}'>{mol}</text>")
        yy = y + 22
        for _m, mode, geom, *_ in mine:
            c = _ecol(mol, mode)
            left.append(
                f"<rect x='20' y='{yy}' width='130' height='15' rx='3' "
                f"fill='{c}' fill-opacity='.22'/>"
                f"<text x='27' y='{yy + 11}' class='chip' fill='{c}'>mode {mode}</text>"
                f"<text x='143' y='{yy + 11}' class='chip' fill='{c}' "
                f"text-anchor='end'>{geom:.2f}</text>")
            yy += 20
        y += h + 10

    right, yy = [], 26
    for i, (mol, mode, geom, *_) in enumerate(ranked, start=1):
        c = _ecol(mol, mode)
        right.append(
            f"<rect x='250' y='{yy}' width='168' height='17' rx='3' "
            f"fill='{c}' fill-opacity='.22'/>"
            f"<text x='243' y='{yy + 12}' class='chip' fill='var(--muted)' "
            f"text-anchor='end'>{i}</text>"
            f"<text x='257' y='{yy + 12}' class='chip' fill='{c}'>"
            f"{mol} &middot; mode {mode}</text>"
            f"<text x='411' y='{yy + 12}' class='chip' fill='{c}' "
            f"text-anchor='end'>{geom:.2f}</text>")
        yy += 21

    return f"""<svg viewBox="0 0 430 300" class="dia" role="img"
 aria-label="Modes from three molecules pooled into a single ranking">
<text x="10" y="18" class="cap">each molecule's modes</text>
<text x="250" y="18" class="cap">one ranked list</text>
{''.join(left)}
<path d="M175 150 L238 150" stroke="var(--blue)" stroke-width="1.4" fill="none"/>
<path d="M232 145 L240 150 L232 155 Z" fill="var(--blue)"/>
{''.join(right)}
</svg>"""


def _split_counts() -> tuple[str, str]:
    """How many candidates pose splitting turned into how many modes.

    From the screen's own output, and only the two columns that answer it — the
    full frame is ~2.9M rows across 131 files and the page needs none of the rest.
    """
    import glob as _g
    # THIS RUN's screen output. This was hardcoded `nac_v3/`, so the page
    # explaining the CURRENT pipeline stated a superseded run's split: "5,773
    # candidates -> 8,152 modes" against this screen's 561 -> 4,432.
    fs = sorted(_g.glob(str(rp.BLACKSMITH / rp.topic() / "*.csv")))
    if not fs:
        return "", ""
    # PER FILE, not one try around the lot. Half this directory is aggregate
    # output carrying no `mode` column at all, and wrapping the whole concat meant
    # one unusable file returned "no data" for all 131 -- a blank label with the
    # answer sitting in the other 64.
    import pandas as _pd
    parts = []
    for f in fs:
        try:
            parts.append(_pd.read_csv(f, usecols=["parent_ident", "mode"]))
        except Exception:                                  # noqa: BLE001
            continue                                       # aggregate file
    if not parts:
        return "", ""
    d = _pd.concat(parts, ignore_index=True).dropna(subset=["parent_ident", "mode"])
    return (f"{d.parent_ident.nunique():,}",
            f"{len(d.drop_duplicates(['parent_ident', 'mode'])):,}")


def _run_counts() -> dict:
    """How far the real run actually got — swept, elevated, held.

    Read from the outputs rather than typed in, so the arrow on the page cannot
    drift away from the tree the way a hand-written number would.
    """
    import glob as _g
    out = {"swept": 0, "survivors": 0, "md": 0, "held": 0}
    try:
        import pandas as _pd
        sw = _pd.concat([_pd.read_csv(f) for f in rp.sweep_result_files()],
                        ignore_index=True)
        sw = sw[(sw.get("sweep_ps", 0) > 1000) & (sw.status == "ok")]
        out["swept"] = int(sw.parent_ident.nunique())
        # SURVIVORS = the sweep's own rule, a sustained episode. D0076 records
        # that this filter discards the brief approaches n_visits exists to
        # count, so the number is reported as what the sweep DID, not as a
        # statement that the rest are unreactive.
        # SURVIVED = held under the bar, which is what earns a 100 ns run.
        # This counted `n_visits > 0` -- a mode that merely approached the
        # anchor once -- so the arrow claimed a shortlist the cascade does not
        # act on. D0085 made stability the gate.
        try:
            from . import pipeline as _pl
            out["survivors"] = int(len(_pl.survivors()))
        except Exception:                                  # noqa: BLE001
            out["survivors"] = 0
        md = _pd.concat([_pd.read_csv(f) for f in
                         _g.glob(str(rp.residence_dir() / "*.csv"))], ignore_index=True)
        m = md[md.production_ps >= 50000].drop_duplicates("ident", keep="last")
        c = "explicit_ligand_rmsd_nm_max"
        m = m[m[c].notna()]
        out["md"] = int(len(m))
                # THE BAR FROM CONFIG (D0085: 0.35 nm), not the old 1.2 "did not
        # dissociate" reading, which passes essentially everything and made the
        # page's final arrow claim a shortlist three times looser than the rule
        # that actually gates BPMD.
        out["held"] = int((m[c] < tc.md_survivor_rmsd_nm()).sum())
    except Exception:                                      # noqa: BLE001
        pass
    return out


def _stage_survival() -> str:
    """The SAME rows step 4 ranked, carried through the sweep and the 100 ns run.

    Same entries, same blended colours, so a row can be followed from one step to
    the next by eye. Ordered by the 100 ns engagement because that is the rank;
    the sweep column sits beside it as the triage that decided who earned a 100 ns
    run at all, and the two disagree on purpose -- mol B mode 2 sweeps above mol C
    mode 1 and still ranks below it, which is the whole reason the sweep is not
    the ranking.
    """
    # Kept in the step-4 ranking order, numbered the same way, so the list the
    # reader just watched being built is the list they now watch being filtered.
    ordered = sorted(ENTRIES, key=lambda r: -r[2])
    rows = []
    y, ROW = 46, 24
    n_swept = sum(1 for r in ENTRIES if r[3] is not None)

    def x_mark(yy, note):
        return (f"<line x1='246' y1='{yy + 4}' x2='258' y2='{yy + 15}' "
                f"stroke='var(--muted)' stroke-width='1.4' opacity='.65'/>"
                f"<line x1='258' y1='{yy + 4}' x2='246' y2='{yy + 15}' "
                f"stroke='var(--muted)' stroke-width='1.4' opacity='.65'/>"
                f"<text x='268' y='{yy + 13}' class='chip' fill='var(--muted)'>"
                f"{note}</text>")

    for rank, (mol, mode, _geom, sweep, eng, held) in enumerate(ordered, start=1):
        name, col = f"{mol} &middot; mode {mode}", _ecol(mol, mode)
        rows.append(
            f"<text x='34' y='{y + 13}' class='rnum' fill='var(--muted)' "
            f"text-anchor='end'>{rank}</text>"
            f"<rect x='42' y='{y}' width='118' height='18' rx='3' fill='{col}' "
            f"fill-opacity='.16'/>"
            f"<text x='50' y='{y + 13}' class='chip' fill='{col}'>{name}</text>")
        if sweep is None:
            # Never swept: no 10 ns figure, because no 10 ns run happened.
            rows.append(
                f"<text x='232' y='{y + 13}' class='stat' fill='var(--rule)' "
                f"text-anchor='end'>&mdash;</text>"
                + x_mark(y, "ranked too low to sweep"))
        else:
            rows.append(
                f"<text x='232' y='{y + 13}' class='stat' fill='var(--muted)' "
                f"text-anchor='end'>{sweep:.2f}</text>")
            if eng is None:
                rows.append(x_mark(y, "swept, not elevated"))
            else:
                tag, tc = ("held", "var(--good)") if held else ("left", "var(--bad)")
                rows.append(
                    f"<line x1='244' y1='{y + 9}' x2='268' y2='{y + 9}' "
                    f"stroke='{col}' stroke-width='1.3' marker-end='url(#ah2)'/>"
                    f"<text x='320' y='{y + 13}' class='stat' fill='var(--navy)' "
                    f"text-anchor='end'>{eng * 100:.0f}%</text>"
                    f"<rect x='336' y='{y}' width='56' height='18' rx='3' "
                    f"fill='{tc}' fill-opacity='.16'/>"
                    f"<text x='364' y='{y + 13}' class='chip' fill='{tc}' "
                    f"text-anchor='middle'>{tag}</text>")
        y += ROW
    # The dashed line is the SWEEP's reach: everything below it was never swept.
    cut = 46 + n_swept * ROW - 3
    return f"""<svg viewBox="0 0 430 {y + 16}" class="dia" role="img"
 aria-label="Ranked modes with the sweep and 100 ns numbers behind the rank">
<defs><marker id="ah2" markerWidth="7" markerHeight="7" refX="6" refY="3.5"
 orient="auto"><path d="M0 0 L7 3.5 L0 7 z" fill="var(--blue)" opacity=".7"/></marker></defs>
<text x="42" y="20" class="cap">ranked mols</text>
<text x="232" y="20" class="cap" text-anchor="end">{_SWEEP_NS} ns</text>
<text x="320" y="20" class="cap" text-anchor="end">100 ns</text>
<text x="336" y="20" class="cap">outcome</text>
<text x="232" y="33" class="cap2" text-anchor="end">attack-ready</text>
<text x="320" y="33" class="cap2" text-anchor="end">engaged</text>
<!-- the sweep runs DOWN the list: it decides how far down the cut falls -->
<line x1="14" y1="46" x2="14" y2="{cut}" stroke="var(--blue)" stroke-width="1.2"
 opacity=".65" marker-end="url(#ah2)"/>
<text x="10" y="{(46 + cut) / 2:.0f}" class="cap2" fill="var(--blue)"
 transform="rotate(-90 10 {(46 + cut) / 2:.0f})" text-anchor="middle">{_SWEEP_NS} ns sweep</text>
<line x1="26" y1="{cut + 3}" x2="400" y2="{cut + 3}" stroke="var(--rule)"
 stroke-dasharray="3 3"/>
{''.join(rows)}
</svg>"""


def _paired(dots_svg: str, real_svg: str, cap_l: str, cap_r: str) -> str:
    """Schematic beside the real thing, so the abstraction is legible.

    A dot cloud alone asks the reader to take on faith that a dot stands for a
    molecule. Putting the actual poses next to it, in the same colour and at the
    same scale, makes the abstraction checkable instead.
    """
    if not real_svg:
        return dots_svg
    return (f"<div class='pair'>"
            f"<figure>{dots_svg}<figcaption>{cap_l}</figcaption></figure>"
            f"<figure>{real_svg}<figcaption>{cap_r}</figcaption></figure></div>")


def _mode_panel(pts, m) -> str:
    """One mode: its own share of the cloud, plus the representative it elects.

    The transform is FITTED to what this panel actually draws -- the mode's own
    spread and the sulfur it points at -- rather than being one hard-coded offset
    shared by all three. A fixed transform framed the biggest mode and pushed the
    smaller two off their own canvases: dots, pocket and labels all clipped, on a
    diagram whose whole job is to show where a mode sits.
    """
    W, H, PAD = 200.0, 150.0, 20.0
    xs = [m["cx"] - 3.3 * m["sx"], m["cx"] + 3.3 * m["sx"], SG[0]]
    ys = [m["cy"] - 3.3 * m["sy"], m["cy"] + 3.3 * m["sy"], SG[1]]
    lo_x, hi_x, lo_y, hi_y = min(xs), max(xs), min(ys), max(ys)
    sc = min((W - 2 * PAD) / max(1.0, hi_x - lo_x),
             (H - 2 * PAD) / max(1.0, hi_y - lo_y))
    ox = PAD - lo_x * sc + (W - 2 * PAD - (hi_x - lo_x) * sc) / 2
    oy = PAD - lo_y * sc + (H - 2 * PAD - (hi_y - lo_y) * sc) / 2
    cx, cy = ox + m["cx"] * sc, oy + m["cy"] * sc
    ang = math.radians(180 - m["ang"])
    sgx, sgy = ox + SG[0] * sc, oy + SG[1] * sc
    wx = cx + 30 * math.cos(ang)
    wy = cy - 30 * math.sin(ang)
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


#: Sulfopin, the covalent parent every T_3/T_4 molecule is grown from.
SULFOPIN = "CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl"


def _svg(smiles: str, w: int, h: int) -> str:
    """A 2D depiction as a data URI, or "" if it will not parse."""
    import base64
    import re as _re
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem, Draw
    except ImportError:
        return ""
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return ""
    # CoordGen, not Compute2DCoords. The default layout puts visibly wrong angles
    # on substituted centres -- @tt8804 spotted it on Sulfopin's R group, where the
    # tert-butyl came out skewed rather than at clean tetrahedral-looking angles.
    # The STRUCTURE was right (verified identical to the 6VAJ free SMILES and to
    # the pose sidecar, C11H20ClNO3S); only the drawing was. CoordGen is the
    # template-based layout and gives conventional geometry.
    try:
        from rdkit.Chem import rdCoordGen
        rdCoordGen.AddCoords(m)
    except Exception:                                      # noqa: BLE001
        AllChem.Compute2DCoords(m)
    d = Draw.rdMolDraw2D.MolDraw2DSVG(w, h)
    d.drawOptions().bondLineWidth = 1
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    s = _re.sub(r"<\?xml.*?\?>", "", d.GetDrawingText(), flags=_re.S)
    s = _re.sub(r"<!--.*?-->", "", s, flags=_re.S)
    return "data:image/svg+xml;base64," + base64.b64encode(s.encode()).decode()


def _chemspace(n: int = 8) -> str:
    """Sulfopin in the middle, real generated derivatives around it.

    The satellites are actual rows from the newest T_4 frame, not drawn examples,
    so what this shows is the chemistry the pipeline really produced rather than
    an artist's impression of it.
    """
    import glob as _g
    core = _svg(SULFOPIN, 190, 130)
    if not core:
        return ""
    sats = []
    fs = [str(x) for x in rp.frames("T4")]
    if fs:
        try:
            import pandas as _pd
            d = _pd.read_parquet(fs[-1]).drop_duplicates("candidate_id")
            d = d[d.canonical_smiles.notna()]
            # ONE PER WARHEAD CLASS (@tt8804). Taking the first n rows returned
            # nine chloroacetamides, which showed the R-group axis and hid the
            # warhead axis -- and T_4 is a warhead x R-group enumeration, so half
            # the design space was invisible.
            if "warhead_class" in d.columns:
                d = d.drop_duplicates("warhead_class")
            for _, r in d.head(n).iterrows():
                u = _svg(str(r.canonical_smiles), 132, 92)
                if u:
                    sats.append((str(r.get("warhead_class") or r.candidate_id), u))
        except Exception:                                  # noqa: BLE001
            pass
    if not sats:
        return ""
    # RADIAL, WITH ARROWS OUT OF THE MIDDLE (@tt8804). A grid put the parent in a
    # cell like any other and the growth relationship disappeared -- the reader
    # could not tell which molecule everything else came from.
    W, H, CX, CY = 440, 352, 220, 165
    RX, RY = 152, 118
    sw, sh = 76, 54
    parts = []
    n = len(sats)
    for i, (cid, u) in enumerate(sats):
        a = -math.pi / 2 + 2 * math.pi * i / n
        x, y = CX + RX * math.cos(a), CY + RY * math.sin(a)
        # start/end the arrow outside each box so it never crosses a structure
        x0, y0 = CX + 66 * math.cos(a), CY + 52 * math.sin(a)
        x1, y1 = x - (sw / 2 + 5) * math.cos(a), y - (sh / 2 + 5) * math.sin(a)
        parts.append(
            f"<line x1='{x0:.1f}' y1='{y0:.1f}' x2='{x1:.1f}' y2='{y1:.1f}' "
            f"stroke='var(--blue)' stroke-width='1' opacity='.5' "
            f"marker-end='url(#ah)'/>")
        parts.append(
            f"<image href='{u}' x='{x - sw/2:.1f}' y='{y - sh/2:.1f}' "
            f"width='{sw}' height='{sh}' preserveAspectRatio='xMidYMid meet'/>"
            f"<rect x='{x - sw/2:.1f}' y='{y - sh/2:.1f}' width='{sw}' "
            f"height='{sh}' rx='3' fill='none' stroke='var(--rule)'/>"
            f"<text x='{x:.1f}' y='{y + sh/2 + 9:.1f}' class='wlbl' "
            f"text-anchor='middle'>{cid}</text>")
    return f"""
<svg viewBox="0 0 {W} {H}" class="dia" role="img"
 aria-label="Sulfopin at the centre with generated derivatives around it">
<defs><marker id="ah" markerWidth="7" markerHeight="7" refX="6" refY="3.5"
 orient="auto"><path d="M0 0 L7 3.5 L0 7 z" fill="var(--blue)" opacity=".6"/></marker></defs>
{''.join(parts)}
<rect x="{CX-64}" y="{CY-46}" width="128" height="92" rx="5"
 fill="var(--card)" stroke="var(--navy)" stroke-width="1.4"/>
<image href="{core}" x="{CX-60}" y="{CY-44}" width="120" height="74"
 preserveAspectRatio="xMidYMid meet"/>
<text x="{CX}" y="{CY+40}" class="mtag" fill="var(--navy)"
 text-anchor="middle">Sulfopin</text>
<text x="10" y="{H-8}" class="cap">the parent, and {n} of the molecules grown from it</text>
</svg>"""


def _n_candidates() -> str:
    """How many molecules THIS RUN screens.

    It used to sum every candidate in both generation arms -- 7,179 -- which is
    the library, not the run. This screen takes T_4 across three warhead
    families after the chemistry gates: 561. Stating the library where the run
    belongs overstates the funnel's first step by an order of magnitude.
    """
    try:
        from . import pipeline as _pl
        n = len(_pl.scope_idents())
        if n:
            return f"{n:,}"
    except Exception:                                      # noqa: BLE001
        pass
    return _library_size()


def _library_size() -> str:
    """Every molecule in both arms — context, not the run."""
    import glob as _g
    n = 0
    for tier in ("T4", "T3"):
        fs = [str(x) for x in rp.frames(tier)]
        if not fs:
            continue
        try:
            import pandas as _pd
            n += _pd.read_parquet(fs[-1], columns=["candidate_id"]) \
                    .candidate_id.nunique()
        except Exception:                                  # noqa: BLE001
            pass
    return f"{n:,}" if n else "thousands of"


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


def build(title: str = "DWI Derivative Screen", built: str = "") -> str:
    rng = random.Random(SEED)
    pts = _points(rng)
    chem = _chemspace()
    n_cand = _n_candidates()
    rc = _run_counts()
    n_split, n_modes = _split_counts()
    # THE SPEC, FROM CONFIG. These were prose ("up to five pieces", "the 2 A we
    # claim as pose accuracy"), so retuning the splitter left the page stating
    # the old rule with authority.
    _c = tc.load()
    _st2 = (_c.get("splitting", {}) or {}).get("stage2", {}) or {}
    _cut = _st2.get("cut_diameter_a", 2.0)
    _maxsub = _st2.get("max_sub", 5)
    _minsz = _st2.get("min_mode_size", 12)
    _sp1 = str((_c.get("splitting", {}) or {}).get("stage1", "")).replace("_", " ")
    _sp2 = (f"whole-molecule RMSD, cut {_cut}&nbsp;&Aring;, max {_maxsub} "
            f"sub-modes, min {_minsz} poses"
            if _st2.get("enabled") else "disabled")

    # Real poses, drawn flat. Static by design (@tt8804): the interactive viewer
    # cost 2.5 MB and a hairball; what the page needs is a picture you can put
    # beside the schematic and compare.
    real, real_all, real_by, real_one, real_mess = {}, "", {}, "", ""
    real_med = {}
    try:
        _ex = example_molecule()
        _p = _read_sdf(Path(ALLPOSES) / f"{_ex}.sdf") if _ex else []
    except Exception:                                      # noqa: BLE001
        _p = []
    if _p:
        # THE PALETTE LENGTH IS THE LIMIT, and the code has to say so.
        # `zip(_order, PALETTE)` truncates silently, so with more modes than
        # colours `_cols` was short and the `_cols[m]` below raised
        # `KeyError: '101'` -- which killed `pipeline.html` outright and, with
        # it, `mdprio_combine`, so the MD results page stopped building. Same
        # family as catalogue #19: a fixed-size constant sized for a smaller
        # workload, failing by breaking rather than by raising anything legible.
        #
        # Five overlaid clouds is already the readable maximum, so the fix is to
        # take the five MOST POPULATED modes rather than to cycle colours and
        # produce twenty indistinguishable layers. The count is reported below
        # so the panel never implies it is showing every mode.
        _PALETTE = ["#0072ce", "#7b5ea7", "#c2703d", "#0f7a54", "#b3261e"]
        _all_modes = [m for m, _ in sorted(
            collections.Counter(q["mode"] for q in _p).items(),
            key=lambda kv: -kv[1])]
        _n_modes = len(_all_modes)
        _order = _all_modes[:len(_PALETTE)]
        _cols = dict(zip(_order, _PALETTE))
        _b = _basis(_p)
        real_all = "".join(
            _pose_svg(_p, _b, only=m, colour=_cols[m], w=430, h=300,
                      stroke=0.5, op=0.42).replace("<svg", "<svg style='position:absolute;inset:0'", 1)
            for m in _order)
        # One background layer carrying the pocket AND the Cys113 anchor, under
        # the per-mode clouds -- so the combined panel has the same context as
        # every single-mode panel instead of floating free.
        _bg = _pose_svg(_p, _b, w=430, h=300, pocket=True, only_bg=True) \
            .replace("<svg", "<svg style='position:absolute;inset:0'", 1)
        real_all = (f"<div class='ovl' style='padding-bottom:{300 / 430 * 100:.1f}%'>"
                    f"{_bg}{real_all}</div>")
        # Step 1's picture: every pose, one colour, in the pocket -- the mess as
        # it actually is, before anything has been grouped.
        real_mess = _pose_svg(_p, _b, colour="#4a6885", w=430, h=300,
                              stroke=0.45, op=0.30, pocket=True)
        # ALL THREE MEDOIDS, EACH IN ITS OWN MODE COLOUR. One medoid in one colour
        # said "the cloud reduces to this pose"; the actual claim is that it
        # reduces to one pose PER MODE, and the three sit in different places.
        _meds = "".join(
            _pose_svg(_p, _b, colour=_cols[m], w=430, h=300, stroke=1.7, op=1.0,
                      one=_medoid(_p, m))
            .replace("<svg", "<svg style='position:absolute;inset:0'", 1)
            for m in _order)
        real_one = (f"<div class='ovl' style='padding-bottom:{300 / 430 * 100:.1f}%'>"
                    f"{_bg}{_meds}</div>")
        # Every mode: its cloud in the pocket, with its own medoid picked out on
        # top, and a second panel showing that medoid alone.
        for m in _order:
            _mm = _medoid(_p, m)
            real_by[m] = _pose_svg(_p, _b, only=m, colour=_cols[m], w=200, h=150,
                                   stroke=0.55, op=0.35, pocket=True, highlight=_mm)
            real_med[m] = _pose_svg(_p, _b, colour=_cols[m], w=200, h=150,
                                    stroke=1.6, op=1.0, pocket=True, one=_mm)
        # `modes` is what is DRAWN; `modes_total` is what exists. They differed
        # silently once the palette cap bit, and a panel that says "modes: 5"
        # for a molecule with 213 of them is a number computed from the display
        # rather than from the data.
        real = {"n": len(_p), "modes": len(_order), "modes_total": _n_modes,
                "counts": dict(collections.Counter(q["mode"] for q in _p)),
                "order": _order}
    # Built outside the f-string: a dict literal inside an f-string expression is
    # read as a set of a set, which is unhashable and fails at build time.
    _c = real.get("counts") or {}
    counts_txt = ("" if not _c else
                  " &mdash; " + " / ".join(f"mode {k}: {v}" for k, v in _c.items()))

    _rk = real.get("order") or []
    panels = "".join(
        f"<div class='mcard'>"
        f"<div class='trio'>"
        f"<figure>{_mode_panel(pts, m)}<figcaption>schematic</figcaption></figure>"
        f"<figure>{real_by.get(_rk[i], '') if i < len(_rk) else ''}"
        f"<figcaption>real poses &mdash; medoid picked out</figcaption></figure>"
        f"<figure>{real_med.get(_rk[i], '') if i < len(_rk) else ''}"
        f"<figcaption>the medoid alone</figcaption></figure></div>"
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
        for i, m in enumerate(MODES))

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
        path.append(f"{40 + i * 5.33:.1f},{176 - v * 62:.1f}")
    trace = " ".join(path)

    return f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>
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
 font-variant-numeric:tabular-nums;max-width:1464px}}
h1{{font-size:1.2rem;color:var(--navy);margin:0 0 3px}}
.sub{{color:var(--muted);margin:0 0 4px}}
.step{{display:grid;grid-template-columns:minmax(0,704px) 1fr;gap:26px;align-items:start;
 padding:26px 0;border-top:2.5px solid var(--rule)}}
@media(max-width:1040px){{.step{{grid-template-columns:1fr}}}}
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
.sglbl{{font:600 8.5px var(--mono);fill:#8a6d00}}
.wlbl{{font:600 8px var(--mono);fill:var(--muted)}}
.chip{{font:10.5px var(--mono)}}
.rnum{{font:700 11px var(--mono)}}
.stat{{font:700 11.5px var(--mono)}}
.cap2{{font:9px var(--sans);fill:var(--muted)}}
/* One mode per ROW (@tt8804): three columns squeezed the diagrams to thumbnails.
   Stacked, each mode gets the schematic, the real poses and its numbers side by
   side at a size you can actually read. */
.modes{{display:grid;grid-template-columns:1fr;gap:9px;margin-top:12px}}
.mcard{{border:1px solid var(--rule);border-radius:6px;padding:9px 11px;
 background:var(--raise);display:grid;
 grid-template-columns:minmax(0,540px) minmax(0,300px);gap:20px;align-items:center;
 justify-content:start}}
@media(max-width:820px){{.mcard{{grid-template-columns:1fr}}}}
.mcard .crit{{margin-top:0}}
table{{border-collapse:collapse;width:100%;font-size:12.5px}}
.crit{{margin-top:7px;font-size:14px}}
.crit th{{font-size:13px !important;color:var(--ink) !important}}
.crit td.n{{font-size:15.5px;font-weight:700;color:var(--navy);padding:3px 6px}}
.crit .v{{font-size:12px}}
.crit th{{text-align:left;font-weight:500;color:var(--muted);padding:2px 0;
 font-size:11.5px}}
.crit td{{padding:2px 0}}
td.n,th.n{{text-align:right;font-family:var(--mono)}}
.v{{text-align:right;font-size:10.5px;font-weight:700;padding-left:8px}}
.v.ok{{color:var(--good)}} .v.no{{color:var(--bad)}}
.mnote{{font-size:12.5px;color:var(--muted);margin:7px 0 0;line-height:1.4}}
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
/* Chemical space: the parent in the middle, real generated derivatives around it. */
.cs{{display:grid;grid-template-columns:206px 1fr;gap:14px;align-items:center;
 border:1px solid var(--rule);border-radius:6px;padding:12px;background:var(--card)}}
@media(max-width:640px){{.cs{{grid-template-columns:1fr}}}}
.csmid{{text-align:center;border-right:1px solid var(--rule);padding-right:12px}}
@media(max-width:640px){{.csmid{{border-right:0;border-bottom:1px solid var(--rule);
 padding:0 0 10px}}}}
.csmid img{{width:100%;height:auto;background:#fff;border-radius:3px}}
.csmid span{{display:block;font:600 10.5px var(--mono);color:var(--navy);margin-top:4px}}
.csgrid{{display:grid;grid-template-columns:repeat(4,1fr);gap:8px}}
@media(max-width:900px){{.csgrid{{grid-template-columns:repeat(2,1fr)}}}}
.sat{{text-align:center}}
.sat img{{width:100%;height:auto;background:#fff;border:1px solid var(--rule);
 border-radius:3px}}
.sat span{{display:block;font:9.5px var(--mono);color:var(--muted);margin-top:2px;
 overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}
.arrow{{text-align:center;padding:2px 0 0}}
.arrow span{{font:600 11px var(--mono);color:var(--blue);letter-spacing:.04em}}
.arrow.down{{padding:10px 0 8px;text-align:left}}
.arrow.down span{{font-size:10.5px;line-height:1.4;display:inline-block}}
.tcap{{font:600 10.5px var(--mono);color:var(--muted);margin:0 0 5px;letter-spacing:.04em;text-transform:uppercase}}
/* Schematic beside the real poses, same width and same colours, so a dot in one
   can be matched to a shape in the other without being told to. */
.pair{{display:grid;grid-template-columns:1fr 1fr;gap:10px;align-items:start}}
/* dots -> cloud -> one pose, left to right. */
.trio{{display:grid;grid-template-columns:repeat(2,1fr);gap:8px 10px;align-items:start}}
.trio figure{{margin:0;min-width:0}}
.trio figcaption{{font-size:9px;color:var(--muted);text-align:center;
 margin-top:3px;line-height:1.3}}
@media(max-width:760px){{.trio{{grid-template-columns:1fr}}}}
.pair figure{{margin:0;min-width:0}}
.pair figcaption{{font-size:10px;color:var(--muted);text-align:center;
 margin-top:3px;line-height:1.3}}
.psvg{{width:100%;height:auto;background:var(--card);border:1px solid var(--rule);
 border-radius:6px;display:block}}
.ovl{{position:relative;width:100%;height:0;background:var(--card);
 border:1px solid var(--rule);border-radius:6px;overflow:hidden}}
.ovl .psvg{{border:0;border-radius:0;background:none;height:100%}}
.foot{{margin-top:30px;padding-top:14px;border-top:1px solid var(--rule);
 font-size:12.5px;color:var(--muted)}}

/* CLICK-THROUGH, ONE STAGE AT A TIME (@tt8804). The page was a single long
   scroll, so the cascade read as one wall rather than as five decisions -- and
   the pose-splitting picture, which is the part people ask about, sat halfway
   down where nobody stopped. */
.slide{{display:none;animation:fade .18s ease-out}}
.slide.on{{display:block}}
@keyframes fade{{from{{opacity:0;transform:translateY(4px)}}to{{opacity:1}}}}
#deck{{position:sticky;top:0;z-index:5;display:flex;align-items:center;gap:10px;
 padding:10px 0 14px;background:var(--bg);border-bottom:1px solid var(--rule);
 margin-bottom:22px;flex-wrap:wrap}}
#deck button{{font:600 12px var(--sans);padding:.35rem .75rem;border-radius:99px;
 border:1px solid var(--rule);background:var(--card);color:var(--ink);cursor:pointer}}
#deck button:hover:not(:disabled){{background:var(--accent);color:#fff;border-color:var(--accent)}}
#deck button:disabled{{opacity:.35;cursor:default}}
.dots{{display:flex;gap:6px;margin-left:auto;flex-wrap:wrap}}
.dot{{font:600 11px var(--sans);padding:.3rem .6rem;border-radius:99px;
 border:1px solid var(--rule);background:var(--card);color:var(--muted);cursor:pointer}}
.dot.on{{background:var(--ink);color:var(--bg);border-color:var(--ink)}}
.dot:hover{{border-color:var(--accent);color:var(--accent)}}
/* The spec beside the picture: what the splitter is configured to do, so the
   slide cannot describe a rule the run is not using. */
.kv2{{display:grid;grid-template-columns:auto 1fr;gap:2px 12px;font-size:12px;
 margin:0 0 12px;padding:8px 10px;background:var(--card);
 border:1px solid var(--rule);border-radius:5px}}
.kv2 dt{{font-family:var(--mono);color:var(--muted)}}
.kv2 dd{{margin:0}}
</style></head><body>

<h1>{title}</h1>
<p class="sub">Ranking schematic &middot; Timothy Wu{f" &middot; {built}" if built else ""}</p>


<div id="deck">
 <button id="prev" onclick="go(CUR-1)">&larr; back</button>
 <button id="next" onclick="go(CUR+1)">next &rarr;</button>
 <span id="where" style="font:600 12px var(--sans);color:var(--muted)"></span>
 <span class="dots" id="dots"></span>
</div>
<section class="slide" data-i="0">
<div class="step">
 <div>{chem}</div>
 <div><p class="n0">Step 0 &middot; chemical space</p>
  <h2>Grow molecules around a known binder</h2>
  <p>We start from <strong>Sulfopin</strong>, which is known to react with Cys113,
  and generate variations around it &mdash; swapping the group it carries while
  keeping the warhead that does the chemistry.</p>
  <p>That gives thousands of candidates. Everything below is how we narrow them.</p></div>
</div>

<div class="arrow"><span>{n_cand} candidates &darr;</span></div>
</section>
<section class="slide" data-i="1">
<div class="step">
 <div>{_paired(_stage1(pts), real_mess, 'schematic — one dot per pose',
               f"real — all {real.get('n', 0)} poses in the pocket")}</div>
 <div><p class="n0">Step 1 &middot; dock</p>
  <h2>500 poses, one molecule</h2>
  <p>We dock each molecule into the Pin1 pocket <strong>500 separate times</strong>,
  anchored on the catalytic <strong>Cys113</strong> &mdash; the sulfur every
  distance in this project is measured to. Every run searches independently, so we
  get a mess of possible placements, not one answer.</p>
  <p>500 because we measured it: ~300 runs covers 95% of the poses. 500 leaves margin.</p>
  <p>The search works. The right pose is somewhere in this cloud <strong>93.3%</strong>
  of the time. It gets lost later.</p></div>
</div>
</section>
<section class="slide" data-i="2">
<div class="step">
 <div><div class="trio">
   <figure>{_stage2(pts)}<figcaption>pass 1 &mdash; by warhead</figcaption></figure>
   <figure>{_stage2b(pts)}<figcaption>pass 2 &mdash; by shape, inside one mode</figcaption></figure>
   <figure>{real_all}<figcaption>real &mdash; all {real.get('n', 0)} poses, by mode</figcaption></figure>
   <figure>{real_one}<figcaption>one pose per mode &mdash; the medoids</figcaption></figure>
  </div></div>
 <div><p class="n0">Step 2 &middot; pose splitting</p>
  <h2>Pose splitting &mdash; the mess is several binding modes</h2>
  <p class="mnote" style="margin:0 0 10px">Two passes: by warhead, then by shape.</p>
  <dl class="kv2">
   <dt>stage 1</dt><dd>{_sp1}</dd>
   <dt>stage 2</dt><dd>{_sp2}</dd>
   <dt>this run</dt><dd>{n_split} molecules &rarr; {n_modes} modes</dd>
  </dl>
  <p>We group the poses into <strong>modes</strong> by where the reactive atom sits
  and which way the warhead points.</p>
  <p>Not by docking energy, which we know carries no signal here.</p>
  <p><strong>Then a second pass, on whole-molecule shape.</strong> The first pass
  is deliberately blind to everything but the warhead, so two poses can place the
  reactive atom identically and hang the rest of the molecule 5&nbsp;&Aring; apart
  and still be one mode &mdash; and one representative would have to stand for
  both. It splits again wherever a mode is wider than the 2&nbsp;&Aring; we claim
  as pose accuracy ({_cut}&nbsp;&Aring;), up to {_maxsub} pieces. Sub-splits are lettered
  &mdash; <strong>1a, 1b</strong> &mdash; and each is a row of its own, ranked and
  simulated on its own merit.</p>
  <p>Measured on 82 Pin1 crystal structures: carrying one representative recovers
  the true pose <strong>22%</strong> of the time, four recovers <strong>39%</strong>
  (14 cases gained, none lost, <em>p</em> = 1&times;10<sup>&minus;4</sup>), and
  past five there is nothing left to gain. Sulfopin is why it exists: its 456 poses
  formed one mode holding both a pose 1.4&nbsp;&Aring; from its crystal structure
  and the pose we kept at 5.1&nbsp;&Aring;.</p></div>
</div>

<div class="arrow"><span>{n_split} candidates &rarr; {n_modes} binding modes &darr;</span></div>
</section>
<section class="slide" data-i="3">
<div class="step wide">
 <div class="full"><p class="n0">Step 3 &middot; criteria</p>
  <h2>Can this mode actually react?</h2>
  <p>Two checks, both must pass. <strong>Distance</strong>: the reactive carbon has
  to sit <code>2.8&ndash;4.2 &Aring;</code> from the sulfur. Closer means the bond
  already formed; further means no reaction.</p>
  <p><strong>Angle</strong>: it has to come in roughly head-on.</p>
  <p>Mode 3 shows why you need both — right distance, wrong angle. Distance alone
  would have passed it.</p>
  <div class="modes">{panels}</div>
  <p class="tcap" style="margin-top:14px">the three modes, scored</p>
  <table class="rank">
  <tr><th class="n">#</th><th>mode</th><th class="n">poses</th>
      <th class="n">d &Aring;</th><th class="n">angle</th><th class="n">ready</th></tr>
  {rank_rows}</table></div>
</div>
</section>
<section class="slide" data-i="4">
<div class="step">
 <div>{_stage_pool()}
  <p class="mnote" style="margin-top:10px">Ranked on attack-readiness, not on
  docking energy — energy correlates with reaction competence at
  &rho;&nbsp;=&nbsp;+0.009 across 115,300 poses, which is noise.</p></div>
 <div><p class="n0">Step 4 &middot; rank</p>
  <h2>Modes compete, not molecules</h2>
  <p>We rank on geometry, not on how many poses landed in a mode. Every mode from
  every molecule goes into <strong>one ranked list</strong>, so the modes of a
  single molecule are not kept together &mdash; mol B here lands at ranks
  <strong>1, 3 and 6</strong>.</p>
  <p>Picking this way finds the right pose <strong>93.3%</strong> of the time.
  Picking by docking energy: <strong>60.0%</strong>. That is picking the right
  <em>mode</em>; picking the right pose <em>inside</em> it is what step 2's second
  pass is for.</p>
  <p>Each row carries <em>which</em> mode was picked, and every ranking is stamped
  <code>rank_validated = False</code> &mdash; an ordering we produced, not proof
  the top molecules bind.</p></div>
</div>
</section>
<section class="slide" data-i="5">
<div class="step">
 <div>{_stage_survival()}
 <svg viewBox="0 0 700 196" class="tl" role="img"
  aria-label="A {_SWEEP_NS} ns sweep with attack-ready episodes, then the 100 ns RMSD trace">
  <text x="40" y="18" class="cap">{_SWEEP_NS} ns sweep &middot; attack-ready episodes</text>
  <rect x="40" y="26" width="630" height="26" rx="2" fill="var(--cav-out)"/>
  {sweep}
  <line x1="40" y1="60" x2="670" y2="60" stroke="var(--rule)"/>
  <text x="40" y="74" class="cap">0</text><text x="640" y="74" class="cap">{_SWEEP_NS} ns</text>
  <!-- The RMSD caption used to sit at y=104, straight through the trace it
       labels. The band now starts below it -- baseline 176 -- so nothing drawn
       reaches the caption line. -->
  <text x="40" y="100" class="cap">100 ns MD &middot; ligand RMSD</text>
  <polyline points="{trace}" fill="none" stroke="var(--blue)" stroke-width="1.5"/>
  <line x1="40" y1="{176 - 1.2 * 62:.1f}" x2="670" y2="{176 - 1.2 * 62:.1f}"
        stroke="var(--bad)" stroke-width="1" stroke-dasharray="4 3"/>
  <text x="670" y="{176 - 1.2 * 62 - 4:.1f}" class="cap" text-anchor="end"
        style="fill:var(--bad)">1.2 nm &mdash; above this it has left</text>
  <line x1="40" y1="178" x2="670" y2="178" stroke="var(--rule)"/>
  <text x="40" y="192" class="cap">0</text><text x="634" y="192" class="cap">100 ns</text>
 </svg></div>
 <div><p class="n0">Step 5 &middot; sweep, then MD</p>
  <h2>Does it hold up once things move?</h2>
  <p>A docked pose is frozen. The <strong>{_SWEEP_NS} ns sweep</strong> lets it move and asks
  how much of the time it still looks ready to react. This is <em>triage</em> — it
  picks what is worth a long run.</p>
  <p>Those go to <strong>100 ns MD</strong>, which asks a different question: does
  the molecule stay on target at all?</p>
  <p>The results GUI ranks on <strong>max ligand RMSD</strong> over that run
  &mdash; how far it ever got from where it started, lowest first &mdash; with
  engagement and held/left shown beside it. <strong>No sweep reading appears in
  that ranking at all:</strong> triage sitting next to a result reads as a second,
  competing score, so it lives on each molecule's own page instead.</p>
  <p>A molecule can be attack-ready and still leave. It can also sit there for
  100 ns facing the wrong way.</p></div>
</div>

<div class="arrow"><span>{rc['swept']} swept &rarr; {rc['survivors']} survived the sweep
 &rarr; {rc['md']} ran 100&nbsp;ns &rarr; {rc['held']} still on target &darr;</span></div>
</section>
<section class="slide" data-i="6">
<div class="step wide">
 <div class="full"><p class="n0">Elevation</p>
  <h2>What comes out</h2>
  <p>What survives 100 ns is handed to a chemist. <strong>Sulfopin sits in the
  candidate frame</strong>, not in a reference list &mdash; it is docked, split,
  ranked and swept by the identical path as everything else, so the screen has to
  place a known nanomolar inhibitor without being told what it is.</p>
  <p><strong>Two views, and the first one is the pipeline's real first step.</strong>
  <em>Ranking</em> shows every molecule and every mode the screen scored, with its
  pose, before anything is simulated &mdash; that is where you decide what earns a
  sweep. <em>Results</em> shows the sweep and the 100&nbsp;ns runs, for the subset
  that got one.</p>
  <p class="mnote">Synthesis is roughly one compound a week, so the deliverable is
  a short list, not a score.</p></div>
</div>

<p class="foot">Parameters and measured results shown here come from the 2.2.0
framework and its decision records. The controls that test whether this criterion
recognises chemistry known to react are on the
<a href="controls.html">controls page</a> — read them next; they qualify everything
above.</p>
</section>
<script>
// The stage names come from each slide's own "Step N - label" line, so the
// stepper cannot drift from the slide it points at.
const SLIDES = Array.from(document.querySelectorAll('.slide'));
const LABELS = SLIDES.map(function(s, i){{
  const n = s.querySelector('.n0');
  if (!n) return 'step ' + i;
  const t = n.textContent.split('\u00b7');
  return (t.length > 1 ? t[1] : t[0]).trim();
}});
let CUR = 0;
function go(i){{
  if (i < 0 || i >= SLIDES.length) return;
  CUR = i;
  SLIDES.forEach(function(s, k){{ s.classList.toggle('on', k === i); }});
  document.getElementById('prev').disabled = (i === 0);
  document.getElementById('next').disabled = (i === SLIDES.length - 1);
  document.getElementById('where').textContent =
    'stage ' + (i + 1) + ' of ' + SLIDES.length;
  document.querySelectorAll('#dots .dot').forEach(function(d, k){{
    d.classList.toggle('on', k === i); }});
  window.scrollTo({{top: 0, behavior: 'instant'}});
  // Deep-linkable: a stage can be sent to someone directly.
  history.replaceState(null, '', '#' + (i + 1));
}}
document.getElementById('dots').innerHTML = LABELS.map(function(l, i){{
  return '<button class="dot" onclick="go(' + i + ')">' + (i+1) + '. ' + l + '</button>';
}}).join('');
document.addEventListener('keydown', function(e){{
  if (e.key === 'ArrowRight' || e.key === ' ') {{ e.preventDefault(); go(CUR+1); }}
  else if (e.key === 'ArrowLeft') {{ e.preventDefault(); go(CUR-1); }}
}});
go(Math.max(0, Math.min(SLIDES.length - 1, (parseInt(location.hash.slice(1)) || 1) - 1)));
</script>
</body></html>"""
