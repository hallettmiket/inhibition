"""
Purpose: one self-contained HTML for a chemist — structure, SMILES, MD movie, plots.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-10
Input: --candidates <ident...> (each needs a finished 100 ns trajectory)
Output: 00_outputs/blacksmith/shortlist/shortlist_<N>.html

For sending outside the project. It carries no ranking, no gate verdict and no
pipeline commentary — a reader who does not work on this repo cannot check those
claims and does not need them to look at a trajectory. What it carries is the
molecule, its SMILES as selectable text, the 100 ns movie, the RMSD plots, and the
measured values behind them.

SELF-CONTAINED ON PURPOSE. The per-molecule reports are ~9.5 MB each because the
movie frames are base64 in the page; four of them inlined is large, and that is
the price of a file that opens by double-clicking with no server, no directory and
no network. 3Dmol.js is vendored ONCE for the whole document rather than per
viewer.
"""

from __future__ import annotations

import argparse
import base64
import glob
import json
import logging
import re
import sys
from datetime import date
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import gui_shell as gs                      # noqa: E402
from shared import run_paths as rp                      # noqa: E402
from shared import md_movie as mov                      # noqa: E402
from shared import outputs as sout                      # noqa: E402
from shared import report_theme as rt                   # noqa: E402

log = logging.getLogger("shortlist")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
MD = Path("/data/lab_vm/modifiable/inhibition/md_residence_3ikd")
BOUND_NM = 1.2


def _er():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "elevation_report", REPO / "scripts" / "elevation_report.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["elevation_report"] = m
    spec.loader.exec_module(m)
    return m


def smiles_of(ident: str) -> str | None:
    for sub, stem in (("04_t4_combinatorial", "D4"), ("03_t3_reinvent", "D3")):
        fs = sorted(glob.glob(f"/data/lab_vm/append_only/inhibition/{sub}/{stem}_*.parquet"),
                    key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))
        if not fs:
            continue
        fr = pd.read_parquet(fs[-1]).drop_duplicates("candidate_id").set_index("candidate_id")
        if ident in fr.index:
            return str(fr.loc[ident, "canonical_smiles"])
    p = B / f"pose_sidecars/{ident}.json"
    if p.is_file():
        import json
        return json.loads(p.read_text()).get("canonical_smiles")
    return None


def depiction(smi: str, w: int = 340, h: int = 210) -> str:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Draw, rdCoordGen
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return ""
    rdCoordGen.AddCoords(m)
    d = Draw.rdMolDraw2D.MolDraw2DSVG(w, h)
    d.drawOptions().bondLineWidth = 1
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    svg = re.sub(r"<\?xml.*?\?>", "", d.GetDrawingText(), flags=re.S)
    svg = re.sub(r"<!--.*?-->", "", svg, flags=re.S)
    return "data:image/svg+xml;base64," + base64.b64encode(svg.encode()).decode()


#: Heavy-atom contact cutoff, Angstrom. A residue counts as contacting the ligand
#: in a frame if any of its heavy atoms is within this of any ligand heavy atom.
CONTACT_A = 4.5
#: Polar contact: N/O to N/O within this. Called "polar", not "hydrogen bond" --
#: the movie PDB carries no hydrogens, so donor geometry cannot be checked and
#: calling it an H-bond would claim more than was measured.
POLAR_A = 3.5
#: Surface shell, Angstrom: residues with a heavy atom this close to the ligand
#: get a molecular surface. Same reasoning as the movie's shell -- a whole-protein
#: mesh is the expensive call and the far side of the protein is not the subject.
SURF_SHELL_A = 8.0

#: One key, used by both figures, so a colour cannot mean two things across the
#: page. Wording is deliberately about what was MEASURED: "polar", not H-bond;
#: "no polar partner", not hydrophobic contact in the thermodynamic sense.
_LEGEND_KEYS = [
    ("#b3261e", False, "catalytic Cys113",
     f"Sγ to the electrophilic carbon — the bond the screen is for"),
    ("#0f7a54", True, "polar",
     f"N/O to N/O within {POLAR_A} Å"),
    ("#4a6885", False, "close contact",
     f"heavy atoms within {CONTACT_A} Å, no polar partner"),
]

#: The key's own styling, beside the key itself. A second page rendering
#: `_LEGEND_3D` needs these rules, and a page that has the markup but not the
#: styles shows three unlabelled grey lines -- which is worse than no key. Kept
#: here rather than in the page template so there is one definition to import.
KEY3_CSS = """
.key3{display:flex;flex-wrap:wrap;gap:.35rem 1.6rem;margin:.55rem 0 .1rem;
  font:12px var(--sans)}
.key3 .k{display:flex;align-items:center;gap:.45rem;font-weight:600}
.key3 .k i{width:24px;height:3px;border-radius:2px;display:inline-block}
.key3 .k em{font-style:normal;font-weight:400;color:var(--muted);
  margin-left:.4rem}
table.occ{border-collapse:collapse;margin:.7rem 0 .2rem;font-size:12.5px}
table.occ th,table.occ td{padding:.24rem .8rem .24rem 0;text-align:left;
  border-bottom:1px solid var(--rule)}
table.occ th{font:600 10px var(--sans);color:var(--muted);text-transform:uppercase}
table.occ td{font-family:var(--mono)}
h2.ih{font:600 14px var(--sans);margin:.2rem 0 .3rem}
"""

_LEGEND_3D = ('<div class="key3">'
              + "".join(
                  f'<span class="k"><i style="background:{c};'
                  f'{"opacity:.55;" if d else ""}"></i>{n}'
                  f'<em>{w}</em></span>' for c, d, n, w in _LEGEND_KEYS)
              + '</div>')


def _legend_svg(x: float, y: float) -> str:
    """The same key as the 3D view, drawn into the SVG so the figure travels."""
    out = [f"<text x='{x:.0f}' y='{y:.0f}' class='kh'>interactions</text>"]
    for i, (col, dash, name, what) in enumerate(_LEGEND_KEYS):
        yy = y + 17 + i * 15
        dd = " stroke-dasharray='5 4'" if dash else ""
        out.append(
            f"<line x1='{x:.0f}' y1='{yy - 4:.0f}' x2='{x + 26:.0f}' "
            f"y2='{yy - 4:.0f}' stroke='{col}' stroke-width='2.4'{dd}/>"
            f"<text x='{x + 33:.0f}' y='{yy:.0f}' class='kt' fill='{col}'>"
            f"{name}</text>"
            f"<text x='{x + 33 + 7.2 * len(name):.0f}' y='{yy:.0f}' "
            f"class='kw'>{what}</text>")
    return "".join(out)


def ligand_mol(movie_pdb: Path, smi: str):
    """The ligand as an RDKit mol whose atom order matches the PDB's.

    Built FROM the PDB block, then given bond orders from the SMILES template, so
    atom index i here is atom i in the trajectory. That correspondence is what
    lets a contact be attributed to a specific ATOM rather than to the molecule
    as a whole -- which is the difference between an interaction diagram and a
    list of nearby residues.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    first = movie_pdb.read_text().split("ENDMDL")[0]
    blk = "\n".join(l for l in first.splitlines()
                    if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "MOL")
    m = Chem.MolFromPDBBlock(blk, removeHs=True, sanitize=False)
    if m is None:
        return None
    tpl = Chem.MolFromSmiles(smi) if smi else None
    if tpl is not None and tpl.GetNumAtoms() == m.GetNumAtoms():
        try:
            m = AllChem.AssignBondOrdersFromTemplate(tpl, m)
        except Exception:                                  # noqa: BLE001
            pass
    try:
        Chem.SanitizeMol(m)
    except Exception:                                      # noqa: BLE001
        pass
    return m


def contacts(movie_pdb: Path, max_frames: int = 40) -> tuple[list, int]:
    """Per-residue contact frequency across the trajectory.

    Over FRAMES, not one snapshot: a single frame says which residues happened to
    be near the ligand at one instant, which for a molecule that moves is close to
    arbitrary. Frequency over the run is the thing a chemist can act on.

    Returns (rows, n_frames) where each row is
    (resname, resid, fraction_of_frames, polar_fraction).
    """
    import numpy as np
    text = movie_pdb.read_text()
    models = [m for m in text.split("MODEL")[1:]] or [text]
    step = max(1, len(models) // max_frames)
    models = models[::step]

    seen: dict[tuple, list] = {}
    for mdl in models:
        lig, prot = [], []
        for ln in mdl.splitlines():
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            el = ln[76:78].strip() or ln[12:16].strip()[:1]
            if el == "H":
                continue
            try:
                xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                continue
            rn, ri = ln[17:20].strip(), ln[22:26].strip()
            (lig if rn == "MOL" else prot).append((xyz, rn, ri, el))
        if not lig or not prot:
            continue
        L = np.array([a[0] for a in lig])
        P = np.array([a[0] for a in prot])
        d = np.linalg.norm(P[:, None, :] - L[None, :, :], axis=2)
        near = d.min(axis=1) <= CONTACT_A
        lig_pol = np.array([a[3] in ("N", "O") for a in lig])
        hit_res, pol_res = set(), set()
        for i, ok in enumerate(near):
            if not ok:
                continue
            _, rn, ri, el = prot[i]
            hit_res.add((rn, ri))
            if el in ("N", "O") and lig_pol.any():
                if d[i][lig_pol].min() <= POLAR_A:
                    pol_res.add((rn, ri))
        for k in hit_res:
            seen.setdefault(k, [0, 0, {}])
            seen[k][0] += 1
            if k in pol_res:
                seen[k][1] += 1
        # WHICH ATOM. For each contacting residue, the ligand atom it is nearest
        # to in this frame; the mode over frames is the atom the interaction is
        # attributed to. A residue's closest atom can move, so one frame is not
        # enough to name it.
        for i, ok in enumerate(near):
            if not ok:
                continue
            _, rn, ri, _el = prot[i]
            j = int(d[i].argmin())
            tally = seen[(rn, ri)][2]
            tally[j] = tally.get(j, 0) + 1
    n = len(models)
    rows = [(rn, ri, c / n, p / n,
             max(t, key=t.get) if t else None)
            for (rn, ri), (c, p, t) in seen.items()]
    rows.sort(key=lambda r: -r[2])
    return rows, n


def representative_frame(movie_pdb: Path):
    """The frame whose ligand sits closest to the ligand's mean position.

    Not frame 1, which is where the pose STARTED, and not the last, which is
    wherever it happened to stop. A medoid frame is the one a reader is least
    likely to be misled by, and the interaction lines drawn on it are the ones
    that hold for most of the run.
    """
    import numpy as np
    models = movie_pdb.read_text().split("ENDMDL")
    cent, keep = [], []
    for mdl in models:
        lig = [l for l in mdl.splitlines()
               if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "MOL"
               and (l[76:78].strip() or l[12:16].strip()[:1]) != "H"]
        if not lig:
            continue
        xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lig])
        cent.append(xyz.mean(axis=0))
        keep.append(mdl)
    if not keep:
        return None, 0
    c = np.array(cent)
    return keep[int(np.argmin(((c - c.mean(axis=0)) ** 2).sum(axis=1)))], len(keep)


def reactive_atom_index(mol, warhead_class: str) -> int | None:
    """The electrophilic carbon, taken from the warhead class's own SMARTS.

    NOT the ligand atom nearest Cys113, which is what a plain contact search
    returns and which for `bdhi_c5` is the BROMIDE. That is not wrong as a
    distance -- the bromide sits on the carbon under attack -- but drawing the
    Cys113 line to it says the halogen is the interaction, when the reaction is
    S(gamma) displacing that bromide from the carbon. The class table's
    `reactive_atom_smarts` puts the electrophilic carbon first, the same
    convention `shared/covalent_adduct.py` relies on.

    The movie PDB carries no hydrogens and `ligand_mol` builds the molecule from
    that same PDB block, so an RDKit atom index is an index into the heavy-atom
    list the contact search used.
    """
    from rdkit import Chem
    if mol is None or not warhead_class:
        return None
    fs = sorted(glob.glob(str(REPO / "data" / "reference" / "warhead_classes_*.csv")))
    if not fs:
        return None
    d = pd.read_csv(fs[-1])
    hit = d[d.class_id == warhead_class]
    if hit.empty:
        return None
    sma = str(hit.iloc[0].reactive_atom_smarts or "")
    patt = Chem.MolFromSmarts(sma) if sma else None
    if patt is None:
        return None
    m = mol.GetSubstructMatches(patt)
    return int(m[0][0]) if m else None


def interaction_3d(movie_pdb: Path, rows: list, elem_id: str,
                   cys_resi: int = 63, offset: int = 50,
                   rx_atom: int | None = None) -> str:
    """The same contacts, drawn on the real 3D pose.

    Everything here is measured: the residues sit where they sit, and each dashed
    line joins the actual pair of atoms -- one ligand, one protein -- that are
    closest in this frame. The 2D map projects; this does not.
    """
    import numpy as np
    frame, nfr = representative_frame(movie_pdb)
    if frame is None:
        return ""
    keep = [r for r in rows if r[2] >= 0.20 and r[4] is not None][:12]
    if not keep:
        return ""

    lig, prot = [], []
    for l in frame.splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[:1])
        if el == "H":
            continue
        rec = ((float(l[30:38]), float(l[38:46]), float(l[46:54])),
               l[17:20].strip(), l[22:26].strip(), el)
        (lig if rec[1] == "MOL" else prot).append(rec)
    if not lig or not prot:
        return ""

    # the measured pair per residue: its atom nearest that residue's contact atom
    links, resis = [], []
    for rn, ri, frac, pol, ai in keep:
        if ai >= len(lig):
            continue
        is_cys = (rn == "CYS" and str(ri) == str(cys_resi))
        # CYS113 IS DRAWN AS THE REACTION, NOT AS A CONTACT. For every other
        # residue the pair is "closest atom to closest atom", which is the right
        # question. For the catalytic cysteine the right question is the attack
        # vector: S(gamma) to the electrophilic carbon. Left to the generic rule
        # this molecule drew Cys113 to its bromide, which is the leaving group.
        pick = rx_atom if (is_cys and rx_atom is not None
                           and rx_atom < len(lig)) else ai
        la = np.array(lig[pick][0])
        cand = [(np.linalg.norm(np.array(p[0]) - la), p) for p in prot
                if p[1] == rn and p[2] == ri
                and (not is_cys or p[3] == "S" or rx_atom is None)]
        if not cand:
            continue
        dist, pa = min(cand, key=lambda t: t[0])
        try:
            shown = f"{rn}{int(ri) + offset}"
        except ValueError:
            shown = f"{rn}{ri}"
        col = "0xb3261e" if is_cys else ("0x0f7a54" if pol > 0.2 else "0x4a6885")
        tag = (f"{shown} Sγ→C  {dist:.1f} A" if is_cys and rx_atom is not None
               else f"{shown}  {dist:.1f} A")
        links.append({"a": list(la), "b": list(pa[0]), "c": col,
                      "t": tag, "d": pol > 0.2})
        resis.append(int(ri))

    # THE POCKET WALL. Whole-protein VDW is the expensive call in a 3Dmol viewer
    # and nobody looks at the far side of the protein; the residues with a heavy
    # atom within SURF_SHELL_A of the ligand are what actually forms the pocket.
    LX = np.array([a[0] for a in lig])
    pocket = sorted({int(p[2]) for p in prot if p[2].lstrip("-").isdigit()
                     and np.linalg.norm(LX - np.array(p[0]), axis=1).min()
                     <= SURF_SHELL_A})

    pdb = "\n".join(l for l in frame.splitlines()
                     if l.startswith(("ATOM", "HETATM")))
    return f"""
<div class="glwrap"><div class="glbox">
<div id="{elem_id}" style="position:absolute;inset:0"></div></div>
<div class="vctl"><label class="sfx"><input type="checkbox" id="{elem_id}-surf"
 checked> pocket surface</label>
<span class="hint">{len(pocket)} residues within {SURF_SHELL_A:.0f} &#8491; of the
ligand &#183; Cys113 and the ligand are left uncovered</span></div>
{_LEGEND_3D}
<p class="p3cap">Representative frame — the ligand's medoid position over
{nfr} frames. Every line joins the actual closest pair of atoms in this frame,
with the distance; residues are where they really are. Cys113 is drawn as the
attack vector — S&gamma; to the electrophilic carbon — not as its nearest
contact, which for a halide-displacement warhead is the leaving group.</p></div>
<script type="text/plain" id="{elem_id}-pdb">{pdb}</script>
<script>
(function(){{
  const M = window.$3Dmol || window['3Dmol'];
  const L = {json.dumps(links)}, RES = {json.dumps(sorted(set(resis)))};
  const POCKET = {json.dumps(pocket)}, CYS = {int(cys_resi)};
  // Same lazy boot as the movie: a closed <details> has no height, and a viewer
  // built into a zero-height box draws nothing.
  let built = false;
  function boot(){{
    if (built) return; built = true;
    requestAnimationFrame(function(){{ requestAnimationFrame(function(){{
      const v = M.createViewer(document.getElementById('{elem_id}'),
                               {{backgroundColor:'#eef1f6'}});
      v.addModel(document.getElementById('{elem_id}-pdb').textContent, 'pdb');
      v.setStyle({{}}, {{cartoon:{{color:'#c3ccd8', opacity:0.55}}}});
      // Contact residues are NEUTRAL. They were greenCarbon, and green is the
      // legend's word for "polar" -- a residue drawn green for being a contact
      // and a line drawn green for being polar cannot share a page.
      v.setStyle({{resi: RES}}, {{stick:{{radius:0.14, colorscheme:'whiteCarbon'}},
                                 cartoon:{{color:'#c3ccd8', opacity:0.55}}}});
      // Cys113 in its own sticks, carbons in the key's red, sulfur left at its
      // element colour so the atom under attack is the one you can pick out.
      const CC = Object.assign({{}}, (M.elementColors || {{}}).defaultColors || {{}},
                               {{C: 0xb3261e}});
      v.setStyle({{resi: [CYS]}},
                 {{stick:{{radius:0.26, colorscheme:{{prop:'elem', map: CC}}}},
                  cartoon:{{color:'#c3ccd8', opacity:0.55}}}});
      v.setStyle({{resn:'MOL'}}, {{stick:{{radius:0.22, colorscheme:'yellowCarbon'}}}});
      L.forEach(function(k){{
        v.addCylinder({{start:{{x:k.a[0],y:k.a[1],z:k.a[2]}},
                       end:{{x:k.b[0],y:k.b[1],z:k.b[2]}},
                       radius:0.045, color:k.c, dashed:k.d, fromCap:1, toCap:1}});
        v.addLabel(k.t, {{position:{{x:(k.a[0]+k.b[0])/2, y:(k.a[1]+k.b[1])/2,
                                    z:(k.a[2]+k.b[2])/2}},
                         fontSize:10, fontColor:k.c, backgroundColor:'white',
                         backgroundOpacity:0.72, borderThickness:0}});
      }});
      // THE POCKET MESH. Built over the shell only, and never over the ligand or
      // Cys113 -- a surface drawn on top of them hides the two things the figure
      // exists to show. Translucent, so the sticks and the measured lines read
      // through it. removeSurface before re-adding: addSurface STACKS meshes.
      let surf = null;
      const chk = document.getElementById('{elem_id}-surf');
      function setSurf() {{
        if (surf) {{ try {{ v.removeSurface(surf.surfid); }} catch (e) {{}} surf = null; }}
        if (!chk || chk.checked) {{
          surf = v.addSurface(M.SurfaceType.VDW,
            {{opacity: 0.62, color: '#b9c7db'}},
            {{and: [{{resi: POCKET}},
                   {{not: {{or: [{{resn:'MOL'}}, {{resi: [CYS]}}]}}}}]}});
        }}
        v.render();
      }}
      if (chk) chk.addEventListener('change', setSurf);
      setSurf();
      v.zoomTo({{resn:'MOL'}}); v.zoom(0.55); v.resize();
      // 3Dmol draws NOTHING until render() is called. The labels are DOM
      // overlays and appear without it, which is what made an unrendered
      // viewer look like a viewer with a missing molecule.
      v.render();
    }}); }});
  }}
  (function(){{
    const host = document.getElementById('{elem_id}');
    const det = host && host.closest ? host.closest('details') : null;
    if (det) {{
      if (det.open) window.addEventListener('load', boot);
      det.addEventListener('toggle', function(){{ if (det.open) boot(); }});
    }} else {{ window.addEventListener('load', boot); }}
  }})();
}})();
</script>"""


def contact_distances(movie_pdb: Path, rows: list) -> dict:
    """Closest measured atom-pair distance per residue, in the medoid frame.

    The same frame and the same pairs the 3D view draws, so the number beside a
    residue in the flat map is the number on the line in the 3D one. Two figures
    of the same contact disagreeing by a tenth of an Angstrom is the kind of
    thing that costs a reader an hour.
    """
    import numpy as np
    frame, _ = representative_frame(movie_pdb)
    if frame is None:
        return {}
    lig, prot = [], []
    for l in frame.splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if (l[76:78].strip() or l[12:16].strip()[:1]) == "H":
            continue
        rec = ((float(l[30:38]), float(l[38:46]), float(l[46:54])),
               l[17:20].strip(), l[22:26].strip())
        (lig if rec[1] == "MOL" else prot).append(rec)
    out = {}
    for rn, ri, _f, _p, ai in rows:
        if ai is None or ai >= len(lig):
            continue
        la = np.array(lig[ai][0])
        cand = [np.linalg.norm(np.array(p[0]) - la) for p in prot
                if p[1] == rn and p[2] == ri]
        if cand:
            out[(rn, ri)] = float(min(cand))
    return out


def interaction_map(mol, rows: list, n_frames: int, cys_resi: int = 63,
                    offset: int = 50, dist_of: dict | None = None,
                    rx_atom: int | None = None) -> str:
    """A real interaction diagram: each residue drawn against the atom it contacts.

    The ligand is rendered by RDKit, and RDKit is then asked where it PUT each
    atom (`GetDrawCoords`). Residues are placed on that same canvas, along the
    ray from the molecule's centre through their contact atom, and joined to that
    atom by a line. The previous version arranged residues on a ring in arbitrary
    order and joined nothing to anything -- it looked like an interaction map and
    carried none of the information one has.

    Still not a LigPlot in one respect, stated on the figure: the residue's
    position is a projection along that ray, not its real 3D position. The ATOM
    it is joined to is measured.
    """
    import math
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Draw, rdCoordGen
    RDLogger.DisableLog("rdApp.*")
    if mol is None:
        return ""
    dist_of = dist_of or {}
    keep = [r for r in rows if r[2] >= 0.20 and r[4] is not None][:12]
    if not keep:
        return ""

    # THE MOLECULE GETS THE MIDDLE, NOT THE WHOLE CANVAS. Drawn at full size it
    # filled the frame and the residues had nowhere to go but on top of it, which
    # is what made the last version unreadable and gave no sense of separation.
    # It is rendered into an inner box and placed centred, leaving an annulus.
    W, H = 980, 720
    IW, IH = 470, 340
    OX, OY = (W - IW) / 2, (H - IH) / 2
    m = Chem.Mol(mol)
    m.RemoveAllConformers()
    rdCoordGen.AddCoords(m)
    d = Draw.rdMolDraw2D.MolDraw2DSVG(IW, IH)
    o = d.drawOptions()
    o.bondLineWidth = 2
    o.additionalAtomLabelPadding = 0.15
    # the ligand occupies the middle; the ring outside it is for residues
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    inner = svg[svg.index(">", svg.index("<svg")) + 1: svg.rindex("</svg>")]
    inner = re.sub(r"<rect[^>]*>", "", inner, count=1)     # drop its white backdrop
    inner = f"<g transform='translate({OX:.0f},{OY:.0f})'>{inner}</g>"

    pos = {}
    for i in range(m.GetNumAtoms()):
        pt = d.GetDrawCoords(i)
        pos[i] = (pt.x + OX, pt.y + OY)     # same shift as the drawing above
    cx = sum(p[0] for p in pos.values()) / len(pos)
    cy = sum(p[1] for p in pos.values()) / len(pos)

    # EVERY BOX IS SIZED TO ITS OWN TEXT. A fixed 104px rect fitted "ALA124" and
    # not "3.2 A - 78% of frames", so the second line ran out through the border
    # and over whatever was behind it -- which read as labels being covered.
    # Collision is tested rectangle-against-rectangle for the same reason: a
    # single radius cannot describe boxes of different widths.
    def _label(rn, ri, frac, dist):
        try:
            shown = f"{rn}{int(ri) + offset}"
        except ValueError:
            shown = f"{rn}{ri}"
        sub = (f"{dist:.1f} &#8491; &#183; {frac*100:.0f}%" if dist
               else f"{frac*100:.0f}% of frames")
        plain = sub.replace("&#8491;", "A").replace("&#183;", "-")
        w = max(len(shown) * 7.4, len(plain) * 5.9) + 22
        return shown, sub, max(w, 76.0)

    placed = []
    for rn, ri, frac, pol, ai in sorted(keep, key=lambda r: -r[2]):
        is_cys = (rn == "CYS" and str(ri) == str(cys_resi))
        # Cys113 points at the electrophilic carbon, not at whichever atom
        # happens to be nearest -- see reactive_atom_index.
        pick = rx_atom if (is_cys and rx_atom is not None) else ai
        ax, ay = pos.get(pick, pos.get(ai, (cx, cy)))
        shown, sub, bw = _label(rn, ri, frac, dist_of.get((rn, ri)))
        bh = 32.0
        vx, vy = ax - cx, ay - cy
        L = math.hypot(vx, vy) or 1.0
        vx, vy = vx / L, vy / L
        r = 150.0
        for _ in range(40):
            x, y = ax + vx * r, ay + vy * r
            x = min(max(x, bw / 2 + 8), W - bw / 2 - 8)
            y = min(max(y, bh / 2 + 6), H - bh / 2 - 8)
            if all(abs(x - px) > (bw + pw) / 2 + 8 or abs(y - py) > bh + 5
                   for px, py, pw, *_ in placed):
                break
            r += 22.0
        placed.append((x, y, bw, ax, ay, rn, ri, frac, pol, shown, sub, is_cys))

    parts = []
    for x, y, bw, ax, ay, rn, ri, frac, pol, shown, sub, is_cys in placed:
        col = "#b3261e" if is_cys else ("#0f7a54" if pol > 0.2 else "#4a6885")
        dash = " stroke-dasharray='5 4'" if pol > 0.2 else ""
        parts.append(
            f"<line x1='{ax:.0f}' y1='{ay:.0f}' x2='{x:.0f}' y2='{y:.0f}' "
            f"stroke='{col}' stroke-width='{0.9 + 2.4 * frac:.1f}' "
            f"opacity='.55'{dash}/>"
            f"<circle cx='{ax:.0f}' cy='{ay:.0f}' r='3.4' fill='{col}'/>"
            f"<rect x='{x - bw/2:.0f}' y='{y-16:.0f}' width='{bw:.0f}' "
            f"height='32' rx='6' fill='#ffffff' fill-opacity='.94' "
            f"stroke='{col}' stroke-width='1.2'/>"
            f"<text x='{x:.0f}' y='{y-3:.0f}' class='rl' fill='{col}'>"
            f"{shown}{' S&#947;&#8594;C' if is_cys and rx_atom is not None else ''}"
            f"</text>"
            f"<text x='{x:.0f}' y='{y+10:.0f}' class='rf' fill='{col}'>"
            f"{sub}</text>")

    LEGY = H + 16
    return f"""<svg viewBox="0 0 {W} {H + 98}" class="imap" role="img"
 aria-label="each contacting residue joined to the ligand atom it contacts">
<style>.rl{{font:600 11.5px ui-monospace,monospace;text-anchor:middle}}
.rf{{font:8.5px ui-monospace,monospace;text-anchor:middle;opacity:.85}}
.ik{{font:10.5px Helvetica,Arial,sans-serif;fill:#5b6b80}}
.kh{{font:600 10.5px Helvetica,Arial,sans-serif;fill:#3c4a5c;
  letter-spacing:.06em;text-transform:uppercase}}
.kt{{font:600 10.5px ui-monospace,monospace}}
.kw{{font:10.5px Helvetica,Arial,sans-serif;fill:#5b6b80}}</style>
{inner}
{''.join(parts)}
{_legend_svg(12, LEGY)}
<text x="470" y="{LEGY + 17}" class="ik">Line width = fraction of the
{n_frames} frames in contact. Crystal numbering.</text>
<text x="470" y="{LEGY + 32}" class="ik">Each residue is joined to the ligand ATOM
it is nearest to in the most frames &#8212; that atom is measured.</text>
<text x="470" y="{LEGY + 47}" class="ik">The residue's POSITION on the page is
projected along that direction, not its real 3D position.</text>
</svg>"""


def md_row(ident: str):
    parts = []
    for f in glob.glob(str(rp.residence_dir() / "*.csv")):   # this run only (#74)
        try:
            parts.append(pd.read_csv(f))
        except Exception:                                # noqa: BLE001
            continue
    if not parts:
        return None
    d = pd.concat(parts, ignore_index=True)
    d = d[(d.ident.astype(str) == ident) & (d.get("production_ps", 0) >= 50000)]
    if "status" in d.columns:
        d = d[d.status.astype(str).str.startswith("ok")]
    e = "explicit_frac_frames_engaged"
    if e in d.columns:
        d = d[d[e].notna()]
    return None if d.empty else d.iloc[-1]


def sweep_row(ident: str):
    parts = []
    for f in sorted(glob.glob(str(rp.sweep_dir() / "attack_sweep_*.csv"))):  # (#74)
        try:
            parts.append(pd.read_csv(f))
        except Exception:                                # noqa: BLE001
            continue
    if not parts:
        return None
    d = pd.concat(parts, ignore_index=True)
    d = d[(d.parent_ident.astype(str) == ident) & (d.status == "ok")]
    return None if d.empty else d.sort_values("frac_attack_ready").iloc[-1]


def classes() -> dict:
    out = {}
    for tier, score in (("T4", "conditional_eb"), ("T3", "enrichment_conditional")):
        fs = sorted(glob.glob(str(B / f"rank_v2/rank_v2_{tier}_{score}_*.csv")))
        if fs:
            d = pd.read_csv(fs[-1]).drop_duplicates("parent_ident")
            out.update(dict(zip(d.parent_ident, d.warhead_class)))
    return out


def block(ident: str, er, three: str, cls: dict,
          rep_dir: Path | None = None, suffix: str = "",
          heading: str | None = None) -> str:
    """One molecule's section. `rep_dir` overrides the default replicate-1 run,
    and `suffix` keeps element ids unique when the same molecule appears more
    than once in a document -- getElementById returns the FIRST match, so two
    sections sharing ids means every control drives the first viewer."""
    rep = Path(rep_dir) if rep_dir else MD / ident / "md" / "rep1"
    if not rep.is_dir():
        log.warning("%s: no trajectory at %s", ident, rep)
        return ""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "mdprio_report", REPO / "scripts" / "mdprio_report.py")
    mp = importlib.util.module_from_spec(spec)
    sys.modules["mdprio_report"] = mp
    spec.loader.exec_module(mp)

    total_ns = mp.prod_ns(rep)
    s = mp.series(rep, er, total_ns)
    res = mp.residence(s)
    if res["status"] != "ok":
        log.warning("%s: %s", ident, res["status"])
        return ""

    mpdb = rep / "movie.pdb"
    if not mpdb.is_file():
        mov.build_movie_pdb(rep, mpdb, total_ps=total_ns * 1000.0)
    movie = ""
    nacs = None
    if mpdb.is_file():
        pdb_txt, dsg, labels, lpos = er.surface_payload(mpdb)
        movie = mov.viewer_html(pdb_txt, dsg, labels, lpos, "", elem_id=f"gl_{ident}{suffix}")
        nacs = mp.nac_series(ident, rep, mpdb, total_ns)
    img = mp.figure(ident, s, res, er, nacs)

    # THE STRUCTURE THE CHEMIST CAN OPEN. First model of the fitted movie: the
    # protein and ligand as simulated, PBC-repaired and CA-fitted, one frame.
    # Offered as a download rather than a path, because the recipient has no
    # access to this filesystem.
    i3d, pdb_href, pdb_bytes = "", "", 0
    if mpdb.is_file():
        raw = mpdb.read_text()
        first = raw.split("ENDMDL")[0]
        if not first.lstrip().startswith(("MODEL", "ATOM", "HETATM", "TITLE", "REMARK")):
            first = raw
        frame1 = first.replace("MODEL", "REMARK MODEL", 1).rstrip() + "\nEND\n"
        pdb_bytes = len(frame1.encode())
        pdb_href = ("data:chemical/x-pdb;base64,"
                    + base64.b64encode(frame1.encode()).decode())
        try:
            rows_c, nfr = contacts(mpdb)
            lm = ligand_mol(mpdb, smiles_of(ident) or "")
            dmap = contact_distances(mpdb, rows_c)
            rx = reactive_atom_index(lm, cls.get(ident, ""))
            if rx is None:
                log.warning("%s: no reactive atom for class %r; Cys113 will be "
                            "drawn to its nearest atom", ident, cls.get(ident))
            # The flat map is no longer emitted (@tt8804): the 3D view carries the
            # same contacts in the real geometry, and two figures of one thing
            # invite a reader to look for a difference that is only projection.
            # interaction_map() is kept -- it is the printable version.
            i3d = interaction_3d(mpdb, rows_c, f"i3_{ident}{suffix}", rx_atom=rx)
        except Exception as exc:                          # noqa: BLE001
            log.warning("%s: interaction map unavailable: %s", ident, exc)

    smi = smiles_of(ident) or ""
    svg = depiction(smi) if smi else ""
    m, sw = md_row(ident), sweep_row(ident)

    rows = [("trajectory", f"{res['length_ns']:.1f} ns, {res['n_frames']:,} frames"),
            ("warhead class", cls.get(ident, "unclassified")),
            ("mean ligand RMSD", f"{res['rmsd_mean_nm']:.3f} nm"),
            ("max ligand RMSD", f"{res['rmsd_max_nm']:.3f} nm"),
            ("final ligand RMSD", f"{res['rmsd_final_nm']:.3f} nm"),
            ("residence fraction", f"{res['residence_frac']:.3f}")]
    if m is not None and pd.notna(m.get("explicit_frac_frames_engaged")):
        rows.insert(2, ("target engagement, 100 ns",
                        f"{float(m['explicit_frac_frames_engaged'])*100:.2f}%"))
    if res.get("left_at_ns") is not None:
        rows.append(("left the pocket at", f"{res['left_at_ns']:.1f} ns"))
    if sw is not None:
        rows.append((f"attack-ready, {gs.sweep_label()} sweep",
                     f"{float(sw.frac_attack_ready):.4f}"))
        rows.append(("median C&ndash;S&gamma; distance",
                     f"{float(sw.median_dist_a):.2f} &Aring;"))
        rows.append(("median attack angle", f"{float(sw.median_angle_deg):.1f}&deg;"))
    facts = "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in rows)
    dl = (f'<a class="dl" download="{ident}{suffix}_md.pdb" href="{pdb_href}">'
          f'Download the MD structure (PDB, {pdb_bytes/1024:.0f} KB)</a>'
          if pdb_href else "")

    return f"""
<section class="mol">
  <h2>{heading or ident}</h2>
  <div class="top">
    <div class="struct">{f'<img alt="" src="{svg}">' if svg else ''}</div>
    <div class="side">
      <label for="s_{ident}{suffix}">SMILES</label>
      <textarea id="s_{ident}{suffix}" readonly rows="3" onclick="this.select()">{smi}</textarea>
      <table class="kv">{facts}</table>
      {dl}
    </div>
  </div>
  <details class="panel"><summary>RMSD plots
    <span class="hint">ligand RMSD, warhead&ndash;Cys113 distance, attack angle</span></summary>
    <div class="pbody">
      <img class="plots" src="data:image/png;base64,{img}" alt="RMSD, distance and angle traces">
    </div></details>
  <details class="panel"><summary>MD movie
    <span class="hint">{res['length_ns']:.0f} ns, surface by charge, ligand in yellow</span></summary>
    <div class="pbody">{movie}</div></details>
  <details class="panel"><summary>Interactions in the 3D pose
    <span class="hint">real positions, each line a measured atom pair</span></summary>
    <div class="pbody">{i3d}</div></details>
</section>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--name", default="Shortlist")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    # Version from the CHANGELOG, via the GUI's own parser rather than a second
    # copy of the logic -- one source, so the report and the GUI cannot disagree
    # about which release produced the numbers.
    import importlib.util as _u
    _sp = _u.spec_from_file_location("mdprio_combine", REPO / "scripts" / "mdprio_combine.py")
    _mc = _u.module_from_spec(_sp); _sp.loader.exec_module(_mc)
    ver, code = _mc._version()

    er = _er()
    three = (REPO / "scripts/.cache_3dmol-min.js").read_text()
    cls = classes()
    blocks = [block(c, er, three, cls) for c in args.candidates]
    blocks = [b for b in blocks if b]
    if not blocks:
        raise SystemExit("nothing to report")

    title = f"{date.today().isoformat()} {args.name}"
    byline_ver = " ".join(x for x in (f"version {ver}" if ver else "",
                                      f"\u201c{code}\u201d" if code else "") if x)
    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{rt.CSS}{mov.VIEWER_CSS}
/* CENTRED ON THE BODY, not on each child. rt.CSS centres direct children
   individually, which leaves a page of differently-sized blocks looking
   left-anchored; giving the body itself the measure centres the whole column and
   every child then fills it. */
body{{max-width:1180px;margin:0 auto;padding:0 30px 70px}}
body>*{{max-width:none;padding-left:0;padding-right:0}}
section.mol{{border-top:4px solid var(--rule);padding-top:1.6rem;margin-top:2rem}}
section.mol:first-of-type{{border-top:0;margin-top:1rem}}
header.mast{{border-bottom:4px solid var(--rule);padding-bottom:1rem}}
h2{{font-family:var(--mono);font-size:1.05rem;color:var(--navy);margin:0 0 .8rem}}
h3{{font-size:.66rem;font-family:var(--mono);letter-spacing:.14em;
  text-transform:uppercase;color:var(--blue);margin:1.4rem 0 .5rem}}
.top{{display:grid;grid-template-columns:360px 1fr;gap:22px;align-items:start}}
@media(max-width:820px){{.top{{grid-template-columns:1fr}}}}
.struct img{{width:100%;height:auto;background:#fff;border:1px solid var(--rule);
  border-radius:5px}}
.side label{{font:600 .6rem var(--mono);letter-spacing:.12em;text-transform:uppercase;
  color:var(--muted);display:block;margin-bottom:.3rem}}
textarea{{width:100%;font-family:var(--mono);font-size:12.5px;padding:.5rem;
  border:1px solid var(--rule);border-radius:4px;background:var(--card);
  color:var(--ink);resize:vertical;margin-bottom:.9rem}}
table.kv{{border-collapse:collapse;width:100%;font-size:13px}}
table.kv th{{text-align:left;font-weight:500;color:var(--muted);padding:3px 14px 3px 0;
  white-space:nowrap}}
table.kv td{{font-family:var(--mono);padding:3px 0}}
img.plots{{width:100%;height:auto;border:1px solid var(--rule);border-radius:5px;
  background:#fff}}
svg.imap{{width:100%;height:auto;background:var(--card);border:1px solid var(--rule);
  border-radius:5px}}
label.sfx{{display:flex;align-items:center;gap:.4rem;font:600 12px var(--sans);
  cursor:pointer;user-select:none}}
{KEY3_CSS}
a.dl{{display:inline-block;margin-top:.7rem;font:600 12px var(--sans);
  color:var(--blue);text-decoration:none;border:1px solid var(--blue);
  border-radius:4px;padding:.35rem .7rem}}
a.dl:hover{{background:var(--blue-pale)}}
</style></head><body>
<header class="mast"><h1>{title}</h1>
<p class="standfirst">Timothy Wu &middot; {byline_ver}</p></header>
<script>{three}</script>
{''.join(blocks)}
</body></html>"""

    dest = sout.Topic("blacksmith", "shortlist").write("shortlist", ".html")
    dest.write_text(page)
    side = B / "shortlist" / "shortlist.html"
    side.write_text(page)
    print(f"\n  {len(blocks)} molecules -> {side}  ({len(page)/1048576:.1f} MB)")
    print(f"  versioned: {dest}")


if __name__ == "__main__":
    main()
