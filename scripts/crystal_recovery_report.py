#!/usr/bin/env python3
"""
Purpose: does the screen generate the crystal pose, and does it keep it?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-11
Input: --ident (docked via dock_reference_modes) --crystal <sdf>
Output: 00_outputs/blacksmith/crystal_recovery/crystal_recovery_<ident>_<N>.html

TWO SEPARATE QUESTIONS, AND THEY HAVE DIFFERENT ANSWERS (@tt8804, #39/#56).

  1. Does docking SAMPLE the crystal pose?   -- best RMSD over the whole cloud
  2. Does the pipeline KEEP it?              -- RMSD of the exported representative

A screen can pass 1 and fail 2, and if it does, no amount of downstream
simulation recovers it: everything after the export sees one pose per mode.

THE COMPARISON IS INTERNALLY CONTROLLED. The crystal ligand is the covalent
ADDUCT -- chlorine gone, 16 heavy atoms against the reactant's 17 -- so an
absolute RMSD to it carries a component from bond formation pulling the ligand
in, and no docked pose can score zero. That penalty applies equally to every pose
in the cloud, so the statement this report leads with is a RELATIVE one: where
the kept pose sits in the distribution of poses that were available. That
quantity is unaffected by the adduct/reactant difference.

Atoms are matched by maximum common substructure, computed once and reused for
every pose, so the same atoms are compared throughout.
"""

from __future__ import annotations

import argparse
import base64
import glob
import html
import io
import json
import logging
import sys
from datetime import date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                         # noqa: E402
from shared import report_theme as rt                      # noqa: E402

log = logging.getLogger("xtal-recovery")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
RECEPTOR = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb")
CYS_RESI = 113


def load(p: Path):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    return [m for m in Chem.SDMolSupplier(str(p), removeHs=True) if m is not None]


def rmsd_series(xtal, poses):
    """In-place RMSD of every pose to the crystal, over the shared MCS.

    NO SUPERPOSITION. Both are already in the receptor's frame, and fitting them
    onto each other would answer "is it the same shape" rather than "is it in the
    same place", which is the question.
    """
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    mcs = rdFMCS.FindMCS([xtal, poses[0]], timeout=30, ringMatchesRingOnly=True)
    patt = Chem.MolFromSmarts(mcs.smartsString)
    mx = xtal.GetSubstructMatch(patt)
    cx = xtal.GetConformer()
    A = np.array([list(cx.GetAtomPosition(j)) for j in mx])
    out = []
    for p in poses:
        mp = p.GetSubstructMatch(patt)
        if not mp or len(mp) != len(mx):
            out.append(np.nan); continue
        cp = p.GetConformer()
        Bq = np.array([list(cp.GetAtomPosition(j)) for j in mp])
        out.append(float(np.sqrt(((A - Bq) ** 2).sum(axis=1).mean())))
    return np.array(out), mcs.numAtoms


def hist_png(vals, kept, best):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.6, 2.9), dpi=150)
    ax.hist(vals[~np.isnan(vals)], bins=40, color="#9fb4cc", edgecolor="white")
    ax.axvline(2.0, color="#5b6b80", ls=":", lw=1.4)
    ax.text(2.03, ax.get_ylim()[1]*0.92, "2 Å", color="#5b6b80", fontsize=8)
    ax.axvline(best, color="#0f7a54", lw=2)
    ax.text(best, ax.get_ylim()[1]*0.72, f"  best {best:.2f}", color="#0f7a54", fontsize=8)
    ax.axvline(kept, color="#b3261e", lw=2)
    ax.text(kept, ax.get_ylim()[1]*0.52, f"  kept {kept:.2f}", color="#b3261e", fontsize=8)
    ax.set_xlabel("in-place RMSD to the crystal ligand (Å), over the shared MCS")
    ax.set_ylabel("docked poses")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png"); plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def pdb_block(mol, resname="MOL"):
    from rdkit import Chem
    m = Chem.Mol(mol)
    for a in m.GetAtoms():
        ri = a.GetPDBResidueInfo()
        if ri is None:
            ri = Chem.AtomPDBResidueInfo(f" {a.GetSymbol():<3s}")
            a.SetMonomerInfo(ri)
        ri.SetResidueName(resname); ri.SetResidueNumber(1); ri.SetIsHeteroAtom(True)
    return "\n".join(l for l in Chem.MolToPDBBlock(m).splitlines()
                     if l.startswith(("ATOM", "HETATM")))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--ident", required=True)
    ap.add_argument("--crystal", required=True)
    ap.add_argument("--all-poses", default=None)
    ap.add_argument("--representative", default=None)
    ap.add_argument("--name", default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    allp = Path(args.all_poses or B / f"nac_v3_allposes/{args.ident}.sdf")
    repp = Path(args.representative or B / f"nac_v3_poses/{args.ident}.sdf")
    xtal = load(Path(args.crystal))[0]
    poses, reps = load(allp), load(repp)
    log.info("%d poses, %d representatives, crystal %d heavy atoms",
             len(poses), len(reps), xtal.GetNumAtoms())

    vals, n_mcs = rmsd_series(xtal, poses)
    rvals, _ = rmsd_series(xtal, reps)
    ok = vals[~np.isnan(vals)]
    kept = float(np.nanmin(rvals))
    best = float(ok.min())
    pct_better = float((ok < kept).mean() * 100)

    three = (REPO / "scripts/.cache_3dmol-min.js").read_text()
    rec = "\n".join(l for l in RECEPTOR.read_text().splitlines()
                    if l.startswith(("ATOM", "HETATM")))
    # The three structures the argument is about: the crystal, the pose the
    # pipeline kept, and the closest pose it generated and did not keep.
    # The pocket wall for the surface: residues near Cys113's SG. A whole-protein
    # VDW mesh is the expensive call and the far side of the protein is not the
    # subject.
    import numpy as _np
    sg, res = None, {}
    for ln in RECEPTOR.read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        ri, nm = ln[22:26].strip(), ln[12:16].strip()
        try:
            xyz = _np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
        res.setdefault(ri, []).append(xyz)
        if ri == str(CYS_RESI) and nm == "SG":
            sg = xyz
    pocket = sorted(int(r) for r, xs in res.items() if r.lstrip("-").isdigit()
                    and sg is not None
                    and min(float(_np.linalg.norm(x - sg)) for x in xs) <= 12.0)

    closest = poses[int(np.nanargmin(vals))]
    blocks = {"xtal": pdb_block(xtal), "kept": pdb_block(reps[int(np.nanargmin(rvals))]),
              "closest": pdb_block(closest)}

    ver = ""
    try:
        import importlib.util as _u
        sp = _u.spec_from_file_location("mc", REPO / "scripts/mdprio_combine.py")
        mc = _u.module_from_spec(sp); sp.loader.exec_module(mc)
        v, c = mc._version(); ver = f"version {v} “{c}”"
    except Exception:                                      # noqa: BLE001
        pass

    title = f"{date.today().isoformat()} {args.name or args.ident} — crystal pose recovery"
    page = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<script>{three}</script>
<style>{rt.CSS}
body{{max-width:1120px;margin:0 auto;padding:0 30px 70px}}
body>*{{max-width:none}}
.big{{display:flex;gap:2.2rem;flex-wrap:wrap;margin:1.2rem 0 1.6rem}}
.big div b{{display:block;font:600 1.9rem var(--mono);line-height:1.1}}
.big div span{{font-size:12px;color:var(--muted)}}
.glbox{{position:relative;width:100%;height:480px;background:#eef1f6;
 border:1px solid var(--rule);border-radius:5px;overflow:hidden}}
.glbox>div{{position:absolute;inset:0}} .glbox canvas{{position:absolute;top:0;left:0}}
.key{{display:flex;gap:1.2rem;flex-wrap:wrap;margin:.7rem 0 0;font:600 12px var(--sans)}}
.key label{{display:flex;align-items:center;gap:.35rem;cursor:pointer;user-select:none}}
.key i{{width:22px;height:4px;border-radius:2px;display:inline-block;margin-right:.4rem;
 vertical-align:3px}}
img.h{{width:100%;height:auto;border:1px solid var(--rule);border-radius:5px;background:#fff}}
</style></head><body>
<header class="mast"><h1>{html.escape(title)}</h1>
<p class="standfirst">Timothy Wu &middot; {ver}</p></header>

<p><strong>Two questions, and they have different answers.</strong> Does docking
sample the crystal pose, and does the pipeline keep it?</p>

<div class="big">
 <div><b>{best:.2f} Å</b><span>best pose generated, of {len(ok)}</span></div>
 <div><b>{int((ok<=2.0).sum())}</b><span>poses within 2 Å ({(ok<=2.0).mean()*100:.0f}%)</span></div>
 <div><b style="color:var(--bad)">{kept:.2f} Å</b><span>the pose the pipeline KEPT</span></div>
 <div><b style="color:var(--bad)">{pct_better:.0f}%</b><span>of the cloud is closer than the kept pose</span></div>
</div>

<p>The screen <strong>does</strong> generate the crystal pose &mdash; the closest
of {len(ok)} docked poses sits {best:.2f} Å from it, and
{int((ok<=2.0).sum())} are inside the 2 Å bar normally used to call a pose
correct. The pose the pipeline exported as this mode&rsquo;s representative sits
at <strong>{kept:.2f} Å</strong>, which is worse than
<strong>{pct_better:.0f}%</strong> of the poses it had to choose from.</p>

<img class="h" alt="RMSD distribution" src="data:image/png;base64,{hist_png(vals, kept, best)}">

<h2>The three structures</h2>
<div class="glbox"><div id="gl"></div></div>
<div class="key">
 <label><input type="checkbox" id="c-xtal" checked>
  <i style="background:#b3261e"></i>crystal (6VAJ, adduct)</label>
 <label><input type="checkbox" id="c-closest" checked>
  <i style="background:#0f7a54"></i>closest generated pose ({best:.2f} Å)</label>
 <label><input type="checkbox" id="c-kept" checked>
  <i style="background:#8a6d1f"></i>pose the pipeline kept ({kept:.2f} Å)</label>
 <label><input type="checkbox" id="c-surf" checked>
  <i style="background:#b9c7db"></i>pocket surface</label>
 <label><input type="checkbox" id="c-cys" checked>
  <i style="background:#e8c33a"></i>Cys113 sticks</label>
</div>

<p class="note" style="margin-top:1.4rem">Atoms are matched by maximum common
substructure &mdash; <strong>{n_mcs} atoms</strong> &mdash; computed once and
reused for every pose. <strong>No superposition:</strong> both are already in the
receptor&rsquo;s frame, and fitting them onto each other would ask whether they
are the same shape rather than whether they are in the same place.</p>

<p class="note"><strong>What the absolute number is worth.</strong> The crystal
ligand is the covalent <em>adduct</em> &mdash; chlorine gone, {xtal.GetNumAtoms()}
heavy atoms against the reactant&rsquo;s {poses[0].GetNumAtoms()} &mdash; so bond
formation has pulled it toward Cys113 and no docked pose can score zero. That
penalty applies equally to every pose in the cloud, which is why the claim above
is a <em>relative</em> one: where the kept pose sits among the poses that were
available. That comparison is internally controlled and unaffected by the
adduct/reactant difference. <strong>n = 1 molecule.</strong></p>

<script type="text/plain" id="rec">{rec}</script>
<script>
// EVERY ELEMENT THIS READS IS DECLARED ABOVE IT. The receptor block used to sit
// after this script, so getElementById('rec') returned null, the block threw on
// its first line, and the box rendered empty with no error a reader could see --
// the fourth distinct cause of a blank viewer in this project in one day.
// tests/test_viewer_html.py now asserts this ordering.
window.addEventListener('DOMContentLoaded', function(){{
// TWO ANIMATION FRAMES BEFORE THE VIEWER IS BUILT. On DOMContentLoaded the box
// has its markup but not its final laid-out size, so 3Dmol sizes its canvas from
// stale dimensions: the scene renders into a corner of an oversized canvas and
// the mouse handlers hit-test against the wrong rectangle, which is why the view
// could not be rotated. Every other viewer in this repo already waits; this one
// did not.
requestAnimationFrame(function(){{ requestAnimationFrame(function(){{
const M = window.$3Dmol || window['3Dmol'];
const D = {json.dumps(blocks)};
const v = M.createViewer(document.getElementById('gl'), {{backgroundColor:'#eef1f6'}});
v.addModel(document.getElementById('rec').textContent, 'pdb');
// Carbons carry the identity colour; every other element keeps its conventional
// one, so the chemistry still reads while the three poses stay distinguishable.
function cs(c){{ return {{prop:'elem', map: Object.assign({{}},
  (M.elementColors||{{}}).defaultColors||{{}}, {{C:c}})}}; }}
const LIG = [['xtal',0xb3261e],['closest',0x0f7a54],['kept',0x8a6d1f]];
LIG.forEach(function(kv){{ v.addModel(D[kv[0]], 'pdb'); }});
let SURF = null;
function restyle(){{
  v.setStyle({{}}, {{cartoon:{{color:'#c3ccd8', opacity:0.5}}}});
  if (document.getElementById('c-cys').checked){{
    v.setStyle({{resi:[{CYS_RESI}]}}, {{stick:{{radius:0.26, colorscheme:'default'}},
                                       cartoon:{{color:'#c3ccd8',opacity:0.5}}}});
    v.addStyle({{resi:[{CYS_RESI}], atom:'SG'}}, {{sphere:{{radius:0.6}}}});
  }}
  LIG.forEach(function(kv, i){{
    const on = document.getElementById('c-' + kv[0]).checked;
    v.setStyle({{model: i+1}}, on
      ? {{stick:{{radius: kv[0]==='xtal' ? 0.22 : 0.17, colorscheme: cs(kv[1])}}}}
      : {{}});
  }});
  // The pocket wall, shell only and never over Cys113 or a ligand -- a mesh on
  // top of those hides what the figure is for. `and` rather than a bare `not`,
  // so the selection cannot leak onto the ligand models.
  if (SURF){{ try {{ v.removeSurface(SURF.surfid); }} catch(e){{}} SURF = null; }}
  if (document.getElementById('c-surf').checked){{
    SURF = v.addSurface(M.SurfaceType.VDW, {{opacity:0.6, color:'#b9c7db'}},
      {{and:[{{model:0}}, {{resi:{json.dumps(pocket)}}}, {{not:{{resi:[{CYS_RESI}]}}}}]}});
  }}
  v.render();
}}
['c-xtal','c-closest','c-kept','c-surf','c-cys'].forEach(function(id){{
  document.getElementById(id).addEventListener('change', restyle);
}});
restyle();
// resize BEFORE framing, so zoomTo computes against the real viewport.
v.resize();
v.zoomTo({{resn:'MOL'}});
v.zoom(0.75);
v.render();
// And re-frame if the window changes, or the canvas keeps a stale viewport.
window.addEventListener('resize', function(){{ v.resize(); v.render(); }});
}}); }});
}});
</script>
</body></html>"""

    dest = sout.Topic("blacksmith", "crystal_recovery").write(
        f"crystal_recovery_{args.ident}", ".html")
    dest.write_text(page)
    print(f"\n  best {best:.2f} A | kept {kept:.2f} A | {pct_better:.0f}% of the "
          f"cloud is closer than the kept pose\n  -> {dest}")


if __name__ == "__main__":
    main()
