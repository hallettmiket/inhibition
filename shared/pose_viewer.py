"""One pose in the pocket, drawn the same way wherever it appears.

@tt8804: "the sweep should show the structure and poses in the viewer just like
the ranked list".

WHAT IS SHARED HERE IS THE CONVENTIONS, NOT THE INTERACTION. The ranking view's
viewer compares several modes of one molecule at once -- overlay, sub-mode
colouring, per-mode checkboxes. The sweep view inspects ONE mode: the thing that
was simulated. Those are different instruments and forcing them into one widget
would make both worse.

What must NOT differ is how a pose is drawn: which residues form the pocket
shell, that Cys113 keeps conventional element colours with a sphere on its
sulfur, that carbons carry the mode colour and every other element does not,
that hydrogens are never shown. A reader moving between two pages has to be able
to trust that the same picture means the same thing. Those rules live here.

THE FOUR WAYS A 3Dmol PANEL COMES UP BLANK, all of which this file has already
cost someone: `render()` never called; a container with no height; the library
loaded after the code that uses it; the data element parsed after the code that
reads it. `mount_js` closes over none of those -- it is called explicitly, after
the DOM exists, and it measures its container first.
"""

from __future__ import annotations

#: Cartoon, pocket mesh and the anchor. Identical to the ranking view's, because
#: a surface that means "pocket" on one page and something else on another is
#: worse than no surface.
CSS = """
.pvbox{position:relative;width:100%;height:340px;background:#eef1f6;
 border:1px solid var(--rule);border-radius:4px;overflow:hidden}
:root[data-theme="dark"] .pvbox{background:#0b1016}
.pvbox>div{position:absolute;inset:0}
.pvbox canvas{position:absolute;top:0;left:0}
.pvempty{display:flex;align-items:center;justify-content:center;height:100%;
 color:var(--muted);font-size:12px;text-align:center;padding:0 2rem}
.pvstruct{width:100%;max-height:150px;object-fit:contain;background:#fff;
 border:1px solid var(--rule);border-radius:4px;margin-bottom:8px}
.pvctl{display:flex;flex-wrap:wrap;gap:.35rem 1.1rem;padding:.5rem .1rem 0;
 font-size:12px;color:var(--muted)}
.pvctl label{display:flex;align-items:center;gap:.3rem;cursor:pointer}
"""


def mount_js(cys_resi: int, pocket: str) -> str:
    """JS defining `mountPose(elemId, pdbText, mode)`.

    `pdbText` is a multi-MODEL PDB whose MODEL numbers ARE mode numbers -- the
    same file the ranking view reads -- so a pose is selected by IDENTITY and
    never by position in the file (#53's mistake, one layer down).
    """
    return """
const PV_CYS = %d, PV_POCKET = %s;
let PV = null, PV_SURF = null;
function pvLib(){ return window.$3Dmol || window.$3dmol || null; }
// CARBONS CARRY THE MODE COLOUR; every other element keeps its conventional one.
// Colouring a whole ligand by mode hides the chemistry a chemist reads it by --
// the sulfur, the halogen, the oxygens must look the same in every view.
function pvCarbon(col){
  const M = pvLib();
  return {prop:'elem',
          map: Object.assign({}, (M.elementColors||{}).defaultColors||{}, {C: col})};
}
function mountPose(elemId, pdbText, mode, recTxt){
  const M = pvLib(); if (!M) return;
  const host = document.getElementById(elemId); if (!host) return;
  if (!PV) PV = M.createViewer(host, {backgroundColor:'#eef1f6'});
  PV.clear(); PV_SURF = null;
  PV.addModel(recTxt, 'pdb');
  PV.setStyle({}, {cartoon:{color:'#c3ccd8', opacity:0.5}});
  // Cys113 in element colours, with a sphere on SG: at stick radius a lone
  // sulfur is easy to lose against the cartoon, and it is the atom the whole
  // screen is aimed at.
  PV.setStyle({resi:[PV_CYS]}, {stick:{radius:0.28, colorscheme:'default'},
                                cartoon:{color:'#c3ccd8', opacity:0.5}});
  PV.addStyle({resi:[PV_CYS], atom:'SG'}, {sphere:{radius:0.62}});
  const modes = [];
  (pdbText.match(/^MODEL\\s+(-?\\d+)/gm) || []).forEach(function(l){
    modes.push(parseInt(l.replace(/^MODEL\\s+/, ''), 10)); });
  PV.addModelsAsFrames(pdbText, 'pdb');
  let shown = modes.indexOf(mode);
  modes.forEach(function(mo, i){
    PV.setStyle({model: i+1}, (mo === mode)
      ? {stick:{radius:0.22, colorscheme: pvCarbon(0x0072ce)}} : {});
  });
  if (document.getElementById('pv-surf') &&
      document.getElementById('pv-surf').checked){
    // The pocket shell only, and never over Cys113 or the ligand: a mesh on top
    // of those hides the two things the panel exists to show.
    PV_SURF = PV.addSurface(M.SurfaceType.VDW, {opacity:0.62, color:'#b9c7db'},
      {and:[{model:0}, {resi:PV_POCKET}, {not:{resi:[PV_CYS]}}]});
  }
  PV.zoomTo(shown >= 0 ? {model: shown+1} : {resn:'MOL'});
  PV.zoom(0.5);
  // Size is read AFTER layout, twice: a canvas built before the grid settles
  // comes up off-centre and unrotatable, which has cost two rounds of "I can't
  // move it".
  requestAnimationFrame(function(){
    requestAnimationFrame(function(){ PV.resize(); PV.render(); }); });
  PV.render();                      // 3Dmol draws NOTHING without this
  return modes;
}
""" % (cys_resi, pocket)
