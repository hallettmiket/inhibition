"""
Purpose: turn a finished 100 ns trajectory into an embeddable 3Dmol movie.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: a completed md_residence rep directory (prod.tpr + prod.xtc + fit.ndx)
Output: a multi-model PDB, and the HTML/JS block that renders it

@tt8804 asked for surface, red/blue charge colouring, and labelled key residues.
`elevation_report` already builds exactly that, but its viewer is written inline
into one 500-line f-string tied to the elevation experiment's own variables, so
it cannot be called for a single MD-priority molecule. This carries the parts
that generalise.

WHAT IS REUSED RATHER THAN REWRITTEN. `elevation_report.surface_payload` does the
data preparation -- formal charge into the B-factor column, key-residue label
anchors, and the warhead-to-SG distance computed FROM THE RENDERED COORDINATES so
the number on screen cannot disagree with the structure on screen. It also
refuses outright if the residue numbering does not match what the labels claim,
which is the check that stops a viewer confidently labelling the wrong residue.

TRAJECTORY PREPARATION IS NOT COSMETIC. Three passes, in this order:
  1. `-pbc whole` repairs molecules split across the periodic box. Skipping it
     draws bonds stretching across the whole cell.
  2. `-center` on the complex, `-pbc mol`, so the protein does not wander out of
     view.
  3. `-fit rot+trans` on CA atoms, so what you see is the LIGAND moving relative
     to the protein rather than the whole system tumbling. Without it every
     trajectory looks unstable.
"""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path

GMX = Path("/data/lab_vm/envs/dwi_gromacs_cuda/bin.AVX2_256/gmx")
log = logging.getLogger("md-movie")

#: Frames in the movie. 10,001 frames of a 100 ns run is far more than a viewer
#: can show and would make the page enormous; ~150 is smooth at a readable rate.
N_FRAMES = 120


def _count_frames(rep: Path, tpr: Path, xtc: Path) -> int:
    """How many frames the trajectory actually stored, from gmx check."""
    r = subprocess.run([str(GMX), "check", "-f", str(xtc)], cwd=rep,
                       capture_output=True, text=True, timeout=900)
    for ln in (r.stderr + r.stdout).splitlines():
        if ln.strip().startswith("Step"):
            parts = ln.split()
            if len(parts) >= 2 and parts[1].isdigit():
                return int(parts[1])
    return 0


def build_movie_pdb(rep: Path, dest: Path, total_ps: float = 100_000.0,
                    n_frames: int = N_FRAMES) -> Path | None:
    """PBC-corrected, CA-fitted multi-model PDB of protein + ligand."""
    tpr, xtc, ndx = rep / "prod.tpr", rep / "prod.xtc", rep / "fit.ndx"
    if not (tpr.is_file() and xtc.is_file()):
        log.warning("no prod.tpr/prod.xtc in %s", rep)
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".whole.xtc")

    def run(args, stdin_groups):
        return subprocess.run([str(GMX)] + args, cwd=rep, input=stdin_groups,
                              capture_output=True, text=True, timeout=3600)

    # 1. repair molecules broken across the box, subsampling with -skip.
    #
    # NOT `-dt`: the trajectory is saved every 10 ps, and `-dt` keeps only frames
    # whose TIME is a multiple of the value given. A computed dt of 666 ps
    # therefore matched almost nothing and produced 31 frames instead of 150.
    # `-skip` takes every Nth stored frame and cannot fall out of step with the
    # save interval.
    n_stored = _count_frames(rep, tpr, xtc)
    skip = max(1, round(n_stored / n_frames)) if n_stored else 67
    r = run(["trjconv", "-s", "prod.tpr", "-f", "prod.xtc", "-o", str(tmp),
             "-pbc", "whole", "-skip", str(skip)], "System\n")
    if r.returncode != 0:
        log.warning("trjconv (whole) failed: %s", r.stderr[-300:])
        return None

    # An index carrying BOTH the CA fit group and protein+ligand.
    #
    # The run's own fit.ndx holds only `name_CA` and the ligand, and passing it
    # to trjconv REPLACES the default groups -- so "System" stops existing and
    # the selection fails. Writing a combined index is the fix; selecting
    # "System" for output would also drag 24,774 water atoms into the viewer.
    # HEAVY ATOMS ONLY, SELECTED BY NAME.
    #
    # With hydrogens the page reached 21.5 MB, over the 16 MB publish limit --
    # and hydrogens are invisible under a cartoon and a VdW surface, so they were
    # pure payload.
    #
    # `gmx select` rather than `make_ndx`: the make_ndx version referred to the
    # new group by NUMBER (`2 | 17`), and the number depends on how many default
    # groups the system happens to generate. It worked for one molecule and threw
    # a fatal index error on the next, whose ion content differs. A selection
    # expression names what it wants and cannot drift.
    mndx = rep / "movie.ndx"
    if not mndx.is_file():
        r = run(["select", "-s", "prod.tpr", "-on", str(mndx), "-select",
                 '"C-alpha" name CA; '
                 '"complex_heavy" (group "Protein" or resname MOL) '
                 'and not name "H*"'], "")
        if r.returncode != 0 or not mndx.is_file():
            log.warning("gmx select failed: %s", r.stderr[-300:])
            tmp.unlink(missing_ok=True)
            return None

    # 2. STOP THE LIGAND JUMPING PERIODIC IMAGES.
    #
    # This pass was missing and it mattered. `-pbc whole` only repairs molecules
    # broken across the box; it does not stop a ligand hopping to a neighbouring
    # image between frames. Without it the warhead-to-SG distance computed from
    # these frames reached 97 A on a trajectory whose RMSD says the ligand never
    # left -- a plot that would have contradicted the residence number on the
    # same page. `-pbc nojump` cannot be combined with `-fit`, hence a pass of
    # its own.
    tmp2 = dest.with_suffix(".nojump.xtc")
    r = run(["trjconv", "-s", "prod.tpr", "-f", str(tmp), "-o", str(tmp2),
             "-pbc", "nojump"], "System\n")
    tmp.unlink(missing_ok=True)
    if r.returncode != 0:
        log.warning("trjconv (nojump) failed: %s", r.stderr[-300:])
        return None

    # 3. centre on the protein and fit on CA, so the LIGAND is what moves.
    #    With both flags trjconv asks for fit group, then centre group, then
    #    output group, in that order.
    r = run(["trjconv", "-s", "prod.tpr", "-f", str(tmp2), "-n", str(mndx),
             "-o", str(dest), "-fit", "rot+trans", "-center"],
            "C-alpha\nC-alpha\ncomplex_heavy\n")
    tmp2.unlink(missing_ok=True)
    if r.returncode != 0 or not dest.is_file():
        log.warning("trjconv (fit) failed: %s", r.stderr[-300:])
        return None
    n = dest.read_text().count("MODEL")
    log.info("movie: %d frames -> %s", n, dest.name)
    return dest if n else None


#: The MD system renumbers from 1, so the crystal's Cys113 is residue 63 in the
#: rendered movie. `elevation_report.PIN1_OFFSET` is the same constant; it is
#: repeated here rather than imported because this module is deliberately free of
#: script-level imports. Selecting residue 113 instead picks up a GLUTAMATE -- a
#: long side chain with two carboxyl oxygens that looks nothing like a cysteine,
#: which is exactly how the wrong residue got drawn and shipped.
PIN1_OFFSET = 50
CYS113_RESI = 113 - PIN1_OFFSET

#: Radius of the surface shell around Cys113, Angstrom. The mesh is rebuilt every
#: frame, so it has to be small enough to be built ~10x a second; this covers the
#: pocket the movie exists to show and leaves the rest to the cartoon.
SURF_SHELL_A = 14


def viewer_html(pdb_text: str, dist: list, labels: list, positions: list,
                three_js: str, nac_lo: float | None = None,
                nac_hi: float | None = None,
                total_ps: float | None = None, fate: str | None = None,
                elem_id: str = "gl") -> str:
    """A self-contained 3Dmol block: surface, charge colouring, labels, slider.

    EVERY CONTROL ID IS NAMESPACED BY `elem_id`. They used to be bare -- `frame`,
    `play`, `surf` -- which is fine for one viewer per page and breaks the moment
    two share a document: `getElementById` returns the first match, so all the
    sliders drive the first movie and the rest are inert. The shortlist report
    puts four viewers on one page, so the ids have to be unique per viewer.
    """
    # THE BAND COMES FROM THE GATE, not from defaults in this signature.
    # `nac_lo=2.8, nac_hi=4.2` were keyword defaults no caller ever overrode --
    # a pin that could not announce itself (catalogue #32/#35) -- so the viewer
    # said "(in attack window)" for a pose at 4.0 A that the sweep page scored
    # as not engaged. Passing an explicit value still works; passing nothing now
    # means "ask the criterion" rather than "assume the screen's window".
    from shared import nac_criterion as _nac
    _dlo, _dhi = _nac.attack_ready_window()
    _lo = _dlo if nac_lo is None else float(nac_lo)
    _hi = _dhi if nac_hi is None else float(nac_hi)
    return f"""
<div class="glwrap">
  <div class="glbox"><div id="{elem_id}"></div></div>
  <div class="glctl">
    <button id="{elem_id}-play">&#9654; play</button>
    <input id="{elem_id}-frame" type="range" min="0" max="{max(0, len(dist) - 1)}" value="0">
    <span id="{elem_id}-ftxt" class="mono"></span>
    <span id="{elem_id}-sstat" class="mono"></span>
    <label><input id="{elem_id}-surf" type="checkbox" checked> surface</label>
    <label><input id="{elem_id}-labs" type="checkbox" checked> labels</label>
  </div>
</div>
<script>{three_js}</script>
<script type="text/plain" id="pdbdata-{elem_id}">{pdb_text}</script>
<script>
(function(){{
  // The UMD footer only guarantees the global as 3Dmol; $3Dmol is the
  // conventional alias and is not exported by every build. Bind whichever exists.
  const M = window.$3Dmol || window['3Dmol'];
  const DSG = {json.dumps(dist)}, LABELS = {json.dumps(labels)},
        LPOS = {json.dumps(positions)};
  const LO = {_lo}, HI = {_hi};
  // SIMULATION TIME PER FRAME. A frame index says nothing about when in the run
  // it happened, and runs are no longer the same length -- frame 60 of a 1.2 ns
  // sweep and frame 60 of a 10 ns one are 8 ns apart.
  const TOTAL_PS = {("null" if total_ps is None else float(total_ps))};
  const FATE = {json.dumps(fate) if fate else "null"};
  const box = document.getElementById('{elem_id}');
  const raw = document.getElementById('pdbdata-{elem_id}').textContent;
  let viewer = null, frame = 0, timer = null, surf = null, surfFrame = -1;

  // charge -> colour. 0 is a mid grey that stands off a white page; negative
  // ramps to red and positive to blue, so the charged patches still read.
  const GREY = [0x6e, 0x76, 0x82], NEG = [0x9e, 0x24, 0x1a], POS = [0x14, 0x3d, 0x7a];
  function chargeColour(atom) {{
    const b = Math.max(-1, Math.min(1, atom.b || 0));
    const t = Math.abs(b), to = b < 0 ? NEG : POS;
    const c = GREY.map(function(g, i) {{ return Math.round(g + (to[i] - g) * t); }});
    return (c[0] << 16) | (c[1] << 8) | c[2];
  }}

  function surfaceOn() {{
    const c = document.getElementById('{elem_id}-surf');
    return !c || c.checked;
  }}
  // THE SURFACE IS REBUILT ON EVERY FRAME (@tt8804), so it moves with the
  // backbone under it and the ligand can never clip through a wall belonging to
  // a different frame.
  //
  // What makes that affordable is the SHELL: the mesh covers only residues
  // within {SURF_SHELL_A} A of Cys113, not the whole protein. A full-protein VDW
  // mesh is the expensive call in this viewer and rebuilding it 126 times is not
  // watchable; the pocket shell is a small fraction of the atoms and is the only
  // part anyone looks at. Everything outside it is carried by the cartoon, which
  // tracks every frame for free.
  //
  // removeSurface FIRST -- addSurface stacks meshes, and stacking one per frame
  // is how the viewer dies.
  function buildSurface() {{
    if (!viewer) return;
    if (surf) {{ try {{ viewer.removeSurface(surf.surfid); }} catch (e) {{}} surf = null; }}
    if (!surfaceOn()) {{ surfFrame = frame; return; }}
    surf = viewer.addSurface(M.SurfaceType.VDW,
      {{opacity: 0.98, colorfunc: chargeColour}},
      {{not: {{or: [{{resn: 'MOL'}}, {{resi: {CYS113_RESI}}}]}}}});
    surfFrame = frame;
  }}
  // Nothing to warn about any more: the mesh is rebuilt with the frame, so it
  // can never describe a different one. The readout just says what it covers.
  // THE SURFACE IS NEVER HIDDEN (@tt8804). Rebuilding it per frame could not keep
  // up with playback and it flickered out; hiding it while it lagged left the
  // frame with no surface at all. It now stays on screen the whole time and is
  // REFRESHED WHEN PLAYBACK STOPS -- on pause, and on releasing the slider. While
  // the video runs the mesh is one frame's shell over a moving cartoon, which is
  // the compromise being chosen deliberately; the readout says so.
  function markStale() {{
    const el = document.getElementById('{elem_id}-sstat');
    if (!el) return;
    const stale = surfFrame !== frame;
    el.textContent = (surfaceOn() && stale) ? 'surface: frame ' + surfFrame
                                            + ' — refreshes on pause' : '';
    el.className = (surfaceOn() && stale) ? 'mono stale' : 'mono';
  }}

  function draw() {{
    viewer.setFrame(frame).then(function() {{
      markStale();
      viewer.setStyle({{}}, {{cartoon: {{color: '#2f3742', opacity: 1.0}}}});
      viewer.setStyle({{resn: 'MOL'}},
                      {{stick: {{radius: 0.26, colorscheme: 'yellowCarbon'}}}});
      // CYS113 IN STICKS. It is the atom every distance in this project is
      // measured to, and under a 90%-opaque surface it was invisible. Drawn in
      // green carbons so it reads as protein rather than as a second ligand,
      // and kept in front of the surface so the approach can be seen.
      viewer.setStyle({{resi: {CYS113_RESI}}},
                      {{stick: {{radius: 0.3, colorscheme: 'greenCarbon'}},
                        cartoon: {{color: '#2f3742', opacity: 1.0}}}});
      // the sulfur itself -- every distance in this project is measured to it
      viewer.addStyle({{resi: {CYS113_RESI}, atom: 'SG'}},
                      {{sphere: {{radius: 0.75, color: '#f0c000'}}}});
      viewer.removeAllLabels();
      if (document.getElementById('{elem_id}-labs').checked) {{
        LABELS.forEach(function(L, i) {{
          const p = LPOS[frame] && LPOS[frame][i];
          if (!p) return;
          viewer.addLabel(L.text, {{position: {{x: p[0], y: p[1], z: p[2]}},
            backgroundColor: L.kind === 'reactive' ? '#b3261e' : '#003087',
            backgroundOpacity: 0.82, fontColor: 'white', fontSize: 11,
            borderThickness: 0}});
        }});
      }}
      const d = DSG[frame];
      document.getElementById('{elem_id}-ftxt').textContent =
        (TOTAL_PS == null || DSG.length < 2
           ? frame + ' / ' + (DSG.length - 1)
           : (TOTAL_PS * frame / (DSG.length - 1) / 1000).toFixed(2) + ' ns'
             + ' / ' + (TOTAL_PS / 1000).toFixed(TOTAL_PS >= 10000 ? 0 : 1) + ' ns')
        + '   warhead\\u2192SG ' +
        (d == null ? '\\u2014' : d.toFixed(2) + ' \\u00c5' +
          (d >= LO && d < HI ? '  (attack ready)' : '')) +
        // WHY THE MOVIE ENDS. On the last frame, say whether the run stopped
        // because the molecule left or because it reached the cap -- otherwise
        // the end of every movie looks identical.
        (FATE && frame === DSG.length - 1 ? '   \\u2014 ' + FATE : '');
      viewer.render();
    }});
  }}

  function boot() {{
    if (!M || !box.clientHeight) {{ return setTimeout(boot, 120); }}
    viewer = M.createViewer(box, {{backgroundAlpha: 0}});
    viewer.addModelsAsFrames(raw, 'pdb');
    // Charge lives in the B-factor column, so a red-white-blue gradient over b
    // colours the surface by formal charge rather than by anything cosmetic.
    // NEUTRAL IS GREY, NOT WHITE.
    //
    // The stock RWB gradient maps charge 0 to WHITE, and almost every protein
    // atom is neutral -- so on a white page the surface disappeared entirely and
    // the structure read as a blank silhouette. A red-GREY-blue ramp keeps the
    // charge information and leaves the uncharged bulk visible.
    // TWO SURFACES, AND THE LIGAND IS NOT PART OF THE PROTEIN ONE.
    // `hetflag:false` alone did not exclude it -- the ligand is written as ATOM
    // records in the fitted movie -- so the charge surface closed over the
    // warhead and buried the thing the viewer exists to show. The protein is
    // explicitly everything EXCEPT resn MOL, at half opacity so the ligand reads
    // through it; the ligand carries its own surface in the same yellow as its
    // sticks, so it is legible as one object rather than two.
    // CYS113 IS CUT OUT OF THE SURFACE, not just drawn under it. At 98% opacity
    // a stick inside the protein is invisible, so adding the residue in sticks
    // did nothing on its own -- the surface simply covered it. Excluding both the
    // ligand and residue 113 leaves an opening at the reaction site, which is
    // also the one place the reader needs to see into.
    // THE SURFACE IS A MESH, NOT A STYLE. `setFrame` regenerates the cartoon,
    // sticks and spheres from the new coordinates on every frame; a surface is
    // triangulated once, at the coordinates it was asked for, and a frame change
    // never recomputes it. Built once at load it stayed frozen at frame 0 while
    // the cartoon animated underneath -- and the Cys113 cut-out stayed carved
    // where the ligand USED to be, which is the one spot the reader is looking.
    //
    // Rebuilding every frame is correct and too slow to play through: the mesh is
    // the expensive call in this viewer. So it is rebuilt ON RELEASE -- slider
    // let go, playback paused -- and marked STALE in between, because a shell
    // describing a different frame with nothing saying so is precisely the
    // plausible-and-wrong rendering this project keeps producing.
    buildSurface();
    // NO SURFACE ON THE LIGAND (@tt8804: too distracting). It is drawn as
    // yellow sticks only, so the warhead and its approach vector stay readable
    // against the protein surface instead of competing with it.
    viewer.zoomTo({{resn: 'MOL'}});
    viewer.zoom(0.55);
    draw();
    viewer.resize();
  }}

  // `input` fires continuously through the drag, `change` fires once on release.
  // The cheap redraw rides the drag; the expensive mesh waits for the release.
  document.getElementById('{elem_id}-frame').addEventListener('input', function(e) {{
    frame = +e.target.value; draw();
  }});
  document.getElementById('{elem_id}-frame').addEventListener('change', function(e) {{
    frame = +e.target.value; draw();
  }});
  document.getElementById('{elem_id}-frame').addEventListener('change', function(e) {{
    frame = +e.target.value; draw();
    if (!timer) {{ buildSurface(); markStale(); viewer.render(); }}
  }});
  document.getElementById('{elem_id}-play').addEventListener('click', function(e) {{
    if (timer) {{ clearInterval(timer); timer = null; e.target.innerHTML = '&#9654; play';
                  buildSurface(); markStale(); viewer.render(); return; }}
    e.target.innerHTML = '&#10073;&#10073; pause';
    timer = setInterval(function() {{
      frame = (frame + 1) % DSG.length;
      document.getElementById('{elem_id}-frame').value = frame;
      draw();
    }}, 70);
  }});
  document.getElementById('{elem_id}-surf').addEventListener('change', function(e) {{
    if (!viewer) return;
    buildSurface(); markStale(); viewer.render();
  }});
  document.getElementById('{elem_id}-labs').addEventListener('change', draw);
  // 3Dmol absolutely-positions its canvas, so the container must already have a
  // real height when the viewer is created -- building on DOMContentLoaded gave
  // a 0x0 box and a blank viewer.
  // BUILD WHEN THE PANEL OPENS, not on window load. 3Dmol absolutely-positions
  // its canvas, so the container must already have a real height -- and a closed
  // <details> has none. Booting on load inside a collapsed panel gives a 0x0 box
  // and a viewer that renders nothing, which is exactly what happened when the
  // panels were changed to start closed: the movies were there, sized to nothing.
  var built = false;
  function boot1() {{
    if (built) return;
    built = true;
    requestAnimationFrame(function() {{ requestAnimationFrame(boot); }});
  }}
  (function() {{
    var host = document.getElementById('{elem_id}');
    var det = host && host.closest ? host.closest('details') : null;
    if (det) {{
      if (det.open) window.addEventListener('load', boot1);
      det.addEventListener('toggle', function() {{ if (det.open) boot1(); }});
    }} else {{
      window.addEventListener('load', boot1);
    }}
  }})();
}})();
</script>
"""


VIEWER_CSS = """
.glwrap { margin: 1.2rem 0; }
.glbox { position: relative; width: 100%; height: 500px;
         border: 1px solid #c3c9d4; border-radius: 4px; overflow: hidden;
         background: #eef1f6; }
.glbox > div { position: absolute; inset: 0; }
.glctl { display: flex; align-items: center; gap: .9rem; margin-top: .5rem;
         font-size: .85rem; flex-wrap: wrap; }
.glctl input[type=range] { flex: 1 1 220px; }
.glctl button { border: 1px solid #003087; background: #fff; color: #003087;
                border-radius: 3px; padding: .25rem .7rem; cursor: pointer; }
/* The surface is rebuilt on release, so while scrubbing or playing it describes
   an older frame. Say which one, rather than letting it pass as current. */
.glctl .stale { color: #8a5a00; background: #fdf0dc; border-radius: 99px;
                padding: .05rem .5rem; font-size: .78rem; }
"""
