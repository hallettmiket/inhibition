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


def viewer_html(pdb_text: str, dist: list, labels: list, positions: list,
                three_js: str, nac_lo: float = 2.8, nac_hi: float = 4.2,
                elem_id: str = "gl") -> str:
    """A self-contained 3Dmol block: surface, charge colouring, labels, slider."""
    return f"""
<div class="glwrap">
  <div class="glbox"><div id="{elem_id}"></div></div>
  <div class="glctl">
    <button id="play">&#9654; play</button>
    <input id="frame" type="range" min="0" max="{max(0, len(dist) - 1)}" value="0">
    <span id="ftxt" class="mono"></span>
    <label><input id="surf" type="checkbox" checked> surface</label>
    <label><input id="labs" type="checkbox" checked> labels</label>
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
  const LO = {nac_lo}, HI = {nac_hi};
  const box = document.getElementById('{elem_id}');
  const raw = document.getElementById('pdbdata-{elem_id}').textContent;
  let viewer = null, frame = 0, timer = null, surf = null;

  // charge -> colour. 0 is a mid grey that stands off a white page; negative
  // ramps to red and positive to blue, so the charged patches still read.
  const GREY = [0x6e, 0x76, 0x82], NEG = [0x9e, 0x24, 0x1a], POS = [0x14, 0x3d, 0x7a];
  function chargeColour(atom) {{
    const b = Math.max(-1, Math.min(1, atom.b || 0));
    const t = Math.abs(b), to = b < 0 ? NEG : POS;
    const c = GREY.map(function(g, i) {{ return Math.round(g + (to[i] - g) * t); }});
    return (c[0] << 16) | (c[1] << 8) | c[2];
  }}

  function draw() {{
    viewer.setFrame(frame).then(function() {{
      viewer.setStyle({{}}, {{cartoon: {{color: '#2f3742', opacity: 1.0}}}});
      viewer.setStyle({{resn: 'MOL'}},
                      {{stick: {{radius: 0.26, colorscheme: 'yellowCarbon'}}}});
      viewer.removeAllLabels();
      if (document.getElementById('labs').checked) {{
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
      document.getElementById('ftxt').textContent =
        frame + ' / ' + (DSG.length - 1) + '   warhead\\u2192SG ' +
        (d == null ? '\\u2014' : d.toFixed(2) + ' \\u00c5' +
          (d >= LO && d <= HI ? '  (in attack window)' : ''));
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
    surf = viewer.addSurface(M.SurfaceType.VDW,
      {{opacity: 0.94, colorfunc: chargeColour}}, {{hetflag: false}});
    viewer.zoomTo({{resn: 'MOL'}});
    viewer.zoom(0.55);
    draw();
    viewer.resize();
  }}

  document.getElementById('frame').addEventListener('input', function(e) {{
    frame = +e.target.value; draw();
  }});
  document.getElementById('play').addEventListener('click', function(e) {{
    if (timer) {{ clearInterval(timer); timer = null; e.target.innerHTML = '&#9654; play'; return; }}
    e.target.innerHTML = '&#10073;&#10073; pause';
    timer = setInterval(function() {{
      frame = (frame + 1) % DSG.length;
      document.getElementById('frame').value = frame;
      draw();
    }}, 70);
  }});
  document.getElementById('surf').addEventListener('change', function(e) {{
    if (!viewer) return;
    viewer.setSurfaceMaterialStyle(surf.surfid, {{opacity: e.target.checked ? 0.72 : 0}});
    viewer.render();
  }});
  document.getElementById('labs').addEventListener('change', draw);
  // 3Dmol absolutely-positions its canvas, so the container must already have a
  // real height when the viewer is created -- building on DOMContentLoaded gave
  // a 0x0 box and a blank viewer.
  window.addEventListener('load', function() {{
    requestAnimationFrame(function() {{ requestAnimationFrame(boot); }});
  }});
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
"""
