"""
Purpose: Show a docked pose IN THE POCKET, and animate an MD trajectory.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: docked SDF per candidate, the prepared receptor, GROMACS trajectories
Output: py3Dmol HTML for embedding in the Streamlit dossier

WHY THE RECEPTOR IS NOT OPTIONAL (issue #1, T_4 note: "poses were given without
the pocket for some reason??"). A ligand rendered alone is a conformer, not a
pose. Every claim a docked pose makes -- that a warhead reaches Cys113, that a
substituent occupies a subpocket, that a molecule is too large -- is a claim
about the ligand RELATIVE TO the protein, and none of it is checkable without
the protein on screen. So the receptor is drawn first and the ligand into it;
there is no ligand-only path in this module.

CYS113 IS ALWAYS HIGHLIGHTED. It is the residue the covalent approaches target
and the one a reader is looking for. Leaving them to find it in a 163-residue
ribbon is a needless obstacle.

THE MD ANIMATION IS THE PBC-CORRECTED TRAJECTORY, NOT prod.xtc. `whole.xtc` is
what `trjconv -pbc mol -center` produced. Animating the raw trajectory shows the
ligand teleporting across the box every time it crosses a periodic boundary --
which is what made the uncorrected RMSD 24x too large (D0038). A viewer reading
the raw file would see that artefact and reasonably conclude the ligand
dissociated.
"""

from __future__ import annotations

#: Mtime of THIS file at the moment it was imported. Frozen at import, so
#: comparing it with the file's current mtime is the only reliable way to tell
#: that a running process is executing stale code -- Streamlit re-runs the
#: script on every interaction but never re-imports helper modules.
LOADED_MTIME = __import__("os").stat(__file__).st_mtime

import subprocess
from pathlib import Path

RECEPTOR = Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdb")
DATA = Path("/data/lab_vm/append_only/inhibition")
GMX = Path("/data/lab_vm/envs/dwi_gromacs_cuda/bin/gmx")

DOCKING_DIRS = {
    "t1": DATA / "01_t1_de_novo" / "docking",
    "t2": DATA / "02_t2_atra_crem" / "docking",
    "t3": DATA / "03_t3_reinvent" / "docking",
    "t4": DATA / "04_t4_combinatorial" / "docking",
}
GROMACS_DIRS = {
    "t1": DATA / "01_t1_de_novo" / "gromacs",
    "t2": DATA / "02_t2_atra_crem" / "gromacs",
}

CATALYTIC_RESI = 113


class MovieTooLarge(RuntimeError):
    """The rendered movie would exceed what a browser will embed."""


def find_pose(approach: str, candidate_id: str,
              dock_id: str | None = None) -> Path | None:
    """The docked pose file for one candidate, or None.

    THE FILENAME CANNOT BE DERIVED FROM THE CANDIDATE ID, AND GUESSING IT FAILED
    SILENTLY. The covalent approaches name pose files by a SEPARATE `dock_id`
    (`d_31fe9e96d8bb_docked.sdf`) that shares no hash with `candidate_id`
    (`t3_31d6edc305b0`) — the overlap between the two sets is exactly zero
    across 4080 files. An earlier version built the name from `candidate_id`,
    found nothing, and rendered "no docked pose file found", which reads as
    missing data rather than as a lookup bug. Pass `dock_id` from the frame.

    THE TWO STRATA ALSO USE DIFFERENT LAYOUTS AND FORMATS:
      T_3, T_4 (gnina, covalent)  docking/d_<dock_id>_docked.sdf
      T_1, T_2 (Vina, non-covalent)  docking/poses/<candidate_id>_out.pdbqt
    Both are returned; `pose_html` dispatches on the suffix.
    """
    d = DOCKING_DIRS.get(approach)
    if d is None or not d.is_dir():
        return None

    # Non-covalent: Vina writes PDBQT into a poses/ subdirectory.
    poses = d / "poses"
    if poses.is_dir():
        for name in (f"{candidate_id}_out.pdbqt", f"{candidate_id}.pdbqt"):
            p = poses / name
            if p.is_file() and p.stat().st_size > 0:
                return p

    # Covalent: gnina writes SDF named by dock_id.
    if dock_id:
        for name in (f"{dock_id}_docked.sdf", f"{dock_id}.sdf"):
            p = d / name
            if p.is_file() and p.stat().st_size > 0:
                return p
    return None


def find_trajectory(approach: str, candidate_id: str,
                    replicate: int = 1) -> tuple[Path, Path] | None:
    """(PBC-corrected trajectory, its tpr) for one replicate, or None."""
    root = GROMACS_DIRS.get(approach)
    if root is None:
        return None
    wd = root / candidate_id
    for rep_dir in (wd / f"rep{replicate}", wd):
        whole, tpr = rep_dir / "whole.xtc", rep_dir / "prod.tpr"
        if whole.is_file() and tpr.is_file():
            return whole, tpr
    return None


#: How to drive a 3Dmol.js canvas. Written out because the viewer offers no
#: on-screen affordances whatsoever -- no buttons, and nothing indicates that
#: right-drag zooms or Ctrl+drag pans. A reader who tries only left-drag
#: concludes the view is stuck off-centre and unusable.
CONTROLS = """
| action | mouse | trackpad |
|---|---|---|
| **rotate** | left-click + drag | one-finger drag |
| **zoom** | scroll wheel, or right-click + drag | two-finger scroll / pinch |
| **pan** — shift the molecule sideways | middle-click + drag, or **Ctrl** + left-drag | **Ctrl** + one-finger drag |
| **slab** — cut away the front | **Ctrl** + right-drag | — |
| **reset** | change the *framing* control and back | same |
"""


def pose_html(sdf: Path, *, width: int | None = None, height: int = 620,
              surface: bool = True, zoom_on: str = "ligand") -> str:
    """The docked pose inside the receptor, with Cys113 called out.

    `width=None` makes the canvas fill its container. This matters more than it
    sounds: the viewer was fixed at 700 px inside a one-of-four column roughly
    300 px wide, so most of the scene sat off-canvas and no amount of mouse work
    brought it back. The complaint was that the image could not be centred; the
    cause was that the canvas was wider than the space it was drawn into.

    `zoom_on` sets the opening framing -- "ligand" fills the view with the
    binding site, "pocket" backs off to show the surrounding fold, "all" shows
    the whole protein for orientation.
    """
    import py3Dmol

    v = py3Dmol.view(width=width or "100%", height=height)
    v.addModel(RECEPTOR.read_text(), "pdb")
    v.setStyle({"model": 0}, {"cartoon": {"color": "spectrum", "opacity": 0.65}})
    # The catalytic cysteine, drawn as sticks so its SG is locatable by eye.
    v.addStyle({"model": 0, "resi": CATALYTIC_RESI},
               {"stick": {"colorscheme": "yellowCarbon", "radius": 0.22}})
    v.addResLabels({"model": 0, "resi": CATALYTIC_RESI},
                   {"fontSize": 11, "backgroundOpacity": 0.55})
    if surface:
        # Pocket surface only — a whole-protein surface hides the ligand.
        v.addSurface("VDW", {"opacity": 0.55, "color": "lightgrey"},
                     {"model": 0, "resi": list(range(CATALYTIC_RESI - 12,
                                                     CATALYTIC_RESI + 12))})
    # PDBQT is PDB plus partial-charge/type columns; py3Dmol parses it as PDB.
    # Vina writes every mode into one file, so only the FIRST is shown --
    # showing all nine overlaid looks like one impossible molecule.
    if sdf.suffix.lower() == ".pdbqt":
        text = sdf.read_text()
        first = text.split("ENDMDL")[0] + "ENDMDL\n" if "ENDMDL" in text else text
        v.addModel(first, "pdb")
    else:
        v.addModel(sdf.read_text(), "sdf")
    v.setStyle({"model": 1}, {"stick": {"colorscheme": "cyanCarbon",
                                        "radius": 0.16}})
    if zoom_on == "all":
        v.zoomTo()
    elif zoom_on == "pocket":
        v.zoomTo({"model": 0, "resi": list(range(CATALYTIC_RESI - 12,
                                                 CATALYTIC_RESI + 12))})
    else:
        v.zoomTo({"model": 1})
        # Back off slightly: zoomTo on a small ligand crops the pocket walls
        # out of frame, which is the context the pose exists to show.
        v.zoom(0.7)
    return v._make_html()


#: Streamlit embeds the viewer HTML inline. Past roughly this size the browser
#: renders nothing at all and the click looks like a dead button -- which is
#: exactly what a 12 MB movie did. Guarded rather than hoped about.
MAX_EMBED_BYTES = 4_000_000


def trajectory_html(xtc: Path, tpr: Path, *, n_frames: int = 15,
                    width: int = 700, height: int = 480) -> str | None:
    """An animated MD trajectory, or None if it cannot be built or is too big.

    WHY SO FEW FRAMES. `whole.xtc` already contains ONLY the solute -- the
    PBC-correction step wrote protein + ligand and dropped ~27,600 waters and
    ions. What remains is 2,391 atoms, and at 60 frames that is 12 MB of inline
    HTML. Streamlit embeds the viewer directly in the page and the browser
    silently renders nothing above a few MB, so the button appeared dead. Frame
    count is the only lever that does not require re-deriving an index.

    A BACKBONE SELECTION WAS TRIED AND IS WRONG HERE. `gmx select` indexes
    against `prod.tpr` (30,037 atoms) while `whole.xtc` holds 2,391, so every
    index past the solute overruns the trajectory and trjconv aborts with a
    mismatch. Any future selection must be built against the solute subset, not
    the full topology.
    """
    ndx = xtc.parent / "analysis.ndx"
    if not ndx.is_file():
        return None
    out = xtc.parent / f"movie_{n_frames}.pdb"
    if not out.is_file():
        try:
            subprocess.run(
                [str(GMX), "trjconv", "-s", str(tpr), "-f", str(xtc),
                 "-n", str(ndx), "-o", str(out), "-skip",
                 str(max(1, 1000 // n_frames))],
                input="0\n", capture_output=True, text=True, timeout=900,
                cwd=xtc.parent,
                env={"GMX_MAXBACKUP": "-1", "PATH": "/usr/bin:/bin"}, check=True)
        except Exception:  # noqa: BLE001 - a missing movie is not a page failure
            return None
    if not out.is_file() or out.stat().st_size == 0:
        return None
    if out.stat().st_size > MAX_EMBED_BYTES:
        log_size = out.stat().st_size / 1e6
        raise MovieTooLarge(
            f"{out.name} is {log_size:.1f} MB; above ~{MAX_EMBED_BYTES/1e6:.0f} MB "
            "the browser renders nothing and the control looks broken. Reduce "
            "n_frames.")

    import py3Dmol
    v = py3Dmol.view(width=width, height=height)
    v.addModelsAsFrames(out.read_text(), "pdb")
    v.setStyle({"cartoon": {"color": "spectrum", "opacity": 0.6}})
    v.addStyle({"hetflag": True},
               {"stick": {"colorscheme": "cyanCarbon", "radius": 0.18}})
    v.addStyle({"resi": CATALYTIC_RESI},
               {"stick": {"colorscheme": "yellowCarbon", "radius": 0.2}})
    v.zoomTo({"hetflag": True})
    v.animate({"loop": "forward", "interval": 80})
    return v._make_html()


def rmsd_series(approach: str, candidate_id: str) -> dict[int, list[tuple]]:
    """Per-replicate (time_ns, ligand RMSD nm) from each rep's rmsd.xvg.

    Returned per replicate rather than averaged. D0038's whole lesson was that
    one trajectory cannot separate "this ligand leaves" from "this trajectory
    wandered", so the spread between replicates IS the result and collapsing it
    to a mean would hide exactly what the replicates were run to show.
    """
    root = GROMACS_DIRS.get(approach)
    if root is None:
        return {}
    wd = root / candidate_id
    if not wd.is_dir():
        return {}
    out: dict[int, list[tuple]] = {}
    for rep_dir in sorted(wd.glob("rep*")):
        xvg = rep_dir / "rmsd.xvg"
        if not xvg.is_file():
            continue
        rows = []
        for line in xvg.read_text(errors="ignore").splitlines():
            line = line.strip()
            if not line or line[0] in "@#":
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
        if rows:
            try:
                out[int(rep_dir.name[3:])] = rows
            except ValueError:
                continue
    return out
