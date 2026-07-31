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


def find_pose(approach: str, candidate_id: str) -> Path | None:
    """The docked SDF for one candidate, or None.

    Candidate ids carry an approach prefix (`t4_ab12…`) while the docking
    filenames do not (`d_ab12…_docked.sdf`), so both spellings are tried rather
    than assuming one — a mismatch here shows up as "no pose available", which
    looks like missing data instead of a naming bug.
    """
    d = DOCKING_DIRS.get(approach)
    if d is None or not d.is_dir():
        return None
    stem = candidate_id.split("_", 1)[-1]
    for name in (f"{candidate_id}_docked.sdf", f"d_{stem}_docked.sdf",
                 f"{stem}_docked.sdf"):
        p = d / name
        if p.is_file() and p.stat().st_size > 0:
            return p
    hits = sorted(d.glob(f"*{stem}*_docked.sdf"))
    return hits[0] if hits else None


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


def pose_html(sdf: Path, *, width: int = 700, height: int = 480,
              surface: bool = True) -> str:
    """The docked pose inside the receptor, with Cys113 called out."""
    import py3Dmol

    v = py3Dmol.view(width=width, height=height)
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
    v.addModel(sdf.read_text(), "sdf")
    v.setStyle({"model": 1}, {"stick": {"colorscheme": "cyanCarbon",
                                        "radius": 0.16}})
    v.zoomTo({"model": 1})
    return v._make_html()


def trajectory_html(xtc: Path, tpr: Path, *, n_frames: int = 60,
                    width: int = 700, height: int = 480) -> str | None:
    """An animated MD trajectory of the solute, or None if it cannot be built.

    Converts the PBC-corrected trajectory to a multi-model PDB of protein +
    ligand only. Water and ions are dropped: they are ~90% of a 30k-atom system
    and would make the viewer unusable while showing nothing a reader needs.
    """
    out = xtc.parent / f"movie_{n_frames}.pdb"
    if not out.is_file():
        ndx = xtc.parent / "analysis.ndx"
        if not ndx.is_file():
            return None
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
