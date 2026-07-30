"""
Purpose: Measure whether a docked pose survives explicit water, comparably to
         the implicit-solvent run.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-29
Input: a completed GROMACS workdir (prod.tpr, prod.xtc)
Output: ligand RMSD and contact statistics per candidate

THE QUESTION. Under GB implicit solvent two T_1 candidates left the pocket
outright (ligand RMSD 9.0 and 7.3 nm, engaged in 0.07 and 0.14 of frames).
Implicit solvent has no water molecules, so nothing could have held them there.
This tier asks whether that behaviour survives real water, and it is the only
comparison in the project that puts the same molecule under two solvent models.

WHAT IS COMPARABLE AND WHAT IS NOT -- read this before putting the two side by
side. Ligand RMSD after superposing on the protein IS the same quantity in both
runs: nanometres of ligand displacement from the starting pose, protein motion
removed. It can be compared directly.

The CONTACT COUNTS CANNOT. The implicit-solvent metric counted heavy-atom PAIRS
within 0.45 nm; `gmx mindist -on` counts contacts by its own definition and
returns numbers several-fold smaller on the same complex. Reporting them in one
table as though they were one metric would manufacture a difference that is
purely a definition change. The GROMACS contact count is therefore named
distinctly and carries its own tool attribution.

WHY GROMACS' OWN TOOLS. parmed cannot read XTC and MDAnalysis is not installed.
Rather than add a trajectory parser, `gmx rms` and `gmx mindist` do the work and
their XVG output is parsed. This keeps the dependency surface unchanged.
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

GMX_ENV = Path("/data/lab_vm/envs/dwi_gromacs_cuda")
LIGAND_RESNAMES = ("MOL", "LIG")
CONTACT_CUTOFF_NM = 0.45

# The implicit-solvent run's own threshold, reused so "engaged" means the same
# thing in both: a frame retaining at least a quarter of the starting contacts.
ENGAGED_FRACTION = 0.25


class AnalysisError(RuntimeError):
    """The trajectory could not be analysed or the result cannot be trusted."""


def _gmx(*args: str, cwd: Path, log_name: str, stdin: str | None = None) -> str:
    cmd = [str(GMX_ENV / "bin" / "gmx"), *args]
    proc = subprocess.run(cmd, cwd=cwd, input=stdin, capture_output=True,
                          text=True, timeout=7200,
                          env={"GMX_MAXBACKUP": "-1", "PATH": "/usr/bin:/bin"})
    (cwd / log_name).write_text((proc.stdout or "") + "\n" + (proc.stderr or ""),
                                encoding="utf-8")
    if proc.returncode != 0:
        raise AnalysisError(f"gmx {args[0]} failed ({proc.returncode}); "
                            f"see {cwd/log_name}")
    return proc.stdout or ""


def _read_xvg(path: Path) -> np.ndarray:
    """Two-column XVG (time, value) as an (n, 2) array, comments stripped."""
    if not path.is_file():
        raise AnalysisError(f"no XVG output at {path}")
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line or line[0] in "@#":
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    if not rows:
        raise AnalysisError(f"{path.name} contained no data rows")
    return np.asarray(rows)


def _ligand_resname(tpr: Path, wd: Path) -> str:
    """Whichever ligand residue name this system actually uses."""
    for name in LIGAND_RESNAMES:
        try:
            _gmx("select", "-s", str(tpr), "-select",
                 f'resname {name} and not name "H*"', "-on", "lig.ndx",
                 cwd=wd, log_name="select_lig.log")
            if (wd / "lig.ndx").is_file() and (wd / "lig.ndx").stat().st_size > 0:
                return name
        except AnalysisError:
            continue
    raise AnalysisError(
        f"no ligand residue found among {LIGAND_RESNAMES}; the analysis cannot "
        "guess which residue is the ligand")


def analyse(wd: Path) -> dict:
    """Ligand RMSD and contacts over the production trajectory."""
    tpr, xtc = wd / "prod.tpr", wd / "prod.xtc"
    for p in (tpr, xtc):
        if not p.is_file() or p.stat().st_size == 0:
            raise AnalysisError(f"{wd.name}: missing or empty {p.name}")

    resname = _ligand_resname(tpr, wd)
    _gmx("select", "-s", str(tpr), "-select", "name CA", "-on", "fit.ndx",
         cwd=wd, log_name="select_fit.log")
    both = wd / "both.ndx"
    both.write_text((wd / "fit.ndx").read_text() + (wd / "lig.ndx").read_text(),
                    encoding="utf-8")

    # Group 0 = protein CA (superposition), group 1 = ligand heavy atoms.
    _gmx("rms", "-s", str(tpr), "-f", str(xtc), "-n", str(both),
         "-o", "rmsd.xvg", "-tu", "ns", cwd=wd, log_name="rms.log",
         stdin="0\n1\n")
    rmsd = _read_xvg(wd / "rmsd.xvg")

    _gmx("mindist", "-s", str(tpr), "-f", str(xtc), "-n", str(both),
         "-on", "numcont.xvg", "-d", str(CONTACT_CUTOFF_NM), "-tu", "ns",
         cwd=wd, log_name="mindist.log", stdin="1\n0\n")
    cont = _read_xvg(wd / "numcont.xvg")

    if rmsd.shape[0] < 4:
        raise AnalysisError(f"{wd.name}: only {rmsd.shape[0]} frames; refusing "
                            "to summarise a trajectory this short")

    r = rmsd[:, 1]
    c = cont[:, 1]
    start = c[0] if c[0] > 0 else 1.0
    engaged = float((c >= ENGAGED_FRACTION * start).mean())

    return {
        "ligand_resname": resname,
        "n_frames_analysed": int(rmsd.shape[0]),
        "ns_analysed": round(float(rmsd[-1, 0]), 3),
        # Directly comparable with the implicit-solvent run's
        # ligand_rmsd_nm_* fields: same quantity, same units, protein removed.
        "explicit_ligand_rmsd_nm_mean": round(float(r.mean()), 4),
        "explicit_ligand_rmsd_nm_final": round(float(r[-1]), 4),
        "explicit_ligand_rmsd_nm_max": round(float(r.max()), 4),
        # NOT comparable with mean_contacts from the implicit run: different
        # definition, different tool. Named so it cannot be mistaken for it.
        "gmx_contacts_mean": round(float(c.mean()), 1),
        "gmx_contacts_start": int(c[0]),
        "explicit_frac_frames_engaged": round(engaged, 4),
        "contact_cutoff_nm": CONTACT_CUTOFF_NM,
        "contacts_metric": "gmx mindist -on; NOT the heavy-atom pair count "
                           "used by the implicit-solvent tier",
        "solvent": "explicit TIP3P",
        "not_a_ranking": "describes one complex under real water; feeds no gate",
    }
