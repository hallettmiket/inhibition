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


def _prepare_pbc(wd: Path, tpr: Path, xtc: Path, resname: str) -> tuple[Path, Path]:
    """Make molecules whole and centre on the protein BEFORE measuring anything.

    THIS STEP IS NOT OPTIONAL AND ITS ABSENCE IS NOT VISIBLE. A ligand that
    crosses a periodic boundary reappears on the far side of the box, and
    `gmx rms` then measures the distance to its own periodic image. On
    t2_bc8a4b62eb0e that turned a mean ligand RMSD of 0.203 nm into 4.856 nm --
    a molecule that barely moved, reported as having left the pocket.

    The tell was available and went unexamined: the uncorrected maximum was
    9.159 nm against a largest box dimension of 8.13 nm, i.e. a displacement
    larger than the box that contains it. Any "drift" close to half a box
    dimension should be assumed to be imaging until this step is applied.

    Three groups are written in a fixed order because trjconv and rms both take
    them positionally and getting them the wrong way round silently measures
    something else: 0 = protein+ligand (centred and written out), 1 = protein CA
    (superposition), 2 = ligand heavy atoms (the quantity).
    """
    _gmx("select", "-s", str(tpr), "-select",
         f'group "Protein" or resname {resname}', "-on", "solute.ndx",
         cwd=wd, log_name="select_solute.log")
    _gmx("select", "-s", str(tpr), "-select", "name CA", "-on", "fit.ndx",
         cwd=wd, log_name="select_fit.log")
    _gmx("select", "-s", str(tpr), "-select",
         f'resname {resname} and not name "H*"', "-on", "lig.ndx",
         cwd=wd, log_name="select_lig.log")
    ndx = wd / "analysis.ndx"
    ndx.write_text("".join((wd / f).read_text()
                           for f in ("solute.ndx", "fit.ndx", "lig.ndx")),
                   encoding="utf-8")

    whole = wd / "whole.xtc"
    # stdin "1\n0\n": centre on group 1 (CA), output group 0 (protein+ligand).
    _gmx("trjconv", "-s", str(tpr), "-f", str(xtc), "-n", str(ndx),
         "-o", str(whole), "-pbc", "mol", "-center",
         cwd=wd, log_name="trjconv_pbc.log", stdin="1\n0\n")
    if not whole.is_file() or whole.stat().st_size == 0:
        raise AnalysisError(f"{wd.name}: trjconv produced no corrected "
                            "trajectory; refusing to measure the raw one")
    return whole, ndx


def analyse(wd: Path) -> dict:
    """Ligand RMSD and contacts over the PBC-corrected production trajectory."""
    tpr, xtc = wd / "prod.tpr", wd / "prod.xtc"
    for p in (tpr, xtc):
        if not p.is_file() or p.stat().st_size == 0:
            raise AnalysisError(f"{wd.name}: missing or empty {p.name}")

    resname = _ligand_resname(tpr, wd)
    whole, ndx = _prepare_pbc(wd, tpr, xtc, resname)

    # Groups are positional: 1 = CA (fit), 2 = ligand heavy atoms (measured).
    _gmx("rms", "-s", str(tpr), "-f", str(whole), "-n", str(ndx),
         "-o", "rmsd.xvg", "-tu", "ns", cwd=wd, log_name="rms.log",
         stdin="1\n2\n")
    rmsd = _read_xvg(wd / "rmsd.xvg")

    _gmx("mindist", "-s", str(tpr), "-f", str(whole), "-n", str(ndx),
         "-on", "numcont.xvg", "-d", str(CONTACT_CUTOFF_NM), "-tu", "ns",
         cwd=wd, log_name="mindist.log", stdin="2\n1\n")
    cont = _read_xvg(wd / "numcont.xvg")

    if rmsd.shape[0] < 4:
        raise AnalysisError(f"{wd.name}: only {rmsd.shape[0]} frames; refusing "
                            "to summarise a trajectory this short")

    r = rmsd[:, 1]
    c = cont[:, 1]

    # A displacement approaching half a box dimension is the signature of
    # residual imaging, not motion. The correction above should make this
    # impossible; it is checked rather than trusted, because the failure is
    # silent and produced a 24-fold error once already.
    box_half_nm = 3.3
    if float(r.max()) > box_half_nm:
        log.warning("%s: max ligand RMSD %.2f nm exceeds %.1f nm even after PBC "
                    "correction. Treat as suspect and inspect the trajectory "
                    "before reporting it.", wd.name, float(r.max()), box_half_nm)

    start = c[0] if c[0] > 0 else 1.0
    engaged = float((c >= ENGAGED_FRACTION * start).mean())

    return {
        "ligand_resname": resname,
        "n_frames_analysed": int(rmsd.shape[0]),
        "ns_analysed": round(float(rmsd[-1, 0]), 3),
        "pbc_corrected": True,
        "explicit_ligand_rmsd_nm_mean": round(float(r.mean()), 4),
        "explicit_ligand_rmsd_nm_final": round(float(r[-1]), 4),
        "explicit_ligand_rmsd_nm_max": round(float(r.max()), 4),
        "explicit_rmsd_suspect": bool(float(r.max()) > box_half_nm),
        "gmx_contacts_mean": round(float(c.mean()), 1),
        "gmx_contacts_start": int(c[0]),
        "explicit_frac_frames_engaged": round(engaged, 4),
        "contact_cutoff_nm": CONTACT_CUTOFF_NM,
        "contacts_metric": "gmx mindist -on; NOT the heavy-atom pair count "
                           "used by the implicit-solvent tier",
        "solvent": "explicit TIP3P",
        "not_a_ranking": "describes one complex under real water; feeds no gate",
    }


PROT_RMSD_XVG = "rmsd_protein.xvg"


def protein_rmsd(wd: Path) -> Path | None:
    """CA RMSD over the SAME corrected trajectory the ligand RMSD is measured on.

    THE LIGAND NUMBER ALONE CANNOT DISTINGUISH TWO OPPOSITE EVENTS. Ligand RMSD
    is measured after superposing on protein CA, so it answers "where is the
    ligand relative to the fitted protein" -- and it rises both when the ligand
    slides out of a rigid pocket (a real failure) and when the protein itself
    relaxes, carrying a still-bound ligand with it (not a failure). @tt8804:
    "if they change tgt then they are fine".

    Fitting and measuring on the same CA group gives the protein's own
    displacement, so the two traces on one axis separate the cases: parallel
    rise means the complex drifted together, ligand-only rise means it left.

    Cheap and re-derivable: `whole.xtc` and `analysis.ndx` are persisted by
    `analyse()`, so this is one `gmx rms` over an existing trajectory -- no
    re-simulation. Idempotent, so it can be backfilled across a finished
    campaign. Returns None when the inputs are absent (an older run, or one that
    died before the PBC step) rather than raising, because a missing second
    trace should degrade the plot, not fail the report.
    """
    out = wd / PROT_RMSD_XVG
    tpr, whole, ndx = wd / "prod.tpr", wd / "whole.xtc", wd / "analysis.ndx"
    # THE CACHE IS KEYED ON CURRENCY, NOT ON EXISTENCE. This returned the file
    # whenever it merely existed, so after an adaptive extension -- which
    # rewrites `whole.xtc` over the longer trajectory -- the protein trace kept
    # covering the first 1.2 ns while the ligand trace beside it covered 5.2.
    # The plot then showed one line stopping a quarter of the way across and the
    # other continuing, which reads as the protein having been measured for less
    # time on purpose. Same shape as every pinned default in this project: right
    # when written, unable to announce that it no longer is.
    if out.is_file() and out.stat().st_size > 0:
        if not whole.is_file() or out.stat().st_mtime >= whole.stat().st_mtime:
            return out
        log.info("%s: protein RMSD is older than the trajectory; recomputing",
                 wd.name)
    if not all(p.is_file() and p.stat().st_size for p in (tpr, whole, ndx)):
        return None
    try:
        # Group 1 BOTH times: fit on CA, measure CA. Group 2 here would just
        # re-measure the ligand under a different filename.
        _gmx("rms", "-s", str(tpr), "-f", str(whole), "-n", str(ndx),
             "-o", PROT_RMSD_XVG, "-tu", "ns", cwd=wd,
             log_name="rms_protein.log", stdin="1\n1\n")
    except (AnalysisError, subprocess.SubprocessError):
        return None
    return out if out.is_file() and out.stat().st_size else None
