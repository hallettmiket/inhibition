"""
Purpose: Ensemble MM-GBSA -- three-leg dG per MD frame, with an honest error bar.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: an MM-GBSA workdir (three prmtops) + md/traj_nm.npy from the MD tier
Output: per-frame dG, mean, and a SEM widened by the statistical inefficiency

WHY THIS IS NOT A REFINEMENT. Every dG this project reports comes from ONE
minimised structure, so it has no variance and the "+/- uncertainty" column the
plan asks for could only ever have been fabricated. This module produces the
distribution that number needs.

SINGLE TRAJECTORY, AND WHAT THAT CHANGES. All three legs are sampled from the
COMPLEX frames. The existing dG values come from three INDEPENDENT
minimisations. Those are different estimators, not two precisions of one, so
the ensemble dG REPLACES the single-structure dG rather than refining it, and
putting both in one table invites exactly the wrong comparison. Single-
trajectory cancels receptor internal energy almost exactly and therefore has
far lower variance, at the cost of ignoring conformational reorganisation on
binding -- the standard trade, recorded so it is a choice and not an accident.

SANDER, NOT OPENMM. OpenMM generates the trajectory; it does not score it.
Measured on identical minimised structures, OpenMM's three-leg dG diverges from
sander's by +0.34 and +0.89 kcal/mol on two systems but +12.20 on act_000
(Sulfopin, the one active D0032 rests on) and -55.34 on another. On a dG scale
of ~-15 kcal/mol that is a different answer, not a correction. Scoring with
sander keeps every number on the same engine as the rest of the project and
removes the question entirely.

THE MAPPING IS VERIFIED, NOT ASSUMED. The link-atom scheme caps each leg with a
hydrogen that does not exist in the complex, so slicing complex frames is not a
slice. Worse, the cap is not appended: in the ligand leg it is INSERTED at
index 39, so the atom order diverges from the complex partway through and a
positional mapping would silently corrupt every atom after it. Atoms are
therefore mapped by (residue index, atom name), and the map is checked --
element identity at every mapped pair, exactly one unmapped hydrogen per leg,
counts reconciling -- with a hard failure on any violation. Six defects in this
project (D0025, D0028, D0029, D0030, D0033, D0034) came from a value being
derived from a name or a default rather than from the thing itself.
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from shared import mmgbsa as mg

log = logging.getLogger(__name__)

NM_TO_ANGSTROM = 10.0

# Single-point energies on a supplied trajectory. imin=5 makes sander read the
# whole mdcrd and emit one FINAL RESULTS block per frame, so a leg costs ONE
# invocation instead of one per frame (measured: 100 frames in 37 s, 1 core).
SINGLE_POINT_INPUT = f"""single-point energies on each trajectory frame
 &cntrl
  imin=5, maxcyc=1, ntb=0, igb={mg.IGB}, gbsa=1, cut=999.0,
 /
"""


class EnsembleError(RuntimeError):
    """The ensemble could not be built or trusted."""


@dataclass
class LegMap:
    """How one leg's atoms correspond to atoms of the complex.

    `complex_index[i]` is the complex atom that supplies leg atom i, except at
    `cap_index` where the atom does not exist in the complex and is placed
    geometrically instead.
    """

    leg: str
    complex_index: np.ndarray      # int, len == leg atom count; -1 at the cap
    cap_index: int                 # index WITHIN the leg
    cap_anchor_complex: int        # complex atom the cap bonds to
    cap_direction_complex: int     # complex atom that was severed (gives direction)
    cap_bond_length_a: float       # equilibrium length from the leg's own topology

    @property
    def n_atoms(self) -> int:
        return int(self.complex_index.size)


def _load(workdir: Path):
    try:
        import parmed
    except ImportError as exc:  # noqa: BLE001
        raise EnsembleError(
            "parmed is required; run under /data/lab_vm/envs/dwi_amber_md"
        ) from exc
    tops = {}
    for leg in ("complex", "receptor", "ligand"):
        p = workdir / f"{leg}.prmtop"
        if not p.is_file() or p.stat().st_size == 0:
            raise EnsembleError(f"{workdir.name}: no usable {leg}.prmtop")
        tops[leg] = parmed.load_file(str(p))
    return tops


# Both builders in this project name the ligand residue differently:
# shared/mmgbsa.py (covalent) writes LIG, shared/mmgbsa_noncovalent.py writes
# MOL. The name is accepted from either, but nothing is DECIDED from it -- see
# detect_scheme, which reads the topology instead.
LIGAND_RESNAMES = {mg.LIGAND_RESNAME, "MOL"}


def detect_scheme(complex_top, receptor_top, ligand_top) -> str:
    """'covalent' or 'non_covalent', decided from structure, not from a label.

    T_1 and T_2 dock non-covalently: the ligand is a separate molecule, so no
    bond is severed and no cap hydrogen exists. T_3, T_4 and the gate form a
    covalent adduct and need the link-atom construction. Applying the wrong one
    silently misplaces an atom or drops a real bond.

    The two cases are distinguished by two INDEPENDENT counts that must agree:
    the number of bonds crossing the ligand boundary (1 vs 0) and the atom-count
    excess (2 vs 0). Requiring both to agree means a topology that is neither
    shape fails loudly rather than being forced into whichever branch its
    residue happened to be named after.
    """
    lig = complex_top.residues[-1]
    if lig.name not in LIGAND_RESNAMES:
        raise EnsembleError(
            f"last residue is {lig.name!r}, expected one of "
            f"{sorted(LIGAND_RESNAMES)}")
    lig_idx = {a.idx for a in lig.atoms}
    crossing = sum(1 for b in complex_top.bonds
                   if (b.atom1.idx in lig_idx) != (b.atom2.idx in lig_idx))
    excess = (len(receptor_top.atoms) + len(ligand_top.atoms)
              - len(complex_top.atoms))
    if crossing == 1 and excess == 2:
        return "covalent"
    if crossing == 0 and excess == 0:
        return "non_covalent"
    raise EnsembleError(
        f"cannot classify this decomposition: {crossing} bond(s) cross the "
        f"ligand boundary and the atom-count excess is {excess}. A covalent "
        "link-atom split is (1, 2) and a non-covalent split is (0, 0); "
        "anything else is not a decomposition this module can score.")


def find_junction(complex_top) -> tuple:
    """The single bond joining the ligand residue to the protein.

    Located from the BOND LIST, not from a stored residue index or a SMARTS
    pattern. Both of those have produced silent misidentification in this
    project; a bond that crosses the ligand boundary is the thing itself.
    """
    lig = complex_top.residues[-1]
    if lig.name not in LIGAND_RESNAMES:
        raise EnsembleError(
            f"last residue is {lig.name!r}, expected one of "
            f"{sorted(LIGAND_RESNAMES)}")
    lig_idx = {a.idx for a in lig.atoms}
    crossing = [b for b in complex_top.bonds
                if (b.atom1.idx in lig_idx) != (b.atom2.idx in lig_idx)]
    if len(crossing) != 1:
        raise EnsembleError(
            f"expected exactly 1 bond crossing the ligand boundary, found "
            f"{len(crossing)}. The link-atom scheme assumes a single covalent "
            "attachment and cannot decompose anything else.")
    b = crossing[0]
    prot, ligand = (b.atom1, b.atom2) if b.atom1.idx not in lig_idx else (
        b.atom2, b.atom1)
    if prot.element_name != "S":
        raise EnsembleError(
            f"junction protein atom is {prot.element_name}, expected S "
            "(Cys113 SG)")
    return prot, ligand


def _cap_bond_length(leg_top, cap_atom) -> float:
    """Equilibrium bond length for the cap, from the leg's OWN parameters.

    Read rather than hardcoded: a C-H and an S-H cap differ by ~0.25 A, and the
    force field already states both. A literal here would be a second source of
    truth that could drift from the topology actually being scored.
    """
    bonds = [b for b in leg_top.bonds
             if cap_atom in (b.atom1, b.atom2)]
    if len(bonds) != 1:
        raise EnsembleError(
            f"cap atom {cap_atom.name} has {len(bonds)} bonds, expected 1")
    req = getattr(bonds[0].type, "req", None)
    if not req:
        raise EnsembleError(f"cap bond for {cap_atom.name} has no equilibrium "
                            "length in the topology")
    return float(req)


def build_leg_maps(workdir: Path) -> dict[str, LegMap]:
    """Verified correspondence from each leg's atoms to the complex's.

    Every check here exists because its absence would produce a plausible
    number rather than an error.
    """
    tops = _load(workdir)
    cx, rec, lig_top = tops["complex"], tops["receptor"], tops["ligand"]
    scheme = detect_scheme(cx, rec, lig_top)
    n_c, n_r, n_l = len(cx.atoms), len(rec.atoms), len(lig_top.atoms)

    if scheme == "covalent":
        prot_atom, lig_atom = find_junction(cx)
        expect_unmapped = 1
    else:
        prot_atom = lig_atom = None
        expect_unmapped = 0

    lig_res = cx.residues[-1]
    maps: dict[str, LegMap] = {}

    def _cap_fields(leg_top, unmapped: list[int], anchor, direction):
        """Cap descriptor for a leg, or the null descriptor when uncapped."""
        if expect_unmapped == 0:
            return -1, -1, -1, 0.0
        cap = leg_top.atoms[unmapped[0]]
        if cap.element_name != "H":
            raise EnsembleError(f"cap is {cap.element_name}, expected H")
        return (unmapped[0], anchor.idx, direction.idx,
                _cap_bond_length(leg_top, cap))

    # --- ligand leg: match on atom name within the single ligand residue ---
    by_name = {a.name: a for a in lig_res.atoms}
    idx = np.full(n_l, -1, dtype=int)
    unmapped = []
    for i, a in enumerate(lig_top.atoms):
        src = by_name.get(a.name)
        if src is None:
            unmapped.append(i)
            continue
        if src.element_name != a.element_name:
            raise EnsembleError(
                f"ligand leg atom {a.name} is {a.element_name} but the complex "
                f"has {src.element_name} under that name")
        idx[i] = src.idx
    if len(unmapped) != expect_unmapped:
        raise EnsembleError(
            f"ligand leg has {len(unmapped)} atoms absent from the complex, "
            f"expected {expect_unmapped} for a {scheme} decomposition")
    ci, ca, cd, cl = _cap_fields(lig_top, unmapped, lig_atom, prot_atom)
    maps["ligand"] = LegMap(leg="ligand", complex_index=idx, cap_index=ci,
                            cap_anchor_complex=ca, cap_direction_complex=cd,
                            cap_bond_length_a=cl)

    # --- receptor leg: match on (residue index, atom name) ---
    prot_by_key = {(a.residue.idx, a.name): a
                   for a in cx.atoms if a.residue.idx != lig_res.idx}
    idx = np.full(n_r, -1, dtype=int)
    unmapped = []
    for i, a in enumerate(rec.atoms):
        src = prot_by_key.get((a.residue.idx, a.name))
        if src is None:
            unmapped.append(i)
            continue
        if src.element_name != a.element_name:
            raise EnsembleError(
                f"receptor leg atom {a.residue.idx}:{a.name} is "
                f"{a.element_name} but the complex has {src.element_name}")
        idx[i] = src.idx
    if len(unmapped) != expect_unmapped:
        raise EnsembleError(
            f"receptor leg has {len(unmapped)} atoms absent from the complex, "
            f"expected {expect_unmapped} for a {scheme} decomposition")
    ci, ca, cd, cl = _cap_fields(rec, unmapped, prot_atom, lig_atom)
    maps["receptor"] = LegMap(leg="receptor", complex_index=idx, cap_index=ci,
                              cap_anchor_complex=ca, cap_direction_complex=cd,
                              cap_bond_length_a=cl)

    # --- complex leg: identity, but expressed the same way so callers do not
    #     special-case it and quietly skip its checks ---
    maps["complex"] = LegMap(
        leg="complex", complex_index=np.arange(n_c, dtype=int),
        cap_index=-1, cap_anchor_complex=-1, cap_direction_complex=-1,
        cap_bond_length_a=0.0)

    log.debug("%s: maps built (cap at ligand[%d], receptor[%d])",
              workdir.name, maps["ligand"].cap_index,
              maps["receptor"].cap_index)
    return maps


def slice_frame(frame_a: np.ndarray, m: LegMap) -> np.ndarray:
    """One leg's coordinates for one complex frame, in Angstrom.

    The cap hydrogen is placed along the severed bond at the force field's own
    equilibrium length -- the standard link-atom construction.
    """
    out = np.empty((m.n_atoms, 3), dtype=float)
    mapped = m.complex_index >= 0
    out[mapped] = frame_a[m.complex_index[mapped]]
    if m.cap_index >= 0:
        anchor = frame_a[m.cap_anchor_complex]
        direction = frame_a[m.cap_direction_complex] - anchor
        norm = np.linalg.norm(direction)
        if norm < 1e-6:
            raise EnsembleError("severed-bond atoms are coincident; cannot "
                                "orient the cap hydrogen")
        out[m.cap_index] = anchor + direction / norm * m.cap_bond_length_a
    return out


def write_mdcrd(path: Path, frames_a: np.ndarray, title: str = "ensemble") -> None:
    """Amber ASCII trajectory (10F8.3).

    Coordinates are truncated to 0.001 A by the format, worth ~0.1 kcal/mol per
    leg here -- well below the ensemble spread this module exists to measure,
    and the format sander reads without further dependencies.
    """
    lines = [title]
    for frame in frames_a:
        flat = frame.reshape(-1)
        for i in range(0, flat.size, 10):
            lines.append("".join(f"{v:8.3f}" for v in flat[i:i + 10]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def score_leg(workdir: Path, leg: str, mdcrd: Path, out_dir: Path
              ) -> list[float]:
    """Per-frame total energy for one leg, via a single sander invocation."""
    inp = out_dir / f"{leg}_sp.in"
    inp.write_text(SINGLE_POINT_INPUT, encoding="utf-8")
    out = out_dir / f"{leg}_sp.out"
    cmd = [str(mg.AMBER_ENV / "bin" / "sander"), "-O",
           "-i", str(inp), "-p", str(workdir / f"{leg}.prmtop"),
           "-c", str(workdir / f"{leg}.inpcrd"), "-y", str(mdcrd),
           "-o", str(out)]
    proc = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True,
                          timeout=7200)
    if not out.is_file():
        raise EnsembleError(f"{leg}: sander produced no output "
                            f"({proc.stderr[:200]})")
    text = out.read_text(errors="ignore")
    blocks = text.split("FINAL RESULTS")[1:]
    if not blocks:
        raise EnsembleError(f"{leg}: sander emitted no FINAL RESULTS block")
    totals = []
    for b in blocks:
        terms = {}
        for k, v in mg.ENERGY_LINE_RX.findall(b):
            terms.setdefault(k.strip(), float(v))
        totals.append(mg.LegEnergies(terms=terms).total)
    return totals


def statistical_inefficiency(series: np.ndarray) -> float:
    """g = 1 + 2*sum(C(t)), frames per independent sample.

    The SEM is widened by sqrt(g). Without it a correlated trajectory reports an
    error bar too small by precisely the factor that makes reporting it
    worthwhile.
    """
    x = np.asarray(series, dtype=float)
    n = x.size
    if n < 4:
        return float(max(1, n))
    x = x - x.mean()
    var = x.var()
    if var <= 0:
        return 1.0
    g, t = 1.0, 1
    while t < n - 1:
        c = float((x[: n - t] * x[t:]).mean() / var)
        if c <= 0:
            break
        g += 2.0 * c * (1.0 - t / n)
        t += 1
    return max(1.0, g)


def ensemble_dg(workdir: Path, traj_nm: np.ndarray | None = None,
                out_dir: Path | None = None, discard_first: int = 0) -> dict:
    """Ensemble-averaged dG with an uncertainty that reflects correlation."""
    out_dir = out_dir or (workdir / "ensemble")
    out_dir.mkdir(parents=True, exist_ok=True)

    if traj_nm is None:
        p = workdir / "md" / "traj_nm.npy"
        if not p.is_file():
            raise EnsembleError(f"no trajectory at {p}; run the MD tier first")
        traj_nm = np.load(p)
    if discard_first:
        traj_nm = traj_nm[discard_first:]
    if traj_nm.shape[0] < 4:
        raise EnsembleError(f"{traj_nm.shape[0]} frames is too few for an "
                            "ensemble; refusing to report a mean and an SEM "
                            "from it")

    frames_a = np.asarray(traj_nm, dtype=float) * NM_TO_ANGSTROM
    maps = build_leg_maps(workdir)
    tops = _load(workdir)
    scheme = detect_scheme(tops["complex"], tops["receptor"], tops["ligand"])

    per_leg: dict[str, list[float]] = {}
    for leg in ("complex", "receptor", "ligand"):
        m = maps[leg]
        sliced = np.stack([slice_frame(f, m) for f in frames_a])
        mdcrd = out_dir / f"{leg}.mdcrd"
        write_mdcrd(mdcrd, sliced, title=f"{workdir.name} {leg}")
        totals = score_leg(workdir, leg, mdcrd, out_dir)
        if len(totals) != frames_a.shape[0]:
            raise EnsembleError(
                f"{leg}: sander returned {len(totals)} energies for "
                f"{frames_a.shape[0]} frames")
        per_leg[leg] = totals

    dg = (np.asarray(per_leg["complex"]) - np.asarray(per_leg["receptor"])
          - np.asarray(per_leg["ligand"]))

    # PER-FRAME VALUES ARE SAVED, not just the summary. Two reasons. First, the
    # summary cannot be re-interrogated: questions like "how much of the
    # single-structure/ensemble disagreement is sampling rather than the change
    # of estimator" need the distribution, and re-deriving it means re-running
    # sander over every frame. Second, a mean and an SEM conceal shape -- a
    # bimodal dG from a ligand flipping between two poses reports the same two
    # numbers as a tight unimodal one.
    np.save(out_dir / "dg_per_frame.npy", dg.astype(np.float32))
    np.savez(out_dir / "legs_per_frame.npz",
             **{k: np.asarray(v, dtype=np.float32) for k, v in per_leg.items()})

    g = statistical_inefficiency(dg)
    n_eff = max(1.0, dg.size / g)
    sem = float(dg.std(ddof=1) / np.sqrt(n_eff)) if dg.size > 1 else float("nan")

    return {
        "dG_kcal_ensemble": round(float(dg.mean()), 3),
        "dG_sem_kcal": round(sem, 3),
        "dG_sd_kcal": round(float(dg.std(ddof=1)), 3),
        "dG_min_kcal": round(float(dg.min()), 3),
        "dG_max_kcal": round(float(dg.max()), 3),
        "n_frames": int(dg.size),
        "n_effective_samples": round(float(n_eff), 2),
        "statistical_inefficiency": round(float(g), 2),
        "frames_discarded": int(discard_first),
        "decomposition": scheme,
        "scheme": ("single-trajectory 3-leg, sander single points"
                   + (", link-atom cap at Cys113 SG-C" if scheme == "covalent"
                      else ", no bond severed (non-covalent)")),
        "engine": "sander (NOT OpenMM -- see module docstring)",
        "igb": mg.IGB,
        "pb_radii": mg.PB_RADII,
        "ensemble_averaged": True,
        "replaces_not_refines": "single-structure dG uses three INDEPENDENT "
                                "minimisations; this uses one trajectory. "
                                "Different estimators, not two precisions.",
    }
