"""
Purpose: Prepare the shared Pin1 receptor (6VAJ) and both docking boxes.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: raw 6VAJ.pdb under immutable/inhibition/receptor/
Output: prepared_receptor.pdb/.pdbqt, box.json, box_expanded.json, prep_log.json

THE LOAD-BEARING ARTIFACT. All four approaches dock into the receptor this
module produces. Switching receptors — or preparing them differently per
approach — invalidates every cross-approach comparison downstream, so this runs
once and its outputs live in immutable/.

TWO BOXES, NOT ONE (adversary finding M5). The reference ligand QT7 is COVALENT
at Cys113, so a box drawn tightly around it is centred on the warhead
sub-pocket. That is correct for T_3/T_4, which attack exactly that atom. For the
non-covalent approaches it would bias generation and docking toward a sub-pocket
they have no reason to prefer, so T_1/T_2 get a box expanded to cover the full
PPIase pocket.
"""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
import subprocess
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path

import yaml

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = _REPO_ROOT / "config" / "receptor.yaml"

# Solvent and cryoprotectant heteroatoms carry no information for docking and
# would occupy pocket volume the ligand needs. Ions and cofactors are NOT in
# this list — stripping a structural metal silently changes the site.
STRIPPABLE_HET = {
    "HOH", "WAT",                              # water
    "GOL", "EDO", "MPD",                       # glycerol / glycols
    "PEG", "PG4", "PGE", "1PE", "P6G", "2PE",  # PEG fragments (6VAJ carries PG4)
    "SO4", "PO4", "NO3", "CL", "NA", "K",      # buffer ions
    "DMS", "ACT", "FMT", "TRS", "IMD",         # DMSO, acetate, formate, Tris, imidazole
}


class ReceptorPrepError(RuntimeError):
    """Receptor preparation failed or produced an unusable artifact."""


@dataclass
class Box:
    """A docking box: centre and size in Angstroms."""

    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    def to_dict(self) -> dict:
        return asdict(self)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_config(path: Path | str | None = None) -> dict:
    """Load config/receptor.yaml."""
    p = Path(path) if path else DEFAULT_CONFIG
    if not p.is_file():
        raise ReceptorPrepError(f"receptor config not found: {p}")
    return yaml.safe_load(p.read_text(encoding="utf-8"))


def parse_ligand_coords(pdb_path: Path, het_code: str) -> list[tuple[float, float, float]]:
    """Extract the reference ligand's atom coordinates from a PDB file.

    The ligand is stripped before docking, but its coordinates DEFINE the box —
    so they are read out first and retained, not discarded with the ligand.

    Parameters
    ----------
    pdb_path : Path
        Raw PDB file.
    het_code : str
        Three-letter HET code of the reference ligand (QT7 for 6VAJ).

    Returns
    -------
    list of tuple
        (x, y, z) per ligand atom.

    Raises
    ------
    ReceptorPrepError
        If the ligand is absent — that means the wrong structure was fetched.
    """
    coords: list[tuple[float, float, float]] = []
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("HETATM"):
            continue
        if line[17:20].strip() != het_code:
            continue
        # Fixed-column PDB format; slicing is correct here, splitting is not
        # (coordinates run together when a value reaches -100.000).
        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not coords:
        raise ReceptorPrepError(
            f"reference ligand {het_code!r} not found in {pdb_path.name} — "
            "wrong structure, or the ligand code changed upstream."
        )
    return coords


def centroid(coords: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    """Geometric centre of a coordinate list."""
    n = len(coords)
    return (
        sum(c[0] for c in coords) / n,
        sum(c[1] for c in coords) / n,
        sum(c[2] for c in coords) / n,
    )


def find_atom(pdb_path: Path, chain: str, resid: int, atom_name: str
              ) -> tuple[float, float, float] | None:
    """Locate a named protein atom (used for Cys113 SG)."""
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        if line[22:26].strip() != str(resid):
            continue
        if line[12:16].strip() != atom_name:
            continue
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return None


def strip_structure(src: Path, dst: Path, *, ligand_het: str) -> dict[str, int]:
    """Write a receptor-only PDB: reference ligand and solvent removed.

    Returns
    -------
    dict
        Counts of what was kept and dropped, for the prep log.
    """
    kept, dropped_lig, dropped_solvent, dropped_other_het = [], 0, 0, 0
    for line in src.read_text(encoding="utf-8").splitlines():
        if line.startswith("HETATM"):
            res = line[17:20].strip()
            if res == ligand_het:
                dropped_lig += 1
                continue
            if res in STRIPPABLE_HET:
                dropped_solvent += 1
                continue
            # An unrecognized heteroatom is kept and reported rather than
            # silently discarded — it may be a structural cofactor.
            dropped_other_het += 1
            kept.append(line)
            continue
        if line.startswith(("ATOM", "TER", "END")):
            kept.append(line)
    dst.write_text("\n".join(kept) + "\n", encoding="utf-8")
    return {
        "atoms_kept": len(kept),
        "ligand_atoms_removed": dropped_lig,
        "solvent_atoms_removed": dropped_solvent,
        "other_het_atoms_retained": dropped_other_het,
    }


def _run(cmd: list[str]) -> str:
    """Run a subprocess, raising with captured stderr on failure."""
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise ReceptorPrepError(
            f"command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr[:2000]}"
        )
    return proc.stdout


# `reduce` lives in the amber_md env, not cheminf. Searched rather than assumed
# so a caller in either env works.
_REDUCE_CANDIDATES = (
    "/data/lab_vm/envs/dwi_amber_md/bin/reduce",
    "reduce",
)


def _find_reduce() -> str:
    """Locate the `reduce` binary, or raise with what was tried."""
    for c in _REDUCE_CANDIDATES:
        if Path(c).is_file() or shutil.which(c):
            return c
    raise ReceptorPrepError(
        "`reduce` not found (tried: " + ", ".join(_REDUCE_CANDIDATES) + "). "
        "It ships with AmberTools; install the amber_md env."
    )


def protonate(src: Path, dst: Path, ph: float) -> None:
    """Add hydrogens with `reduce`, preserving residue identity.

    DO NOT use ``obabel -p`` here. It renumbers residues from 1, RENAMES them
    (LYS 6 became ALA 1), invents chain IDs, and silently drops residues — on
    6VAJ it discarded 28 of 150 including **Cys113**, the catalytic residue the
    entire covalent campaign targets. The output still looked like a protein and
    would have docked happily, producing plausible and meaningless scores.

    `reduce` (Word et al. 1999) is the Richardson-lab hydrogen placer that
    AmberTools ships. It optimizes His/Asn/Gln flips and OH rotamers and leaves
    the heavy-atom record alone.

    Parameters
    ----------
    src, dst : Path
        Input (stripped) and output (protonated) PDB.
    ph : float
        Recorded for provenance. `reduce -BUILD` protonates at physiological
        pH by its own rules rather than taking a pH argument, so this is not
        passed through — it is logged so the discrepancy is visible rather
        than implied.
    """
    exe = _find_reduce()
    # reduce writes the structure to stdout and its commentary to stderr;
    # a non-zero exit is common even on success, so validate the OUTPUT.
    proc = subprocess.run([exe, "-BUILD", str(src)], capture_output=True, text=True)
    if not proc.stdout.strip():
        raise ReceptorPrepError(
            f"reduce produced no structure:\n{proc.stderr[:1000]}")
    dst.write_text(proc.stdout, encoding="utf-8")
    log.info("protonated with reduce (nominal pH %.1f) -> %s", ph, dst.name)


def residue_map(path: Path) -> dict[tuple[str, str], str]:
    """{(chain, resid): resname} for every ATOM record — the preservation check."""
    out: dict[tuple[str, str], str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM"):
            out[(line[21], line[22:26].strip())] = line[17:20].strip()
    return out


def assert_preserved(reference: Path, produced: Path, *, catalytic: tuple[str, int, str],
                     label: str) -> dict:
    """Fail loudly if preparation altered the protein it was meant to protonate.

    These are POST-CONDITIONS, checked on the output. The original module
    verified the catalytic residue in the *input* and never re-checked the
    result, which is exactly how a corrupt receptor reached disk and got
    hash-pinned into a manifest.

    Raises
    ------
    ReceptorPrepError
        If residues were dropped, renamed, renumbered, or the catalytic residue
        is missing.
    """
    ref, got = residue_map(reference), residue_map(produced)
    chain, resid, resname = catalytic

    problems: list[str] = []
    if len(got) != len(ref):
        problems.append(f"residue count {len(got)} != reference {len(ref)} "
                        f"({len(ref) - len(got)} lost)")
    ref_chains, got_chains = sorted({k[0] for k in ref}), sorted({k[0] for k in got})
    if got_chains != ref_chains:
        problems.append(f"chain set {got_chains} != reference {ref_chains}")
    renamed = [f"{c}:{r} {ref[(c, r)]}->{got[(c, r)]}"
               for (c, r) in ref if (c, r) in got and ref[(c, r)] != got[(c, r)]]
    if renamed:
        problems.append(f"{len(renamed)} residue(s) renamed, e.g. {renamed[:3]}")
    if got.get((chain, str(resid))) != resname:
        problems.append(
            f"catalytic {resname}{resid} missing or renamed in {label} "
            f"(found {got.get((chain, str(resid)))!r})")

    if problems:
        raise ReceptorPrepError(
            f"{label} failed structure-preservation checks:\n  - "
            + "\n  - ".join(problems)
            + "\nThe prepared receptor is NOT usable; nothing should dock against it."
        )
    return {"residues": len(got), "chains": got_chains,
            "catalytic_present": f"{resname}{resid}"}


def to_pdbqt(src: Path, dst: Path) -> None:
    """Convert the prepared receptor to PDBQT for Vina/smina.

    ``-xr`` marks the receptor rigid, which is what the docking protocols in
    every approach assume.
    """
    _run(["obabel", str(src), "-O", str(dst), "-xr"])
    if not dst.is_file() or dst.stat().st_size == 0:
        raise ReceptorPrepError(f"pdbqt conversion produced no output at {dst}")


def prepare(config_path: Path | str | None = None, *, force: bool = False) -> dict:
    """Run the full receptor preparation and write every artifact.

    Parameters
    ----------
    config_path : Path or str, optional
        Path to config/receptor.yaml.
    force : bool, optional
        Overwrite existing prepared artifacts. Off by default: these live under
        immutable/ and every downstream result is tied to them, so replacing
        them silently would invalidate prior runs without a trace.

    Returns
    -------
    dict
        The prep log (also written to disk as JSON).
    """
    cfg = load_config(config_path)
    rc = cfg["receptor"]
    raw = Path(rc["raw_path"])
    if not raw.is_file():
        raise ReceptorPrepError(
            f"raw structure not found: {raw}. Fetch it first:\n"
            f"  curl -sS -o {raw} {rc['source_url']}"
        )

    out_pdb = Path(rc["preparation"]["outputs"]["prepared_pdb"])
    out_pdbqt = Path(rc["preparation"]["outputs"]["prepared_pdbqt"])
    out_log = Path(rc["preparation"]["outputs"]["prep_log"])
    for p in (out_pdb, out_pdbqt, out_log):
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.exists() and not force:
            raise ReceptorPrepError(
                f"{p} already exists. Downstream results are tied to it; pass "
                "force=True only if you intend to invalidate them."
            )

    het = rc["ligand"]["het_code"]
    lig_coords = parse_ligand_coords(raw, het)
    center = centroid(lig_coords)
    log.info("reference ligand %s: %d atoms, centroid %.3f %.3f %.3f",
             het, len(lig_coords), *center)

    cys = cfg.get("receptor", {}).get("catalytic_residue")
    if cys is None:
        # Authoritative source is choreography.yaml; receptor.yaml may omit it.
        chor = yaml.safe_load(
            (_REPO_ROOT / "config" / "choreography.yaml").read_text(encoding="utf-8"))
        cys = chor["target"]["catalytic_residue"]
    cys_sg = find_atom(raw, cys["chain"], cys["resid"], cys["atom"])
    if cys_sg is None:
        raise ReceptorPrepError(
            f"{cys['resname']}{cys['resid']} {cys['atom']} not found in the raw "
            "structure — the covalent approaches have nothing to attack."
        )

    stripped_tmp = out_pdb.with_suffix(".stripped.pdb")
    counts = strip_structure(raw, stripped_tmp, ligand_het=het)
    protonate(stripped_tmp, out_pdb, rc["preparation"]["protonation_ph"])
    to_pdbqt(out_pdb, out_pdbqt)

    # POST-CONDITIONS, checked on the OUTPUTS against the stripped input.
    # The original module verified the catalytic residue in the *input* and
    # never re-checked the result — which is exactly how a receptor missing
    # Cys113 reached disk and got hash-pinned into a manifest.
    catalytic = (cys["chain"], cys["resid"], cys["resname"])
    checks = {
        "prepared_pdb": assert_preserved(stripped_tmp, out_pdb,
                                         catalytic=catalytic, label="prepared_pdb"),
        "prepared_pdbqt": assert_preserved(stripped_tmp, out_pdbqt,
                                           catalytic=catalytic, label="prepared_pdbqt"),
    }
    log.info("post-conditions passed: %s", checks)
    stripped_tmp.unlink(missing_ok=True)

    boxes = {}
    for name, spec in cfg["boxes"].items():
        size = spec["size"]
        box = Box(center[0], center[1], center[2], size[0], size[1], size[2])
        bpath = Path(spec["out_path"])
        payload = {
            **box.to_dict(),
            "name": name,
            "used_by": spec["used_by"],
            "derived_from": f"{het} centroid in {raw.name}",
        }
        bpath.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        boxes[name] = payload
        log.info("wrote %s box -> %s", name, bpath)

    prep_log = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pdb_id": rc["pdb_id"],
        "raw_sha256": _sha256(raw),
        "prepared_pdb_sha256": _sha256(out_pdb),
        "prepared_pdbqt_sha256": _sha256(out_pdbqt),
        "reference_ligand": {"het_code": het, "n_atoms": len(lig_coords),
                             "centroid": list(center)},
        "cys113_sg": list(cys_sg),
        "protonation_ph": rc["preparation"]["protonation_ph"],
        "structure_counts": counts,
        "preservation_checks": checks,
        "protonation_tool": "reduce -BUILD (AmberTools; Word et al. 1999)",
        "boxes": boxes,
    }
    out_log.write_text(json.dumps(prep_log, indent=2) + "\n", encoding="utf-8")
    log.info("receptor prep complete -> %s", out_log)
    return prep_log


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    import sys

    print(json.dumps(prepare(force="--force" in sys.argv), indent=2))
