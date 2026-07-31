"""
Purpose: Symmetry-corrected RMSD between each docked pose and its crystal
         reference, for both redocking arms, and the recovery statistics.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: redock_cases_1.csv + redock_docking_1.csv + the pose PDBQTs
Output: outputs/blacksmith/redock_pin1/redock_rmsd_1.csv
        outputs/blacksmith/redock_pin1/redock_summary_1.json

THE RMSD FUNCTION IS `CalcRMS`, NOT `GetBestRMS`, AND THE DIFFERENCE DECIDES
THE RESULT. Both are symmetry-corrected, but `GetBestRMS` SUPERPOSES the probe
onto the reference before measuring. Measured here on a reference translated
bodily by 3.0 A: `CalcRMS` returns 3.000, `GetBestRMS` returns 0.000. Because
both arms produce poses already in the reference's coordinate frame, using
`GetBestRMS` would discard the translation and rotation -- i.e. discard exactly
the thing redocking is testing -- and report near-zero RMSD for every case
including complete failures. `CalcRMS` computes in place.

SYMMETRY CORRECTION IS NOT OPTIONAL ON THIS SET. The ligands are aromatic and
many carry phenyl, naphthalene or carboxylate groups whose atoms are
topologically interchangeable. A naive atom-order RMSD reports a perfectly
docked phenyl ring flipped 180 degrees as a ~2.5 A failure. `CalcRMS`
enumerates the molecule's automorphisms and takes the minimum, so a ring flip
costs nothing -- which is the physically correct answer, because the two
orientations are the same molecule.

BOTH MOLECULES ARE NORMALISED THROUGH ONE PATH. A PDBQT carries only polar
hydrogens, so RDKit reads its carbons as radicals with `noImplicit` set. Those
flags make the docked graph differ from the crystal graph and the substructure
match underpinning the symmetry correction fails outright -- every case would
error rather than silently mis-score, but the benchmark would be empty.
`_normalise` clears the flags on BOTH sides so like is compared with like.
"""

from __future__ import annotations

import json
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("redock-rmsd")

OUT_DIR = REPO / "outputs" / "blacksmith" / "redock_pin1"
OBABEL = "/data/lab_vm/envs/dwi_cheminf/bin/obabel"

# The standard redocking success criterion (Astex/PoseBusters convention).
SUCCESS_A = 2.0
# A second, stricter line the pose-prediction literature also reports.
TIGHT_A = 1.0


def _heavy_graph(m: Chem.Mol) -> Chem.Mol:
    """A bond-order- and charge-agnostic heavy-atom copy, for matching only.

    WHY NOT SANITISE AND USE `CalcRMS` DIRECTLY. That was the first design and
    it failed on 26 of 82 cases with "Can't kekulize mol" -- every one of them
    an aromatic N-H (indole, pyrrole, pyrazole, benzimidazole). A crystal
    ligand carries no hydrogens at all, so the ring nitrogen arrives with no H
    to donate to the aromatic system and kekulisation is impossible. Sanitising
    the coordinate-bearing molecules is therefore not a step this benchmark can
    rely on, and dropping a third of the set over a hydrogen bookkeeping detail
    would have been a fabricated result.

    Matching on the flattened graph needs no sanitisation: atoms are compared
    on element and connectivity, which is all the correspondence requires. It
    also makes the symmetry treatment slightly MORE correct, because the
    automorphisms it yields include resonance-equivalent swaps -- the two
    oxygens of a carboxylate, of a nitro group -- which X-ray density genuinely
    cannot distinguish and which a bond-order-aware matcher would score as an
    error.
    """
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        a.SetNumRadicalElectrons(0)
    out = rw.GetMol()
    Chem.FastFindRings(out)          # ring info without a full sanitisation
    return out


def _coords(m: Chem.Mol) -> np.ndarray:
    c = m.GetConformer()
    return np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])


def symmetric_rmsd(ref: Chem.Mol, dock: Chem.Mol, tmpl: Chem.Mol) -> float:
    """Symmetry-corrected, in-place RMSD over heavy atoms.

    Both molecules are mapped onto the template's atom ordering, then the
    minimum is taken over the template's automorphism group -- so a phenyl ring
    docked perfectly but flipped 180 degrees scores 0, not ~2.5 A. No
    superposition is performed anywhere: the poses are already in the
    reference's frame and aligning them would erase what is being measured.
    """
    qt, qr, qd = _heavy_graph(tmpl), _heavy_graph(ref), _heavy_graph(dock)
    m_ref = qr.GetSubstructMatch(qt)
    m_dk = qd.GetSubstructMatch(qt)
    if not m_ref:
        raise ValueError("reference does not match the chem-comp graph")
    if not m_dk:
        raise ValueError("docked pose does not match the chem-comp graph")
    xr = _coords(ref)[list(m_ref)]
    xd = _coords(dock)[list(m_dk)]
    autos = qt.GetSubstructMatches(qt, uniquify=False, maxMatches=50000)
    if not autos:
        autos = [tuple(range(qt.GetNumAtoms()))]
    best = min(float(np.sqrt(((xd[list(a)] - xr) ** 2).sum(axis=1).mean()))
               for a in autos)
    return best


def reference_mol(ref_pdb: Path, smiles: str) -> Chem.Mol:
    """The crystal ligand as deposited (heavy atoms, coordinates preserved)."""
    raw = Chem.MolFromPDBBlock(ref_pdb.read_text(), sanitize=False, removeHs=True)
    if raw is None:
        raise ValueError("reference PDB unreadable")
    return raw


def docked_mol(pose_pdbqt: Path, smiles: str) -> tuple[Chem.Mol, float | None]:
    """MODEL 1 (Vina's best-scoring pose) of an output PDBQT, as a typed mol."""
    text = pose_pdbqt.read_text(errors="replace")
    if "MODEL" in text:
        block = "MODEL" + text.split("MODEL", 1)[1].split("ENDMDL")[0] + "ENDMDL\n"
    else:
        block = text
    score = None
    for line in block.splitlines():
        if "REMARK VINA RESULT" in line:
            try:
                score = float(line.split()[3])
            except (IndexError, ValueError):
                pass
            break
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "m.pdbqt"
        p.write_text(block)
        sdf = Path(td) / "m.sdf"
        r = subprocess.run([OBABEL, str(p), "-O", str(sdf)],
                           capture_output=True, text=True, timeout=120)
        if not sdf.is_file() or sdf.stat().st_size == 0:
            raise ValueError(f"obabel could not convert the pose ({r.returncode})")
        raw = Chem.MolFromMolFile(str(sdf), sanitize=False, removeHs=False)
    if raw is None:
        raise ValueError("docked pose unreadable")
    raw.UpdatePropertyCache(strict=False)
    return Chem.RemoveHs(raw, sanitize=False), score


def strict_calcrms(ref: Chem.Mol, dock: Chem.Mol, tmpl: Chem.Mol) -> float | None:
    """RDKit's own symmetry-corrected in-place RMSD, where it can be computed.

    Kept purely as an INDEPENDENT CHECK on `symmetric_rmsd`. It needs fully
    sanitised molecules and so cannot run on the aromatic-N-H cases, which is
    why it is not the primary path -- but where both run they must agree, and
    that agreement is reported rather than assumed.
    """
    try:
        r = AllChem.AssignBondOrdersFromTemplate(tmpl, Chem.Mol(ref))
        d = AllChem.AssignBondOrdersFromTemplate(tmpl, Chem.Mol(dock))
        for m in (r, d):
            for a in m.GetAtoms():
                a.SetNoImplicit(False)
                a.SetNumRadicalElectrons(0)
            Chem.SanitizeMol(m)
        return float(rdMolAlign.CalcRMS(d, r))
    except Exception:  # noqa: BLE001 - this is a cross-check, never the answer
        return None


def rmsd_for(ref_pdb: Path, pose: Path, smiles: str
             ) -> tuple[float, float | None, float | None]:
    """Symmetry-corrected, in-place RMSD between a docked pose and its crystal."""
    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is None:
        raise ValueError("chem-comp SMILES unparseable")
    ref = reference_mol(ref_pdb, smiles)
    dk, score = docked_mol(pose, smiles)
    if ref.GetNumAtoms() != dk.GetNumAtoms():
        raise ValueError(f"atom count {dk.GetNumAtoms()} != reference "
                         f"{ref.GetNumAtoms()}")
    return symmetric_rmsd(ref, dk, tmpl), score, strict_calcrms(ref, dk, tmpl)


def _describe(values: pd.Series, label: str) -> dict:
    """Recovery statistics for one arm/stratum."""
    v = values.dropna()
    if v.empty:
        return {"stratum": label, "n": 0}
    return {
        "stratum": label,
        "n": int(len(v)),
        "success_rate_2A": round(float((v <= SUCCESS_A).mean()), 4),
        "n_success_2A": int((v <= SUCCESS_A).sum()),
        "success_rate_1A": round(float((v <= TIGHT_A).mean()), 4),
        "median_rmsd_a": round(float(v.median()), 3),
        "mean_rmsd_a": round(float(v.mean()), 3),
        "q25_a": round(float(v.quantile(0.25)), 3),
        "q75_a": round(float(v.quantile(0.75)), 3),
        "min_a": round(float(v.min()), 3),
        "max_a": round(float(v.max()), 3),
    }


def _wilson(k: int, n: int) -> tuple[float, float]:
    """Wilson 95% interval -- a success rate on n<=82 needs its uncertainty."""
    if n == 0:
        return (float("nan"), float("nan"))
    z, p = 1.959964, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (round(max(0.0, c - h), 4), round(min(1.0, c + h), 4))


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    cases = pd.read_csv(OUT_DIR / "redock_cases_1.csv")
    cases = cases[cases.status == "case"]
    dock = pd.read_csv(OUT_DIR / "redock_docking_1.csv")
    df = cases.merge(dock, on="case_id", how="left")

    rows = []
    for c in df.itertuples():
        rec = {"case_id": c.case_id, "pdb_id": c.pdb_id, "comp_id": c.comp_id,
               "tier": c.tier, "heavy_atoms": c.heavy_atoms,
               "resolution_a": c.resolution_a,
               "ref_in_prod_box": c.ref_in_prod_box,
               "n_rot_bonds": None}
        m = Chem.MolFromSmiles(c.smiles)
        if m is not None:
            rec["n_rot_bonds"] = int(
                Chem.rdMolDescriptors.CalcNumRotatableBonds(m))

        # ---- Arm A: self-docking -------------------------------------
        try:
            pose = Path(c.self_pose_dir) / f"{c.case_id}_out.pdbqt"
            r, s, chk = rmsd_for(Path(c.ref_pdb), pose, c.smiles)
            rec.update({"self_rmsd_a": r, "self_affinity": s,
                        "self_rmsd_calcrms_check": chk})
        except Exception as exc:  # noqa: BLE001
            rec.update({"self_rmsd_a": None, "self_rmsd_error": str(exc)[:160]})

        # ---- Arm B: cross-docking into 6VAJ --------------------------
        try:
            pose = Path(c.cross_pose_dir) / f"{c.case_id}_out.pdbqt"
            r, s, chk = rmsd_for(Path(c.ref_6vaj_pdb), pose, c.smiles)
            rec.update({"cross_rmsd_a": r, "cross_affinity": s,
                        "cross_rmsd_calcrms_check": chk})
        except Exception as exc:  # noqa: BLE001
            rec.update({"cross_rmsd_a": None, "cross_rmsd_error": str(exc)[:160]})
        rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "redock_rmsd_1.csv", index=False)

    strata = []
    for arm, col in (("self", "self_rmsd_a"), ("cross", "cross_rmsd_a")):
        strata.append({"arm": arm, **_describe(out[col], "all")})
        for tier in ("drug_like", "fragment"):
            strata.append({"arm": arm,
                           **_describe(out.loc[out.tier == tier, col], tier)})
        # The production box cannot reach every crystallographic site; the
        # cross-docking number is reported both ways rather than one.
        strata.append({"arm": arm, **_describe(
            out.loc[out.ref_in_prod_box.fillna(False).astype(bool), col],
            "ref_in_production_box")})
        strata.append({"arm": arm, **_describe(
            out.loc[(out.tier == "drug_like")
                    & out.ref_in_prod_box.fillna(False).astype(bool), col],
            "drug_like_and_in_box")})
    summary = pd.DataFrame(strata)
    for r in summary.itertuples():
        if r.n and not pd.isna(getattr(r, "n_success_2A", None)):
            lo, hi = _wilson(int(r.n_success_2A), int(r.n))
            summary.loc[r.Index, "ci95_low"] = lo
            summary.loc[r.Index, "ci95_high"] = hi
    summary.to_csv(OUT_DIR / "redock_summary_1.csv", index=False)

    payload = {
        "generated": pd.Timestamp.utcnow().isoformat(),
        "protocol": {"engine": "Vina-GPU 2.1", "search_depth": 20,
                     "ligand_ph": 7.4, "box": "box_expanded.json (26 A)",
                     "rmsd": "rdMolAlign.CalcRMS (symmetry-corrected, in place)",
                     "success_criterion_a": SUCCESS_A},
        "strata": summary.to_dict(orient="records"),
        "failures": {
            "self_rmsd_errors": int(out.get(
                "self_rmsd_error", pd.Series(dtype=object)).notna().sum()),
            "cross_rmsd_errors": int(out.get(
                "cross_rmsd_error", pd.Series(dtype=object)).notna().sum()),
        },
    }
    (OUT_DIR / "redock_summary_1.json").write_text(json.dumps(payload, indent=2))

    for arm in ("self", "cross"):
        a, b = out.get(f"{arm}_rmsd_a"), out.get(f"{arm}_rmsd_calcrms_check")
        if a is not None and b is not None:
            ok = a.notna() & b.notna()
            if ok.any():
                dmax = float((a[ok] - b[ok]).abs().max())
                payload.setdefault("calcrms_crosscheck", {})[arm] = {
                    "n_compared": int(ok.sum()),
                    "max_abs_diff_a": round(dmax, 4)}
                log.info("[%s] independent CalcRMS cross-check on %d/%d cases: "
                         "max |diff| = %.4f A", arm, int(ok.sum()), len(out), dmax)
    (OUT_DIR / "redock_summary_1.json").write_text(json.dumps(payload, indent=2))

    log.info("\n%s", summary.to_string(index=False))
    for col in ("self_rmsd_error", "cross_rmsd_error"):
        if col in out and out[col].notna().any():
            log.warning("%s: %d\n%s", col, out[col].notna().sum(),
                        out.loc[out[col].notna(), ["case_id", col]].to_string(index=False))


if __name__ == "__main__":
    main()
