"""
Purpose: The single source of RDKit physicochemical descriptors + QED + SAscore.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: canonical SMILES
Output: a fixed descriptor dict per molecule; a DataFrame for a batch

WHY THIS IS THE ONLY PLACE DESCRIPTORS ARE COMPUTED: the physicochemical axes
(MW, cLogP, TPSA, QED, HBD, HBA, rotatable bonds, fraction sp3, SAscore) are the
one set of numbers directly comparable across all four approaches, because they
come from the same tool with the same definition. That comparability is what the
integration phase leans on given there is deliberately no authoritative
cross-approach numeric join. If T_1 and T_4 each computed "logP" their own way,
the pooled plots in the GUI would be quietly wrong. So every approach imports
this module, and no approach computes these itself.

SAscore is an RDKit *Contrib* script, not core API. On a conda-forge RDKit it
ships under $CONDA_PREFIX/share/RDKit/Contrib. If it is ever missing, vendor
sascorer.py + fpscores.pkl.gz into the repo rather than substituting a different
synthetic-accessibility metric — all four approaches report this axis and a
substitution would break comparability silently.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors

from . import smiles as smi

log = logging.getLogger(__name__)

# The exact column set every approach reports. Order is stable so parquet
# schemas match across approaches.
DESCRIPTOR_COLUMNS: tuple[str, ...] = (
    "MW",
    "HAC",
    "cLogP",
    "TPSA",
    "HBD",
    "HBA",
    "rot_bonds",
    "aromatic_rings",
    "aliphatic_rings",
    "total_rings",
    "frac_sp3",
    "formal_charge",
    "n_stereocenters",
    "molar_refractivity",
    "QED",
    "SAscore",
)

_sascorer = None
_sascore_unavailable = False


def _load_sascorer():
    """Import RDKit's Contrib sascorer, or a vendored copy.

    Returns
    -------
    module or None
        The sascorer module, or None if unavailable (logged once).
    """
    global _sascorer, _sascore_unavailable
    if _sascorer is not None or _sascore_unavailable:
        return _sascorer

    candidates: list[Path] = []
    try:
        from rdkit.Chem import RDConfig

        candidates.append(Path(RDConfig.RDContribDir) / "SA_Score")
    except Exception:  # noqa: BLE001 - RDConfig missing is itself the signal
        pass
    # Vendored fallback, per the plan's install note.
    candidates.append(Path(__file__).resolve().parent / "_vendor" / "SA_Score")

    for d in candidates:
        if (d / "sascorer.py").is_file():
            sys.path.append(str(d))
            try:
                import sascorer  # type: ignore

                _sascorer = sascorer
                log.info("SAscore loaded from %s", d)
                return _sascorer
            except Exception as exc:  # noqa: BLE001
                log.warning("SAscore present at %s but failed to import: %s", d, exc)

    _sascore_unavailable = True
    log.error(
        "SAscore unavailable — checked %s. All four approaches report this axis; "
        "vendor sascorer.py + fpscores.pkl.gz rather than substituting a "
        "different metric.",
        [str(c) for c in candidates],
    )
    return None


def sascore(mol: Chem.Mol) -> float | None:
    """Synthetic accessibility score (Ertl & Schuffenhauer 2009), 1=easy..10=hard."""
    scorer = _load_sascorer()
    if scorer is None:
        return None
    try:
        return float(scorer.calculateScore(mol))
    except Exception as exc:  # noqa: BLE001
        log.debug("SAscore failed: %s", exc)
        return None


def compute(smiles: str) -> dict[str, float | None]:
    """Compute the canonical descriptor set for one molecule.

    Parameters
    ----------
    smiles : str
        Candidate SMILES (canonical or not).

    Returns
    -------
    dict
        Keys are exactly ``DESCRIPTOR_COLUMNS``. All values are None when the
        molecule cannot be parsed — never partially filled, so a downstream
        coverage-aware scorer can tell "not computed" from "computed as zero".

    Examples
    --------
    >>> d = compute("CCO")
    >>> round(d["MW"], 2)
    46.07
    """
    mol = smi.to_mol(smiles)
    if mol is None:
        return {c: None for c in DESCRIPTOR_COLUMNS}

    return {
        "MW": Descriptors.MolWt(mol),
        "HAC": mol.GetNumHeavyAtoms(),
        "cLogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "rot_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "total_rings": rdMolDescriptors.CalcNumRings(mol),
        "frac_sp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "n_stereocenters": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
        "molar_refractivity": Crippen.MolMR(mol),
        "QED": QED.qed(mol),
        "SAscore": sascore(mol),
    }


def compute_frame(df: pd.DataFrame, smiles_col: str = "canonical_smiles") -> pd.DataFrame:
    """Add the descriptor columns to a candidate frame.

    Parameters
    ----------
    df : pandas.DataFrame
        Candidate frame; must carry ``smiles_col``.
    smiles_col : str, optional
        Column holding the SMILES.

    Returns
    -------
    pandas.DataFrame
        A copy of ``df`` with ``DESCRIPTOR_COLUMNS`` added. Row count is
        unchanged — this function never drops candidates, per stamp-don't-delete.
    """
    if smiles_col not in df.columns:
        raise KeyError(f"frame has no {smiles_col!r} column")
    records = [compute(s) for s in df[smiles_col]]
    desc = pd.DataFrame.from_records(records, index=df.index)
    return pd.concat([df.drop(columns=[c for c in DESCRIPTOR_COLUMNS if c in df.columns]),
                      desc], axis=1)


def in_envelope(desc: dict[str, float | None], envelope: dict[str, list]) -> bool:
    """Test a descriptor dict against a developability envelope.

    Deliberately NOT a fixed Lipinski box: the acceptable property window is a
    property of the target and the intended route, not a universal constant, so
    it is read from the target dossier.

    Parameters
    ----------
    desc : dict
        Output of :func:`compute`.
    envelope : dict
        ``{descriptor_name: [min, max]}``. Descriptors absent from the envelope
        are unconstrained.

    Returns
    -------
    bool
        True when every constrained descriptor is inside its window. A None
        value (uncomputable) fails closed.
    """
    for name, bounds in envelope.items():
        val = desc.get(name)
        if val is None:
            return False
        lo, hi = bounds
        if lo is not None and val < lo:
            return False
        if hi is not None and val > hi:
            return False
    return True
