"""
Purpose: Canonical SMILES + InChIKey handling — the choreography's join key.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: SMILES strings from any approach (generated, enumerated, or reference)
Output: canonical SMILES, InChIKeys, deterministic candidate ids

Every approach keys its candidates on canonical SMILES, and the integration
phase pools the four shortlists on that key. If two approaches canonicalize
differently, the same molecule appears twice in the pool and the structural
convergence signal — the most defensible cross-approach evidence there is —
silently degrades. So canonicalization happens here and only here.
"""

from __future__ import annotations

import hashlib
import logging
from typing import Iterable

from rdkit import Chem
from rdkit import RDLogger

log = logging.getLogger(__name__)

# RDKit logs every sanitization failure to stderr. At 10^5 candidates that
# buries real output, and we record failures ourselves in the frame.
RDLogger.DisableLog("rdApp.*")


class SmilesError(ValueError):
    """A SMILES string could not be parsed or canonicalized."""


def to_mol(smiles: str, *, sanitize: bool = True) -> Chem.Mol | None:
    """Parse a SMILES string into an RDKit molecule.

    Parameters
    ----------
    smiles : str
        The SMILES string to parse.
    sanitize : bool, optional
        Run RDKit sanitization (valence, aromaticity, kekulization).

    Returns
    -------
    Chem.Mol or None
        The molecule, or None if it could not be parsed.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    return Chem.MolFromSmiles(smiles, sanitize=sanitize)


def canonical(smiles: str, *, strict: bool = False) -> str | None:
    """Canonicalize a SMILES string.

    Parameters
    ----------
    smiles : str
        Input SMILES, in any valid form.
    strict : bool, optional
        Raise on failure instead of returning None.

    Returns
    -------
    str or None
        Canonical isomeric SMILES, or None when unparseable and not strict.

    Raises
    ------
    SmilesError
        If strict and the molecule cannot be parsed.

    Examples
    --------
    >>> canonical("C(C)O") == canonical("CCO")
    True
    """
    mol = to_mol(smiles)
    if mol is None:
        if strict:
            raise SmilesError(f"cannot parse SMILES: {smiles!r}")
        return None
    # isomericSmiles=True keeps stereochemistry: ATRA is all-trans and the
    # sulfopin core carries a defined stereocentre, so dropping it would merge
    # genuinely distinct candidates.
    return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)


def inchikey(smiles: str) -> str | None:
    """Return the InChIKey for a SMILES string, or None if unparseable.

    InChIKey is the dedup key for CReM frontiers in T_2. Path multiplicity is
    high — the same molecule is reachable by many edit sequences — so an
    undeduped frontier explodes combinatorially on the next breadth-first step.
    """
    mol = to_mol(smiles)
    if mol is None:
        return None
    try:
        key = Chem.MolToInchiKey(mol)
    except Exception as exc:  # RDKit raises bare Exception on InChI failures
        log.debug("InChIKey failed for %r: %s", smiles, exc)
        return None
    # RDKit returns an EMPTY STRING (not None) for molecules it cannot key —
    # notably anything carrying a dummy atom. Treating "" as success let 198
    # malformed candidates share one id, because sha256("") is a perfectly
    # valid-looking hash. Empty is failure.
    return key or None


def candidate_id(smiles: str, *, prefix: str = "") -> str | None:
    """Build a deterministic candidate id from a molecule's InChIKey.

    Deterministic rather than sequential so that a resumed or re-run stage
    produces the same ids, and so ids computed in different approaches for the
    same molecule collide on purpose.

    Parameters
    ----------
    smiles : str
        Candidate SMILES.
    prefix : str, optional
        Short approach tag (e.g. ``"t4"``) prepended to the hash.

    Returns
    -------
    str or None
        ``<prefix>_<12-hex-chars>``, or None if the SMILES is unparseable.
    """
    key = inchikey(smiles)
    if not key:                      # None OR "" — see inchikey() above
        return None
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:12]
    return f"{prefix}_{digest}" if prefix else digest


def canonical_many(
    smiles_list: Iterable[str], *, keep_failures: bool = True
) -> list[str | None]:
    """Canonicalize a batch of SMILES.

    Parameters
    ----------
    smiles_list : iterable of str
        Input SMILES.
    keep_failures : bool, optional
        When True, failures are kept in position as None so the result aligns
        row-for-row with the input frame. When False they are dropped.

    Returns
    -------
    list of (str or None)
        Canonical SMILES, positionally aligned when ``keep_failures``.
    """
    out: list[str | None] = []
    for s in smiles_list:
        c = canonical(s)
        if c is None and not keep_failures:
            continue
        out.append(c)
    return out


def same_molecule(a: str, b: str) -> bool:
    """True when two SMILES denote the same molecule after canonicalization."""
    ca, cb = canonical(a), canonical(b)
    return ca is not None and ca == cb


def has_substructure(smiles: str, pattern_smarts: str) -> bool:
    """Test whether a molecule contains a SMARTS substructure.

    Used for verify-after-expansion: every product of a constrained expansion
    (T_3's core+warhead, T_4's core) must still contain the intended fixed
    substructure. A reaction definition correct in general can still misfire on
    an unusual substrate, so this is checked on every row rather than a sample.
    """
    mol = to_mol(smiles)
    patt = Chem.MolFromSmarts(pattern_smarts)
    if mol is None or patt is None:
        return False
    return mol.HasSubstructMatch(patt)
