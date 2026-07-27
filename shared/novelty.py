"""
Purpose: The novelty axis — 1 - max Tanimoto (ECFP4) to the frozen EXTERNAL set.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: candidate SMILES; the frozen reference set from reference_set.py
Output: a novelty score in [0, 1] per candidate

THE CONTROL THIS ENFORCES (adversary finding B4). Novelty must be measured
against the published, external known-binder set — never against the approach's
own seed. Seed-relative novelty is circular: it measures edit distance from the
starting molecule, which mechanically rewards T_1 (which has no seed, so
everything looks novel) and penalizes T_2 (whose entire job is to stay near
ATRA) for succeeding at its task. The resulting axis would rank approaches, not
molecules.

There is therefore no seed parameter in this module, by design. Do not add one.

ECFP4 = Morgan fingerprint, radius 2 (Rogers & Hahn 2010). Tanimoto on binary
ECFP4 is the standard, well-justified choice for this comparison
(Bajusz, Racz & Heberger 2015).
"""

from __future__ import annotations

import logging
from functools import lru_cache

import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator

from . import reference_set as refset_mod
from . import smiles as smi

log = logging.getLogger(__name__)

ECFP4_RADIUS = 2
ECFP4_NBITS = 2048

_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=ECFP4_RADIUS, fpSize=ECFP4_NBITS
)


def fingerprint(smiles: str):
    """ECFP4 bit vector for a SMILES, or None if unparseable."""
    mol = smi.to_mol(smiles)
    if mol is None:
        return None
    return _generator.GetFingerprint(mol)


@lru_cache(maxsize=1)
def _reference_fingerprints() -> tuple:
    """ECFP4 fingerprints of the frozen external reference set.

    Cached: this is loaded once and reused across ~10^5 candidates. The cache is
    keyed on nothing because the set is frozen for the life of a run — if the
    reference CSVs change, the process must be restarted (and the manifest hash
    will differ, which is the point).
    """
    ref_smiles = refset_mod.master_set()
    fps = []
    for s in ref_smiles:
        fp = fingerprint(s)
        if fp is not None:
            fps.append(fp)
    if not fps:
        raise refset_mod.ReferenceSetError(
            "no usable reference fingerprints — the novelty axis cannot be "
            "computed, and falling back to seed-relative novelty is forbidden "
            "(control B4)."
        )
    log.info("novelty: %d external reference fingerprints loaded", len(fps))
    return tuple(fps)


def max_similarity_to_reference(smiles: str) -> float | None:
    """Maximum ECFP4 Tanimoto between a candidate and the external set.

    Returns
    -------
    float or None
        Max Tanimoto in [0, 1], or None if the candidate is unparseable.
    """
    fp = fingerprint(smiles)
    if fp is None:
        return None
    sims = DataStructs.BulkTanimotoSimilarity(fp, list(_reference_fingerprints()))
    return float(max(sims)) if sims else 0.0


def novelty(smiles: str) -> float | None:
    """Novelty of a candidate against the frozen external reference set.

    Parameters
    ----------
    smiles : str
        Candidate SMILES.

    Returns
    -------
    float or None
        ``1 - max_tanimoto``, in [0, 1]; higher is more novel. None when the
        candidate cannot be parsed.

    Notes
    -----
    This is a weighable axis, not a gate. A candidate close to a known binder is
    not thereby disqualified — for T_2, staying near ATRA is the assignment —
    and a highly novel candidate is not thereby good. The panel weighs it.
    """
    sim = max_similarity_to_reference(smiles)
    return None if sim is None else 1.0 - sim


def novelty_frame(
    df: pd.DataFrame, smiles_col: str = "canonical_smiles", out_col: str = "novelty_external"
) -> pd.DataFrame:
    """Add the novelty column to a candidate frame.

    The column is named ``novelty_external`` rather than ``novelty`` so that any
    seed-relative quantity someone later adds cannot be silently swapped in.
    """
    if smiles_col not in df.columns:
        raise KeyError(f"frame has no {smiles_col!r} column")
    out = df.copy()
    out[out_col] = [novelty(s) for s in out[smiles_col]]
    return out


def pairwise_similarity_matrix(smiles_list: list[str]):
    """Full ECFP4 Tanimoto matrix over a candidate list.

    Used by the integration phase for structural convergence — clustering the
    pooled 40 shortlisted candidates to find molecules (or close analogs)
    surfaced independently by more than one approach. That is the most
    defensible cross-approach signal available, because it relies on no score
    commensurability whatsoever.

    Returns
    -------
    list[list[float]]
        Symmetric matrix of Tanimoto similarities; 0.0 where a SMILES failed.
    """
    fps = [fingerprint(s) for s in smiles_list]
    n = len(fps)
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        if fps[i] is None:
            continue
        for j in range(i, n):
            if fps[j] is None:
                continue
            sim = float(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
            matrix[i][j] = matrix[j][i] = sim
    return matrix
