"""
Purpose: Loader + validator for the frozen Pin1 reference binder set.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: data/reference/pin1_reference_binders_2.csv and
       data/reference/pin1_covalent_cys113_anchors_2.csv
Output: validated frames for the novelty axis and T_4's reactivity window

This module exists to enforce two adversary controls that are easy to violate
by accident:

B4 — the novelty axis must be computed against this EXTERNAL published set,
     never against the seed. Novelty-vs-seed measures edit distance from the
     starting molecule, which mechanically rewards de-novo generation (T_1) and
     penalizes derivative search (T_2) for doing exactly what it was asked to
     do. That is circular, not informative.

B5 — T_4's reactivity window must be anchored on real, wet-lab-validated
     covalent-Cys113 actives, and the project's own computational leads are
     excluded as anchors. A window anchored on your own prediction is not a
     control.

Accordingly `covalent_anchors()` REFUSES to return UNVERIFIED rows by default.
One anchor (Byun 2023 BDHI fragment) still has an SI-only structure and is
gated out until its SMILES is resolved. Reddi 2023 was resolved from Figure 5.
"""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from . import smiles as smi

log = logging.getLogger(__name__)

UNVERIFIED = "UNVERIFIED"

# Repo root, resolved from this file's location so callers need no cwd contract.
_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MASTER = _REPO_ROOT / "data" / "reference" / "pin1_reference_binders_3.csv"
DEFAULT_ANCHORS = _REPO_ROOT / "data" / "reference" / "pin1_covalent_cys113_anchors_2.csv"


class ReferenceSetError(RuntimeError):
    """The reference set is missing, malformed, or used in a disallowed way."""


@dataclass(frozen=True)
class ReferenceSet:
    """The frozen reference set plus the hashes that pin it.

    Attributes
    ----------
    master : pandas.DataFrame
        All validated Pin1 binders (novelty axis, every approach).
    anchors : pandas.DataFrame
        Covalent-at-Cys113 subset (T_4 reactivity window).
    master_sha256 : str
        Hash of the master CSV as read, recorded into every manifest.
    anchors_sha256 : str
        Hash of the anchors CSV as read.
    """

    master: pd.DataFrame
    anchors: pd.DataFrame
    master_sha256: str
    anchors_sha256: str


def _sha256(path: Path) -> str:
    """Return the SHA-256 of a file, for manifest pinning."""
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def load(
    master_path: Path | str | None = None,
    anchors_path: Path | str | None = None,
) -> ReferenceSet:
    """Load and validate both reference CSVs.

    Every non-UNVERIFIED SMILES is required to parse and canonicalize. A
    reference set with a broken structure in it silently corrupts the novelty
    axis for every approach, so this fails loudly at load rather than producing
    quiet nonsense downstream.

    Raises
    ------
    ReferenceSetError
        If a file is missing, a required column is absent, or a non-UNVERIFIED
        SMILES fails to parse.
    """
    mpath = Path(master_path) if master_path else DEFAULT_MASTER
    apath = Path(anchors_path) if anchors_path else DEFAULT_ANCHORS

    for p in (mpath, apath):
        if not p.is_file():
            raise ReferenceSetError(f"reference file not found: {p}")

    master = pd.read_csv(mpath)
    anchors = pd.read_csv(apath)

    _require_columns(master, {"name", "canonical_smiles", "mechanism"}, mpath)
    _require_columns(
        anchors,
        {"anchor_id", "name", "warhead_class", "canonical_smiles", "smiles_status"},
        apath,
    )

    _validate_smiles_column(master, "canonical_smiles", mpath)
    _validate_smiles_column(anchors, "canonical_smiles", apath)

    return ReferenceSet(
        master=master,
        anchors=anchors,
        master_sha256=_sha256(mpath),
        anchors_sha256=_sha256(apath),
    )


def _require_columns(df: pd.DataFrame, needed: set[str], path: Path) -> None:
    """Raise if a required column is missing."""
    missing = needed - set(df.columns)
    if missing:
        raise ReferenceSetError(f"{path.name} missing column(s): {sorted(missing)}")


def _validate_smiles_column(df: pd.DataFrame, col: str, path: Path) -> None:
    """Require every non-UNVERIFIED SMILES in a column to canonicalize."""
    bad: list[str] = []
    for _, row in df.iterrows():
        val = str(row[col])
        if val == UNVERIFIED:
            continue
        if smi.canonical(val) is None:
            bad.append(f"{row.get('name', '?')}: {val!r}")
    if bad:
        raise ReferenceSetError(
            f"{path.name} has {len(bad)} unparseable SMILES: " + "; ".join(bad[:5])
        )


def master_set(refset: ReferenceSet | None = None) -> list[str]:
    """Canonical SMILES of every verified binder — the novelty reference.

    Returns
    -------
    list of str
        Canonical SMILES, UNVERIFIED rows dropped. This list, and nothing else,
        is what `novelty.py` computes Tanimoto against.
    """
    rs = refset or load()
    out: list[str] = []
    for val in rs.master["canonical_smiles"]:
        val = str(val)
        if val == UNVERIFIED:
            continue
        c = smi.canonical(val)
        if c is not None:
            out.append(c)
    return out


def covalent_anchors(
    refset: ReferenceSet | None = None, *, verified_only: bool = True
) -> pd.DataFrame:
    """Covalent-at-Cys113 anchors for T_4's reactivity window.

    Parameters
    ----------
    refset : ReferenceSet, optional
        A preloaded set; loaded from disk when omitted.
    verified_only : bool, optional
        Keep the default True. Setting it False returns SI-only structures
        that have not been confirmed against a public record, which is exactly
        what control B5 exists to prevent.

    Returns
    -------
    pandas.DataFrame
        Anchor rows with usable SMILES.

    Raises
    ------
    ReferenceSetError
        If fewer than two distinct electrophile chemistries survive — a window
        anchored on one chemotype is the n≈1 problem B5 was raised about.
    """
    rs = refset or load()
    df = rs.anchors
    if verified_only:
        df = df[df["smiles_status"] != UNVERIFIED].copy()
        dropped = len(rs.anchors) - len(df)
        if dropped:
            log.warning(
                "reference_set: %d UNVERIFIED anchor(s) withheld from the "
                "reactivity window (SI-only structures; see .provenance.md)",
                dropped,
            )

    n_chemistries = df["warhead_class"].nunique()
    if n_chemistries < 2:
        raise ReferenceSetError(
            f"reactivity window needs >= 2 distinct electrophile chemistries, "
            f"got {n_chemistries}. Refusing to anchor a window on one chemotype."
        )
    return df


def window_caveat(refset: ReferenceSet | None = None) -> str:
    """The honest caveat about window composition, for display and manifests.

    The clean, selective anchors are chloroacetamide-dominated and the extra
    kinetic anchors are promiscuous quinones. Any consumer of the window should
    carry this text rather than treating the window as chemotype-balanced.
    """
    df = covalent_anchors(refset)
    classes = sorted(df["warhead_class"].unique())
    promiscuous = df[df.get("promiscuity_flag", "n") == "y"]["name"].tolist()
    return (
        f"Reactivity window anchored on {len(df)} verified covalent-Cys113 actives "
        f"spanning {len(classes)} electrophile chemistries: {', '.join(classes)}. "
        f"The clean, selective anchors are chloroacetamide-dominated; "
        f"{len(promiscuous)} anchor(s) are promiscuous quinones "
        f"({', '.join(promiscuous) or 'none'}). The window is therefore "
        "chloroacetamide-centric with a reactive-quinone upper tail, NOT a "
        "chemotype-balanced distribution. Do not over-trust a single global cutoff."
    )
