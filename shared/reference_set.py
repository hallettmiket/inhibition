"""
Purpose: Loader + validator for the frozen Pin1 reference binder set.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: data/reference/pin1_reference_binders_<latest>.csv and
       data/reference/pin1_covalent_cys113_anchors_<latest>.csv
       (resolved by glob — see latest_reference(); naming a version here is
        how this header went stale twice already)
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
_REF_DIR = _REPO_ROOT / "data" / "reference"


def latest_reference(stem: str) -> Path:
    """Highest integer-versioned ``<stem>_N.csv`` under data/reference.

    RESOLVED BY GLOB, NOT PINNED BY HAND. Every reference file in this project
    is integer-versioned per the lab data rule, and every hand-written pin to a
    specific version has eventually gone stale while looking fine:

    * `warhead_library.DEFAULT_LIBRARY` sat on `warhead_classes_7.csv` after
      `_8` was written, so the acrylamide precedent correction in `_8` was
      never loaded by anything (recorded in data/ready_to_delete.md).
    * `gen_docs.py` named `pin1_reference_binders_1.csv` as a literal and
      carried a load error on the rendered status page for weeks.
    * This module pinned `_3` while `_4` -- which adds two potent covalent
      actives the gate had been reporting UNDERPOWERED without -- sat unread
      beside it.

    The failure mode is always the same and always silent: the pin still points
    at a file that exists and parses, so nothing raises. Globbing gets the same
    answer with no dependency and cannot go stale.

    Raises
    ------
    ReferenceSetError
        If no versioned file matches, which IS worth failing on -- a missing
        reference set must not be papered over with a default.
    """
    best, best_n = None, -1
    for p in _REF_DIR.glob(f"{stem}_*.csv"):
        tail = p.stem[len(stem) + 1:]
        if tail.isdigit() and int(tail) > best_n:
            best, best_n = p, int(tail)
    if best is None:
        raise ReferenceSetError(
            f"no versioned {stem}_N.csv under {_REF_DIR}")
    return best


def warhead_library() -> "pd.DataFrame":
    """The highest-numbered `warhead_classes_N.csv`, read.

    ONE RESOLVER, because two call sites had each written their own and both
    wrote the same bug: `sorted(glob("warhead_classes_*.csv"))[-1]` sorts
    LEXICALLY, so the moment the library reached `_10` the string `"..._9.csv"`
    sorted last and `_10` became unreachable. `_10` adds exactly one class --
    `cinnamamide`, the aryl Michael acceptor with the ESI-MS-confirmed Cys113
    adduct -- so both call sites quietly saw a library with that class missing.

    `dock_reference_modes` failed loudly on it (its lookup raises on a miss).
    `shortlist_report` returns None on a miss, so there it was silent: a
    cinnamamide molecule simply had no warhead resolved, with no error.

    The lesson is the one the module already teaches for pins, one step on: a
    glob is only dynamic resolution if it ORDERS the way the versions are
    numbered. `latest_reference` compares the integer.
    """
    import pandas as pd
    return pd.read_csv(latest_reference("warhead_classes"))


def load_warhead_row(class_id: str) -> "pd.Series":
    """One warhead class by id, from the latest library.

    Raises rather than returning a default: the class list is an ALLOWLIST, and
    a caller handed a name nobody registered must not proceed on a fallback.
    That is the rule `canonical_class()` already follows, and catalogue #29 is
    what happens when it is not followed -- an unregistered mechanism scored
    0.0, the worst LEGAL value, and nothing raised.
    """
    d = warhead_library()
    hit = d[d.class_id == class_id]
    if hit.empty:
        raise ReferenceSetError(
            f"no warhead class {class_id!r} in "
            f"{latest_reference('warhead_classes').name}; known: "
            f"{sorted(d.class_id)}")
    return hit.iloc[0]


def _default_master() -> Path:
    return latest_reference("pin1_reference_binders")


def _default_anchors() -> Path:
    return latest_reference("pin1_covalent_cys113_anchors")


# Kept as module attributes because callers and tests reference them by name.
# Resolved at import from the newest file on disk rather than written literally.
DEFAULT_MASTER = _default_master()
DEFAULT_ANCHORS = _default_anchors()


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
    _validate_mechanism(master, mpath)

    return ReferenceSet(
        master=master,
        anchors=anchors,
        master_sha256=_sha256(mpath),
        anchors_sha256=_sha256(apath),
    )


#: The ONLY mechanism values the pipeline recognises. The enrichment gate
#: selects its strata by exact equality (`master.mechanism == "non_covalent"`),
#: so a value outside this set does not raise anywhere -- it simply matches no
#: stratum and the active vanishes from every gate that should have used it.
#: Liu-2024-C3 was written as "non_covalent_active_site" and silently dropped,
#: costing the sixth independent chemotype, which is the gate's verdict floor.
#: The log said "5 actives" and nothing said one had been discarded.
VALID_MECHANISMS = frozenset({"covalent_cys113", "non_covalent", "unknown"})


def _validate_mechanism(df: pd.DataFrame, path: Path) -> None:
    """Raise if any row carries a mechanism the gate cannot select on."""
    bad = sorted(set(df["mechanism"].dropna().astype(str)) - VALID_MECHANISMS)
    if bad:
        names = df[df["mechanism"].astype(str).isin(bad)]["name"].tolist()
        raise ReferenceSetError(
            f"{path.name}: unrecognised mechanism value(s) {bad} on {names}. "
            f"Allowed: {sorted(VALID_MECHANISMS)}. A value outside this set is "
            "not an error anywhere downstream -- it matches no stratum, so the "
            "affected binder is silently missing from every gate.")


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
