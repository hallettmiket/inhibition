"""
Purpose: Data-driven warhead-class library for T_4 — narrow now, wide later.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: data/reference/warhead_classes_10.csv
Output: validated warhead classes, filtered by how well-founded each one is

DESIGN INTENT. The set of warhead chemistries is DATA, not code. Going wide —
enumerating exploratory electrophiles beyond the Pin1-precedented ones — is
adding rows to a CSV and re-running, never editing a pipeline. That is the same
pattern seeds.yaml uses for seed molecules, and it is what keeps the
choreography retargetable.

But breadth without provenance is how a library ends up half-dead. In the prior
run 6 of 16 warhead classes collapsed to an inert amide once attached to the
core, invisibly. So every class carries a `structure_status` recording how far
its structure is actually established, and the loader makes you say out loud
which tiers you are willing to enumerate:

  VERIFIED             structure traced to a public record or derived from a
                       verified anchor. Safe to enumerate.
  VERIFIED_CLASS_ONLY  the warhead CHEMOTYPE is verified, but the specific
                       literature compound is not.
  DESIGNED_UNTESTED    chemotype verified AND an attachment regiochemistry has
                       been proposed, but which attachment is right has not been
                       established. Competing regiochemistries are enumerated as
                       SEPARATE classes so the 5b validity gate and the LUMO
                       window decide empirically rather than by intuition.
  NEEDS_DESIGN         the anchor exists but gives no attachable fragment;
                       somebody has to design the attachment. Not enumerable.
  UNVERIFIED           structure unresolved. Never enters the reactivity window.

`enumerable()` defaults to VERIFIED only. Widening is deliberate and explicit.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from . import smiles as smi

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LIBRARY = _REPO_ROOT / "data" / "reference" / "warhead_classes_10.csv"

UNRESOLVED = "UNRESOLVED"

# ---------------------------------------------------------------------------
# Free-text warhead descriptions -> the canonical class_id vocabulary
# ---------------------------------------------------------------------------
# The reference set's `warhead_class` column is PROSE, written per row by
# whoever curated the entry. The warhead library's `class_id` is a controlled
# vocabulary. Counting `nunique()` on the prose is wrong and it errs in the
# dangerous direction -- it OVERCOUNTS:
#
#   "chloroacetamide"                            \  the same warhead,
#   "chloroacetamide (N-methyl peptidomimetic)"  /  counted as two classes
#
# On the current lead-tier actives the prose gives 6 classes and the canonical
# mapping gives 4. Six is the gate's floor, so the sloppy count clears it and
# the correct one does not (D0045). A string-uniqueness count would have
# certified a verdict the evidence does not support.
#
# EXPLICIT, AND LOUD ON ANYTHING NEW. This is a hand-checked mapping of the
# exact strings present, not a fuzzy matcher. A warhead description that is not
# in here raises rather than silently becoming its own class -- which is the
# same failure in a new costume.
_PROSE_TO_CLASS_ID: dict[str, str] = {
    "chloroacetamide": "chloroacetamide",
    "chloroacetamide (N-methyl peptidomimetic)": "chloroacetamide",
    "sulfamate acetamide": "sulfamate_acetamide",
    "cinnamamide (aryl Michael acceptor; acrylamide-class)": "cinnamamide",
    "2-chloro-5-nitropyrimidine (SNAr electrophile, X-ray-confirmed "
    "Cys113 adduct)": "snar_chloroazine",
    "1;4-naphthoquinone (Michael acceptor)": "naphthoquinone_c2",
    # KPT-6566 is self-immolative: it RELEASES a naphthoquinone, so its
    # reactive species is that class. It is promiscuous and dual-mechanism, so
    # it is tiered out of `lead` anyway -- mapped for completeness, not use.
    "aryl-sulfonyl-acetate (self-immolative -> 1;4-naphthoquinone)":
        "naphthoquinone_c2",
    "3-bromo-4;5-dihydroisoxazole (BDHI)": "bdhi_c4",
}


class WarheadClassError(RuntimeError):
    """A warhead description has no canonical class."""


def canonical_class(prose: str | None) -> str | None:
    """Canonical ``class_id`` for a reference set's free-text warhead.

    Returns None for a genuinely absent or unverified warhead -- a
    non-covalent binder has none, and an entry whose warhead has not been
    established should not be assigned one. Callers must decide what an
    unknown means for them; the gate excludes it from the chemotype count
    rather than inventing a class for it.

    Raises
    ------
    WarheadClassError
        If the string is present, non-trivial, and not in the mapping. Adding
        an active with a new warhead description is exactly when someone must
        decide which class it belongs to, so this fails at that moment instead
        of quietly inflating the denominator later.
    """
    if prose is None:
        return None
    s = str(prose).strip()
    if not s or s.lower() in {"nan", "none", UNRESOLVED.lower(), "unverified"}:
        return None
    if s in _PROSE_TO_CLASS_ID:
        return _PROSE_TO_CLASS_ID[s]
    raise WarheadClassError(
        f"warhead description {s!r} has no canonical class_id. Add it to "
        "shared.warhead_library._PROSE_TO_CLASS_ID after deciding which class "
        "it belongs to — do not let it become its own class by default, which "
        "would inflate the chemotype count the gate's floor depends on.")

# Ordered from best-founded to least. Membership in this tuple is the whole
# vocabulary — an unknown status is an error, not a silently-tolerated value.
STATUS_TIERS: tuple[str, ...] = (
    "VERIFIED",
    "VERIFIED_CLASS_ONLY",
    "DESIGNED_UNTESTED",
    "NEEDS_DESIGN",
    "UNVERIFIED",
)

REQUIRED_COLUMNS = {
    "class_id",
    "display_name",
    "warhead_fragment_smiles",
    "mechanism",
    "reactive_atom_smarts",
    "precedent",
    "structure_status",
    "structure_source",
    # v4 (D0022): what gnina must actually be given. `reactive_atom_smarts`
    # describes the UNREACTED warhead and names its leaving group, so it cannot
    # match the adduct; docking needs a pattern valid on the product.
    "adduct_attachment_smarts",
    "has_leaving_group",
}


class WarheadLibraryError(RuntimeError):
    """The warhead library is malformed, or is being used past its evidence."""


def load(path: Path | str | None = None) -> pd.DataFrame:
    """Load and validate the warhead-class library.

    Every class with a resolved fragment SMILES must parse and must carry an
    attachment point, because a "fragment" with nowhere to attach cannot be
    coupled to the core and would fail silently at enumeration time.

    Raises
    ------
    WarheadLibraryError
        On a missing file, missing column, unknown status, unparseable
        fragment, or a fragment with no attachment point.
    """
    p = Path(path) if path else DEFAULT_LIBRARY
    if not p.is_file():
        raise WarheadLibraryError(f"warhead library not found: {p}")

    df = pd.read_csv(p)
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise WarheadLibraryError(f"{p.name} missing column(s): {sorted(missing)}")

    bad_status = set(df["structure_status"]) - set(STATUS_TIERS)
    if bad_status:
        raise WarheadLibraryError(
            f"unknown structure_status {sorted(bad_status)}; allowed: {list(STATUS_TIERS)}"
        )

    if df["class_id"].duplicated().any():
        dupes = df.loc[df["class_id"].duplicated(), "class_id"].tolist()
        raise WarheadLibraryError(f"duplicate class_id(s): {dupes}")

    for _, row in df.iterrows():
        frag = str(row["warhead_fragment_smiles"])
        if frag == UNRESOLVED:
            continue
        if smi.to_mol(frag, sanitize=False) is None:
            raise WarheadLibraryError(
                f"class {row['class_id']!r}: unparseable fragment SMILES {frag!r}"
            )
        n_attach = frag.count("*")
        if n_attach == 0:
            raise WarheadLibraryError(
                f"class {row['class_id']!r}: fragment {frag!r} has no attachment "
                "point '[*]' — it cannot be coupled to the core."
            )
        if n_attach > 1:
            # A second attachment point is never filled by the core coupling,
            # so the product carries a dangling dummy atom: not a real molecule,
            # yet it passes core verification and reaches docking. This happened
            # to 198 sulfamate_acetamide candidates.
            raise WarheadLibraryError(
                f"class {row['class_id']!r}: fragment {frag!r} has {n_attach} "
                "attachment points. Exactly one is required — the coupling fills "
                "only the core bond, so any other '[*]' survives into the product "
                "as a dangling dummy atom. Fix the substituent in the fragment."
            )
    return df


def enumerable(
    df: pd.DataFrame | None = None,
    *,
    allow_statuses: tuple[str, ...] = ("VERIFIED",),
) -> pd.DataFrame:
    """Warhead classes cleared for combinatorial enumeration.

    Parameters
    ----------
    df : pandas.DataFrame, optional
        A preloaded library; loaded from disk when omitted.
    allow_statuses : tuple of str, optional
        Which evidence tiers to enumerate. Defaults to VERIFIED only —
        deliberately narrow. Widening is a caller decision that gets logged.

    Returns
    -------
    pandas.DataFrame
        Classes whose status is allowed AND whose fragment SMILES is resolved.

    Notes
    -----
    A class can be well-founded as chemistry and still be unenumerable, because
    its attachment to the core has not been designed. Those are reported, not
    quietly dropped.
    """
    lib = df if df is not None else load()
    unknown = set(allow_statuses) - set(STATUS_TIERS)
    if unknown:
        raise WarheadLibraryError(f"unknown status filter: {sorted(unknown)}")

    allowed = lib[lib["structure_status"].isin(allow_statuses)]
    resolved = allowed[allowed["warhead_fragment_smiles"] != UNRESOLVED]

    blocked = lib[~lib["class_id"].isin(resolved["class_id"])]
    for _, row in blocked.iterrows():
        log.warning(
            "warhead class %r excluded from enumeration (status=%s, fragment=%s): %s",
            row["class_id"],
            row["structure_status"],
            "resolved" if row["warhead_fragment_smiles"] != UNRESOLVED else UNRESOLVED,
            str(row.get("notes", ""))[:120],
        )

    if allow_statuses != ("VERIFIED",):
        log.warning(
            "enumerating beyond VERIFIED: %s — every widened class needs a "
            "chemist's sign-off on its attachment chemistry before it is trusted.",
            list(allow_statuses),
        )
    return resolved


def window_anchor_classes(df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Classes admissible as reactivity-window anchors.

    Control B5: the window is bounded by real, wet-lab-validated covalent-Cys113
    actives. A class whose structure is unresolved cannot anchor a window that
    is supposed to be empirical, so UNVERIFIED is refused here exactly as
    reference_set.covalent_anchors() refuses it.
    """
    lib = df if df is not None else load()
    anchors = lib[
        (lib["precedent"] == "anchored_verified")
        & (lib["structure_status"].isin(("VERIFIED", "VERIFIED_CLASS_ONLY")))
        & (lib["warhead_fragment_smiles"] != UNRESOLVED)
    ]
    if len(anchors) < 2:
        raise WarheadLibraryError(
            f"only {len(anchors)} usable anchor class(es); a reactivity window "
            "anchored on one chemotype is the n~1 problem control B5 exists to "
            "prevent. Resolve more anchor structures before computing a window."
        )
    return anchors


def status_report(df: pd.DataFrame | None = None) -> str:
    """Human-readable summary of what is and is not usable, and why."""
    lib = df if df is not None else load()
    lines = [f"Warhead library: {len(lib)} class(es) declared"]
    for tier in STATUS_TIERS:
        rows = lib[lib["structure_status"] == tier]
        if rows.empty:
            continue
        ids = ", ".join(rows["class_id"])
        lines.append(f"  {tier:20s} {len(rows)}  [{ids}]")
    enum = enumerable(lib)
    lines.append(f"Enumerable now (VERIFIED only): {len(enum)} -> {list(enum['class_id'])}")
    return "\n".join(lines)
