"""
Purpose: Structural-alert filtering — PAINS / BRENK / NIH, with R-group scoping.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: candidate SMILES; optionally a core SMARTS to scope alerts away from
Output: per-catalog alert counts and names, plus a two-tier verdict

THE TWO-TIER PROBLEM, WHICH IS THE WHOLE REASON THIS MODULE EXISTS.

A covalent warhead IS a reactive electrophile. That is the point of it. Run a
whole-molecule PAINS/BRENK filter over a covalent candidate and it will flag the
warhead as a liability — correctly, by the catalog's own rules, and uselessly,
because the warhead is the mechanism.

So for the covalent approaches alerts are computed twice:

  whole molecule  -> ADVISORY. Recorded, never gating. It will flag warheads.
  R-group alone   -> HARD GATE. The R-group is decoration; a PAINS motif there
                     is a genuine liability with no mechanistic excuse.

Isolating the R-group needs care: cutting a substituent off a scaffold leaves an
open valence, and an unsatisfied radical matches substructure patterns that the
real molecule never presents. The fragment is therefore capped with a neutral
methyl before it is scored (see `isolate_rgroup`).

For the NON-covalent approaches (T_1, T_2) there is no warhead to excuse, so the
whole-molecule result is the one that matters — and removing reactive
electrophiles is exactly the discriminator that makes those approaches
non-covalent.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from functools import lru_cache

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from . import smiles as smi

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

CATALOGS = ("PAINS", "BRENK", "NIH")


class AlertError(RuntimeError):
    """Alert filtering could not be applied."""


@lru_cache(maxsize=None)
def _catalog(name: str) -> FilterCatalog.FilterCatalog:
    """Build (once) an RDKit FilterCatalog by name."""
    params = FilterCatalogParams()
    try:
        cat_enum = getattr(FilterCatalogParams.FilterCatalogs, name)
    except AttributeError as exc:
        raise AlertError(f"unknown filter catalog {name!r}") from exc
    params.AddCatalog(cat_enum)
    return FilterCatalog.FilterCatalog(params)


@dataclass
class AlertResult:
    """Alert counts and names for one molecule."""

    counts: dict[str, int] = field(default_factory=dict)
    names: dict[str, list[str]] = field(default_factory=dict)

    @property
    def total(self) -> int:
        return sum(self.counts.values())

    def to_columns(self, prefix: str = "") -> dict:
        out = {f"{prefix}alert_{c.lower()}": self.counts.get(c, 0) for c in CATALOGS}
        out[f"{prefix}alert_total"] = self.total
        out[f"{prefix}alert_names"] = "|".join(
            n for c in CATALOGS for n in self.names.get(c, []))
        return out


def screen(smiles: str, catalogs: tuple[str, ...] = CATALOGS) -> AlertResult:
    """Count structural alerts on a whole molecule.

    Returns
    -------
    AlertResult
        Empty counts when the SMILES cannot be parsed — callers must check
        validity separately rather than reading "0 alerts" as "clean".
    """
    mol = smi.to_mol(smiles)
    res = AlertResult()
    if mol is None:
        return res
    for c in catalogs:
        matches = _catalog(c).GetMatches(mol)
        res.counts[c] = len(matches)
        res.names[c] = [m.GetDescription() for m in matches]
    return res


@lru_cache(maxsize=1)
def _warhead_patterns() -> tuple[tuple[str, Chem.Mol], ...]:
    """Reactive-atom SMARTS per warhead class, for excluding warhead fragments."""
    from . import warhead_library as wl

    out = []
    for _, r in wl.load().iterrows():
        p = Chem.MolFromSmarts(str(r["reactive_atom_smarts"]))
        if p is not None:
            out.append((str(r["class_id"]), p))
    return tuple(out)


def is_warhead_fragment(frag_smiles: str) -> bool:
    """True if a fragment carries any known warhead reactive-atom motif."""
    m = smi.to_mol(frag_smiles)
    if m is None:
        return False
    return any(m.HasSubstructMatch(p) for _, p in _warhead_patterns())


def isolate_rgroup(smiles: str, core_smarts: str, *, cap: str = "C",
                   drop_warhead_fragments: bool = True) -> str | None:
    """Return the R-group(s) with the core removed and open valences capped.

    Parameters
    ----------
    smiles : str
        The whole candidate.
    core_smarts : str
        SMARTS for the protected core to remove. **Do NOT include `[*]`
        attachment wildcards.** A wildcard MATCHES an atom, so `N([*])([*])C1...`
        consumes the first atom of each substituent and the R-group comes back
        truncated — neopentyl returns as four carbons instead of five. Write the
        core as the atoms it actually is: `N[CH]1CCS(=O)(=O)C1`.
    cap : str
        Atom used to satisfy the bond left by the excision. A neutral methyl by
        default: leaving the valence open produces a radical or a bare dummy,
        and both match substructure patterns the real molecule never presents —
        which would generate alerts that are artifacts of the cutting.
    drop_warhead_fragments : bool
        Discard excised fragments that carry a warhead motif. **Required for the
        two-tier gate to mean anything.** The T_4 core has TWO attachment points
        (warhead and R-group), so `ReplaceCore` returns both substituents. Left
        in, the warhead is scored as if it were decoration and the gate flags
        `alkyl_halide` on every chloroacetamide — reintroducing exactly the
        false positive this module exists to avoid.

    Returns
    -------
    str or None
        Dot-separated SMILES of the capped R-group fragment(s), or None if the
        core is absent (itself a finding: the expansion did not preserve it).
    """
    mol = smi.to_mol(smiles)
    core = Chem.MolFromSmarts(core_smarts)
    if mol is None or core is None:
        return None
    core_match = mol.GetSubstructMatch(core)
    if not core_match:
        return None

    # Warhead atoms are identified on the INTACT molecule, before any cutting.
    # Matching fragments after excision does not work: the T_4 core includes the
    # amide nitrogen, so removing it SPLITS an acrylamide into C=O and C=C, and
    # neither piece still matches the Michael-acceptor pattern. The warhead then
    # reappears in the "R-group" and gets flagged — the exact false positive the
    # two-tier design exists to prevent.
    warhead_atoms: set[int] = set()
    if drop_warhead_fragments:
        for _cls, patt in _warhead_patterns():
            for m in mol.GetSubstructMatches(patt):
                warhead_atoms.update(m)

    core_atoms = set(core_match)
    rw = Chem.RWMol(mol)
    for idx in sorted(core_atoms, reverse=True):
        rw.RemoveAtom(idx)
    # Surviving atoms keep their original indices in this map, so fragments can
    # still be tested against warhead_atoms after the core is gone.
    remaining = [i for i in range(mol.GetNumAtoms()) if i not in core_atoms]
    try:
        frag_indices = Chem.GetMolFrags(rw.GetMol())
    except Exception:  # noqa: BLE001 - a malformed excision is a None result
        return None

    kept: list[str] = []
    for frag in frag_indices:
        original = {remaining[i] for i in frag if i < len(remaining)}
        if drop_warhead_fragments and original & warhead_atoms:
            continue                     # this fragment IS (part of) the warhead
        # MolFragmentToSmiles keeps EVERY atom in the fragment. Building from
        # bond paths instead silently drops terminal atoms — neopentyl came back
        # as isobutane — and a gate scoring a fragment that is missing atoms can
        # miss the alert it exists to catch.
        s = Chem.MolFragmentToSmiles(rw.GetMol(), atomsToUse=list(frag),
                                     canonical=True)
        s = s.replace("[*]", cap)
        if smi.to_mol(s) is not None:
            kept.append(s)
    if not kept:
        return None
    joined = ".".join(kept)
    return joined if smi.to_mol(joined) is not None else None


@dataclass
class TwoTierResult:
    """Whole-molecule (advisory) plus R-group (gating) alert results."""

    whole: AlertResult
    rgroup: AlertResult | None
    rgroup_smiles: str | None
    passes_gate: bool
    reason: str = ""

    def to_columns(self) -> dict:
        out = self.whole.to_columns(prefix="whole_")
        if self.rgroup is not None:
            out.update(self.rgroup.to_columns(prefix="rgroup_"))
        out["rgroup_smiles"] = self.rgroup_smiles or ""
        out["alert_gate_pass"] = self.passes_gate
        out["alert_gate_reason"] = self.reason
        return out


def two_tier(smiles: str, core_smarts: str, *,
             max_rgroup_alerts: int = 0) -> TwoTierResult:
    """Advisory whole-molecule alerts + a gating R-group check (covalent approaches).

    Parameters
    ----------
    max_rgroup_alerts : int
        Alerts tolerated on the isolated R-group. Zero by default: a PAINS or
        BRENK motif on the decoration has no mechanistic justification, unlike
        one on the warhead.
    """
    whole = screen(smiles)
    rg_smiles = isolate_rgroup(smiles, core_smarts)
    if rg_smiles is None:
        # Core absent, or the excision failed. Not a pass: verify-after-expansion
        # exists precisely to catch products that lost the core.
        return TwoTierResult(whole=whole, rgroup=None, rgroup_smiles=None,
                             passes_gate=False,
                             reason="core not found or R-group could not be isolated")
    rg = screen(rg_smiles)
    ok = rg.total <= max_rgroup_alerts
    reason = ("" if ok else
              f"R-group carries {rg.total} alert(s): "
              + ", ".join(n for c in CATALOGS for n in rg.names.get(c, []))[:180])
    return TwoTierResult(whole=whole, rgroup=rg, rgroup_smiles=rg_smiles,
                         passes_gate=ok, reason=reason)


def screen_frame(df: pd.DataFrame, smiles_col: str = "canonical_smiles", *,
                 core_smarts: str | None = None,
                 disqualifying: tuple[str, ...] = ()) -> pd.DataFrame:
    """Add alert columns to a candidate frame, stamping rather than deleting.

    Parameters
    ----------
    core_smarts : str, optional
        When given, the two-tier covalent treatment is applied. When omitted,
        whole-molecule alerts are computed and only ``disqualifying`` catalogs
        can stamp a rejection.
    disqualifying : tuple of str
        Catalog names whose hits reject (non-covalent path). Everything else is
        kept as a weighable label — the spec is explicit that alert counts are
        for the panel to weigh, not silent filters.

    Returns
    -------
    pandas.DataFrame
        A copy with alert columns added and ``rejected_at`` stamped where a gate
        failed. Row count is unchanged: stamp, do not delete.
    """
    if smiles_col not in df.columns:
        raise KeyError(f"frame has no {smiles_col!r} column")
    out = df.copy()
    records = []
    for s in out[smiles_col]:
        if core_smarts:
            records.append(two_tier(s, core_smarts).to_columns())
        else:
            r = screen(s)
            cols = r.to_columns(prefix="whole_")
            bad = sum(r.counts.get(c, 0) for c in disqualifying)
            cols["alert_gate_pass"] = bad == 0
            cols["alert_gate_reason"] = (
                "" if bad == 0 else
                f"{bad} hit(s) in disqualifying catalog(s) {list(disqualifying)}")
            records.append(cols)
    alerts = pd.DataFrame.from_records(records, index=out.index)
    out = pd.concat([out.drop(columns=[c for c in alerts.columns if c in out.columns]),
                     alerts], axis=1)

    if "rejected_at" not in out.columns:
        out["rejected_at"] = pd.NA
    stamp = out["alert_gate_pass"].eq(False) & out["rejected_at"].isna()
    out.loc[stamp, "rejected_at"] = "alerts"
    log.info("alerts: %d/%d stamped rejected (row set unchanged)",
             int(stamp.sum()), len(out))
    return out
