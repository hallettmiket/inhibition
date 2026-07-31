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
  decoration      -> HARD GATE. Decoration is not mechanism; a PAINS motif there
                     is a genuine liability with no excuse.

HOW THE DECORATION'S SHARE IS MEASURED (D0025). NOT by cutting the molecule up.
The molecule is screened INTACT, and each alert is attributed by asking RDKit
which atoms it matched: inside the excused region (core plus warhead) it is the
mechanism, outside it belongs to the decoration, and a match spanning both is
charged to the decoration but reported separately as a boundary hit.

Excision was tried first and is a trap, because the cut itself creates
chemistry. An amide `>N-C(=O)-R` severed from its nitrogen and capped with
hydrogen is a formamide, which BRENK flags as an aldehyde; a thioether capped
the same way is a thiol. Neither group exists in the real molecule. Attribution
breaks no bonds, so it cannot invent a functional group.

The excused region must cover the WHOLE warhead, not just its reactive atoms.
`warhead_fragment_smiles` defines it. Excusing only the narrow reactive-atom
SMARTS leaves `alpha_halo_carbonyl` straddling the boundary on every
chloroacetamide — the two-tier false positive returning through a new door.

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


@lru_cache(maxsize=1)
def _warhead_group_patterns() -> tuple[tuple[str, Chem.Mol], ...]:
    """Whole-warhead SMARTS per class, from `warhead_fragment_smiles`.

    Distinct from `_warhead_patterns`, which returns the narrow reactive-atom
    SMARTS used to tell gnina where to form the bond. For alert attribution the
    whole group matters: excusing `[CH2][Cl]` but not the carbonyl it is alpha
    to leaves `alpha_halo_carbonyl` spanning the boundary and blamed on the
    decoration.

    The fragment's `[*]` attachment marker is kept: as SMARTS it matches any
    atom, which is what is wanted — it binds to whatever the warhead is attached
    to, and that atom is part of the core in every approach that uses this.
    """
    from . import warhead_library as wl

    out = []
    for _, r in wl.load().iterrows():
        frag = str(r.get("warhead_fragment_smiles", "")).strip()
        if not frag or frag == getattr(wl, "UNRESOLVED", "UNRESOLVED"):
            continue
        p = Chem.MolFromSmarts(frag)
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
        **IGNORED — see D0025.** The severed bond is filled with hydrogen
        regardless of what is passed here. That is why this function no longer
        feeds the gate: H-capping turns amides into formamides and thioethers
        into thiols, inventing alerts the molecule never had. The result is
        retained for human inspection of what the decoration looks like, and for
        nothing that decides anything. Use `attribute_alerts` to gate.
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
class AttributedAlerts:
    """Alerts on the INTACT molecule, attributed to core or decoration (D0025).

    Three buckets, because an alert can sit wholly inside the core, wholly
    outside it, or straddle the boundary:

    core      every matched atom is a core atom. Expected and excused — this is
              the warhead the approach is built around.
    rgroup    no matched atom is a core atom. The decoration's own alert.
    boundary  the match spans both. Counted against the decoration, because the
              motif would not exist without it, but reported separately so a
              reviewer can see this was a judgement call rather than a clean hit.
    """

    core: list[str] = field(default_factory=list)
    rgroup: list[str] = field(default_factory=list)
    boundary: list[str] = field(default_factory=list)
    core_found: bool = True

    @property
    def attributable(self) -> int:
        """Alerts the decoration is responsible for."""
        return len(self.rgroup) + len(self.boundary)

    def to_columns(self) -> dict:
        return {
            "core_alert_total": len(self.core),
            "core_alert_names": "|".join(self.core),
            "rgroup_alert_total": self.attributable,
            "rgroup_alert_names": "|".join(self.rgroup + self.boundary),
            "boundary_alert_total": len(self.boundary),
            "boundary_alert_names": "|".join(self.boundary),
        }


def attribute_alerts(smiles: str, core_smarts: str,
                     catalogs: tuple[str, ...] = CATALOGS,
                     *, excuse_warheads: bool = True) -> AttributedAlerts:
    """Screen the INTACT molecule and attribute each alert by its matched atoms.

    THIS DOES NOT CUT THE MOLECULE (D0025). The previous approach excised the
    R-group and capped the severed bond, which changed the fragment's chemistry:
    an amide ``>N-C(=O)-R`` severed from its nitrogen and capped with hydrogen
    becomes a formamide ``H-C(=O)-R``, and BRENK correctly flags an aldehyde
    that exists only in the fragment. The same mechanism turned thioethers into
    thiols. Together those two accounted for 3,014 of T_3's 5,270 rejections.

    Screening the intact molecule and asking RDKit which atoms each alert
    matched removes the failure mode rather than tuning around it: no bond is
    broken, so no functional group can be invented by breaking one.
    """
    out = AttributedAlerts()
    mol = smi.to_mol(smiles)
    patt = Chem.MolFromSmarts(core_smarts)
    if mol is None or patt is None:
        out.core_found = False
        return out
    core_match = mol.GetSubstructMatch(patt)
    if not core_match:
        out.core_found = False
        return out
    core_atoms = set(core_match)

    # THE WARHEAD IS EXCUSED ALONGSIDE THE CORE, and it has to be named
    # separately because it is not inside the core SMARTS. T_4 scopes alerts
    # against `N[CH]1CCS(=O)(=O)C1` — the sulfolane and its nitrogen — while the
    # warhead hangs off that nitrogen, outside the pattern. Attributing purely
    # by the core would blame every chloroacetamide's `alkyl_halide` on the
    # decoration and reject all 1,683 survivors: exactly the false positive the
    # two-tier design exists to prevent, arriving through a new door.
    #
    # Warhead atoms are located on the INTACT molecule, which is the same reason
    # the old excision path did it that way: the T_4 core includes the amide
    # nitrogen, so removing it splits an acrylamide in half and neither piece
    # still matches the Michael-acceptor pattern.
    # The WHOLE warhead group is excused, not just its reactive atoms. The
    # reactive-atom SMARTS is deliberately narrow — chloroacetamide's is
    # `[CH2][Cl]`, two atoms — so excusing only those leaves the carbonyl
    # exposed, and BRENK's `alpha_halo_carbonyl` then straddles the boundary and
    # is charged to the decoration. That motif IS the chloroacetamide warhead.
    # The library's `warhead_fragment_smiles` describes the whole group, so it
    # is what defines the excused region.
    if excuse_warheads:
        for _cls, wpatt in _warhead_group_patterns():
            for m in mol.GetSubstructMatches(wpatt):
                core_atoms.update(m)

    for cat_name in catalogs:
        for fm in _catalog(cat_name).GetFilterMatches(mol):
            atoms = {b for _, b in fm.atomPairs}
            name = fm.filterMatch.GetName()
            if atoms <= core_atoms:
                out.core.append(name)
            elif not (atoms & core_atoms):
                out.rgroup.append(name)
            else:
                out.boundary.append(name)
    return out


@dataclass
class TwoTierResult:
    """Whole-molecule (advisory) plus R-group (gating) alert results."""

    whole: AlertResult
    rgroup: AlertResult | None
    rgroup_smiles: str | None
    passes_gate: bool
    reason: str = ""
    attributed: "AttributedAlerts | None" = None
    excused: list = field(default_factory=list)

    def to_columns(self) -> dict:
        out = self.whole.to_columns(prefix="whole_")
        if self.attributed is not None:
            out.update(self.attributed.to_columns())
        else:
            out.update({"core_alert_total": 0, "core_alert_names": "",
                        "rgroup_alert_total": 0, "rgroup_alert_names": "",
                        "boundary_alert_total": 0, "boundary_alert_names": ""})
        # Retained for human inspection only; the gate no longer reads it (D0025).
        out["rgroup_smiles"] = self.rgroup_smiles or ""
        out["excused_alert_names"] = "|".join(self.excused)
        out["excused_alert_total"] = len(self.excused)
        out["alert_gate_pass"] = self.passes_gate
        out["alert_gate_reason"] = self.reason
        # The two-tier path DOES apply a real gate (max_rgroup_alerts), so a
        # True here means a threshold was tested and cleared. Contrast the
        # whole-molecule path — see screen_frame.
        out["alert_gate_applied"] = True
        return out


def two_tier(smiles: str, core_smarts: str, *,
             max_rgroup_alerts: int = 0,
             excused_alerts: tuple[str, ...] = ()) -> TwoTierResult:
    """Advisory whole-molecule alerts + a gating R-group check (covalent approaches).

    Parameters
    ----------
    max_rgroup_alerts : int
        Alerts tolerated on the isolated R-group. Zero by default: a PAINS or
        BRENK motif on the decoration has no mechanistic justification, unlike
        one on the warhead.
    """
    whole = screen(smiles)
    att = attribute_alerts(smiles, core_smarts)
    if not att.core_found:
        # Core absent. Not a pass: verify-after-expansion exists precisely to
        # catch products that lost the core.
        return TwoTierResult(whole=whole, rgroup=None, rgroup_smiles=None,
                             attributed=att, passes_gate=False,
                             reason="core not found in the molecule")

    # The R-group SMILES is still reported for human inspection, but it is NOT
    # what the gate reads any more (D0025). Excision is lossy at the cut bond;
    # attribution on the intact molecule is not.
    rg_smiles = isolate_rgroup(smiles, core_smarts)

    # EXCUSED ALERTS ARE CARRIED, NOT COUNTED (D0026). Naming an alert here is a
    # decision that this specific liability is understood and accepted for this
    # approach. It still travels with the candidate as `excused_alert_names` so
    # the GUI can show it — the D0019 pattern: flag, do not veto. Raising
    # `max_rgroup_alerts` instead would admit the named liability AND every
    # other one-off alert indiscriminately, which is a different decision
    # wearing the same clothes.
    names = att.rgroup + att.boundary
    excused = [n for n in names if n in excused_alerts]
    counted = [n for n in names if n not in excused_alerts]

    ok = len(counted) <= max_rgroup_alerts
    reason = ("" if ok else
              f"decoration carries {len(counted)} alert(s): "
              + ", ".join(counted)[:180]
              + (f" [{len([n for n in att.boundary if n not in excused_alerts])} "
                 "spanning the core boundary]" if att.boundary else ""))
    rg_result = AlertResult(counts={"attributed": len(counted)},
                            names={"attributed": counted})
    res = TwoTierResult(whole=whole, rgroup=rg_result, rgroup_smiles=rg_smiles,
                        attributed=att, passes_gate=ok, reason=reason)
    res.excused = excused
    return res


def screen_frame(df: pd.DataFrame, smiles_col: str = "canonical_smiles", *,
                 core_smarts: str | None = None,
                 disqualifying: tuple[str, ...] = (),
                 max_rgroup_alerts: int = 0,
                 excused_alerts: tuple[str, ...] = ()) -> pd.DataFrame:
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
            records.append(two_tier(s, core_smarts,
                                    max_rgroup_alerts=max_rgroup_alerts,
                                    excused_alerts=excused_alerts).to_columns())
        else:
            r = screen(s)
            cols = r.to_columns(prefix="whole_")
            # `alert_gate_pass = True` MUST NOT BE EMITTED WHEN NO GATE RAN.
            #
            # `disqualifying` defaults to () and is never passed by any caller,
            # so `bad` was always 0 and this column was `True` for every row —
            # 4,803 of 4,803 for T_1. It read as "passed the alert filter" and
            # meant "no filter was configured". T_1's rank-#10 molecule carries
            # catechol, ortho-hydroquinone and phosphor alerts and was `True`.
            #
            # PI decision 2026-07-31: report, do not gate. T_1's alerts stay
            # advisory — a catechol rule would reject EGCG, a known Pin1 binder,
            # and two earlier alert-derived rules died the same way. But the
            # column must say which it is, so `alert_gate_pass` is NA when no
            # gate ran and `alert_gate_applied` records the fact.
            applied = bool(disqualifying)
            bad = sum(r.counts.get(c, 0) for c in disqualifying)
            cols["alert_gate_applied"] = applied
            cols["alert_gate_pass"] = (bad == 0) if applied else pd.NA
            cols["alert_gate_reason"] = (
                "no disqualifying catalog configured — alerts are ADVISORY, "
                "not gated" if not applied else
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
