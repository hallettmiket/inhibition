"""
Purpose: Re-filter a shortlist from plain-language chemistry constraints.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: a shortlist frame + a constraint string
Output: the filtered frame, plus what each constraint removed

Issue #1, general note 1: "we have a ranked list of candidates but we meet and
determine we dont like chlorines anywhere at all or less conjugation. This is
given to a curator agent that presents a new ranked list."

FILTERING IS NOT RE-RANKING, AND THE DIFFERENCE MATTERS HERE. This removes rows
and preserves the original order; it never computes a new score. The gate has
measured that the underlying ranking does not demonstrably enrich (D0041) and is
partly a size ranking (D0043) -- inventing a fresh score on top of that would
add a second unvalidated ordering to an unvalidated one. Dropping molecules a
chemist has ruled out is safe regardless of whether the score works.

EVERY CONSTRAINT REPORTS ITS OWN COUNT. A filter that silently removes 24 of 25
looks the same as one that removes nothing until you notice the table is empty,
so each rule states what it took.

UNPARSEABLE CONSTRAINTS ARE REFUSED, NOT IGNORED. A typo'd SMARTS that quietly
matches nothing reads as "no molecules have that group", which is the same
output as a working filter finding none -- and the user would act on it.

THE FILTER APPLIES TO EVERY CANDIDATE VIEW, NOT JUST THE ONE IT WAS TYPED INTO
(issue #3.2: "the curation feature does not carry over through the dossier and
rest of gui"). It used to be a text box inside the shortlist panel that wrote
`st.session_state["_curate_spec"]` and was then read by nobody: a chemist who
excluded chlorines still met chlorinated molecules in the dossier, in the
convergence pairs and in the axis medians. Worse than not having the filter,
because the shortlist panel had already told them the compounds were gone.

WHICH PANELS IT MUST *NOT* TOUCH IS PART OF THE FEATURE. `PANEL_SCOPE` below
declares both, with a reason for each exclusion, so "the dossier is filtered but
the synthesizability counts are not" is a decision on the record rather than an
inconsistency someone finds later. The rule that generates the list: a curation
filter expresses which molecules the chemist is willing to CONSIDER, so it
belongs on every view of the candidate set -- and it says nothing whatsoever
about what the pipeline generated, docked or measured, so it must never touch a
count, a rank denominator or a gate statistic that reports on that population.
Filtering those would turn "ranked 3rd of 1,204 docked" into "3rd of 46", which
is not a curated version of the fact; it is a false one.
"""

from __future__ import annotations

#: Mtime of THIS file at the moment it was imported. Frozen at import, so
#: comparing it with the file's current mtime is the only reliable way to tell
#: that a running process is executing stale code -- Streamlit re-runs the
#: script on every interaction but never re-imports helper modules.
LOADED_MTIME = __import__("os").stat(__file__).st_mtime

import re
from dataclasses import dataclass
from typing import Iterable

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

#: Plain-language names a chemist is likely to type, mapped to SMARTS.
NAMED_GROUPS: dict[str, str] = {
    "chlorine": "[Cl]",
    "chloro": "[Cl]",
    "fluorine": "[F]",
    "bromine": "[Br]",
    "iodine": "[I]",
    "halogen": "[F,Cl,Br,I]",
    "nitro": "[N+](=O)[O-]",
    "nitrile": "C#N",
    "azide": "[N-]=[N+]=[N-]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "ester": "[CX3](=O)[OX2H0][#6]",
    "sulfonamide": "[SX4](=O)(=O)[NX3]",
    "phosphonate": "[PX4](=O)([OX2])[OX2]",
    "michael_acceptor": "[CX3]=[CX3][CX3]=O",
    "epoxide": "C1OC1",
    "thiol": "[SX2H]",
    "primary_amine": "[NX3;H2][CX4]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "quaternary_carbon": "[CX4]([#6])([#6])([#6])[#6]",
}

#: Numeric properties that can be bounded.
PROPERTIES = {
    "mw": lambda m: Descriptors.MolWt(m),
    "logp": lambda m: Descriptors.MolLogP(m),
    "heavy_atoms": lambda m: m.GetNumHeavyAtoms(),
    "rotatable_bonds": lambda m: Descriptors.NumRotatableBonds(m),
    "hbd": lambda m: Descriptors.NumHDonors(m),
    "hba": lambda m: Descriptors.NumHAcceptors(m),
    "rings": lambda m: Descriptors.RingCount(m),
    "aromatic_rings": lambda m: Descriptors.NumAromaticRings(m),
    "tpsa": lambda m: Descriptors.TPSA(m),
    # Conjugation, asked for by name in issue #1. Counted as bonds between two
    # sp2 atoms -- a coarse proxy, and named coarsely so nobody reads it as a
    # computed conjugation length.
    "sp2_bonds": lambda m: sum(
        1 for b in m.GetBonds()
        if str(b.GetBeginAtom().GetHybridization()) == "SP2"
        and str(b.GetEndAtom().GetHybridization()) == "SP2"),
}


class ConstraintError(ValueError):
    """A constraint could not be understood, so nothing was filtered."""


@dataclass
class Rule:
    text: str
    kind: str          # "exclude_group" | "require_group" | "bound"
    detail: str
    removed: int = 0


def parse(spec: str) -> list[Rule]:
    """Turn a newline- or semicolon-separated constraint string into rules.

    Accepted forms, one per line:
        no chlorine                 exclude anything matching that group
        no [Cl]                     same, with explicit SMARTS
        require sulfonamide         keep only molecules that have it
        mw < 450                    numeric bound on a listed property
        sp2_bonds <= 8              "less conjugation", coarsely
    """
    rules: list[Rule] = []
    for raw in re.split(r"[\n;]+", spec or ""):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        # SMARTS IS CASE-SENSITIVE AND MUST NOT BE LOWERCASED. `[Cl]` is
        # chlorine; `[cl]` is an aromatic chlorine that matches nothing. Only
        # the leading KEYWORD is matched case-insensitively -- the argument
        # keeps the case the user typed. Lowercasing the whole line turned
        # "no [Cl]" into a pattern that silently matched no molecules, which
        # reads identically to a working filter finding none.
        m = re.match(r"^(no|exclude|remove)\s+(.+)$", line, re.IGNORECASE)
        if m:
            rules.append(Rule(line, "exclude_group", _as_smarts(m.group(2), line)))
            continue
        m = re.match(r"^(require|only|keep)\s+(.+)$", line, re.IGNORECASE)
        if m:
            rules.append(Rule(line, "require_group", _as_smarts(m.group(2), line)))
            continue
        low = line.lower()
        m = re.match(r"^(\w+)\s*(<=|>=|<|>|==)\s*(-?\d+(?:\.\d+)?)$", low)
        if m:
            prop, op, val = m.group(1), m.group(2), m.group(3)
            if prop not in PROPERTIES:
                raise ConstraintError(
                    f"unknown property {prop!r} in {line!r}. "
                    f"Available: {', '.join(sorted(PROPERTIES))}")
            rules.append(Rule(line, "bound", f"{prop}{op}{val}"))
            continue
        raise ConstraintError(
            f"could not parse {line!r}. Try 'no chlorine', "
            "'require sulfonamide', or 'mw < 450'.")
    return rules


def _as_smarts(token: str, line: str) -> str:
    """Resolve a group name or SMARTS. NAMES are case-insensitive, SMARTS is not."""
    token = token.strip()
    key = token.lower().replace(" ", "_").replace("-", "_")
    if key in NAMED_GROUPS:
        return NAMED_GROUPS[key]
    if Chem.MolFromSmarts(token) is not None:
        return token
    raise ConstraintError(
        f"{token!r} in {line!r} is neither a known group name nor valid "
        f"SMARTS. Known names: {', '.join(sorted(NAMED_GROUPS))}")


def apply(df: pd.DataFrame, spec: str,
          smiles_col: str = "canonical_smiles") -> tuple[pd.DataFrame, list[Rule]]:
    """Filter `df`, preserving its existing order. Returns (kept, rules)."""
    rules = parse(spec)
    if not rules:
        return df, []
    if smiles_col not in df.columns:
        raise ConstraintError(f"no {smiles_col!r} column to filter on")

    mols = [Chem.MolFromSmiles(str(s)) for s in df[smiles_col]]
    keep = pd.Series(True, index=df.index)

    for rule in rules:
        before = int(keep.sum())
        if rule.kind in ("exclude_group", "require_group"):
            patt = Chem.MolFromSmarts(rule.detail)
            hit = pd.Series(
                [bool(m and m.HasSubstructMatch(patt)) for m in mols],
                index=df.index)
            keep &= (~hit if rule.kind == "exclude_group" else hit)
        else:
            prop, op, val = re.match(r"^(\w+)(<=|>=|<|>|==)(-?[\d.]+)$",
                                     rule.detail).groups()
            fn, val = PROPERTIES[prop], float(val)
            vals = pd.Series(
                [(fn(m) if m else float("nan")) for m in mols], index=df.index)
            cmp = {"<": vals < val, "<=": vals <= val, ">": vals > val,
                   ">=": vals >= val, "==": vals == val}[op]
            keep &= cmp.fillna(False)
        rule.removed = before - int(keep.sum())

    return df[keep], rules


# --------------------------------------------------------------------------
# where the filter applies, and where applying it would be wrong
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Scope:
    """Whether one GUI panel honours the curation filter, and why."""

    panel: str
    filtered: bool
    why: str


#: Every panel in the app, stated explicitly. A panel absent from this table is
#: a panel nobody decided about, so `scope_for` raises rather than guessing —
#: silently defaulting either way is how the original bug survived.
PANEL_SCOPE: tuple[Scope, ...] = (
    Scope("Shortlists", True,
          "The candidate lists themselves. This is where the constraints are "
          "typed and the only place they used to be honoured."),
    Scope("T₂ seed comparison", False,
          "This panel asks whether the SEED determines what the arm produces, "
          "which is a property of the generated population. A curation filter "
          "says which molecules the chemist would consider, and applying it "
          "here would quietly turn the question into 'which seed best suits my "
          "current constraint' — a different and much weaker claim. The "
          "per-seed medians and pool sizes must describe what CReM actually "
          "generated, for the same reason a rank denominator must."),
    Scope("Near-attack ranking", True,
          "A candidate list, so the chemist's constraints apply for the same "
          "reason they do on Shortlists — a molecule that has been ruled out "
          "should not be offered at the top of a second ranking as though it "
          "were still in contention. The panel's own caveats are about how much "
          "the ORDER can be trusted (the score does not converge, D0068); "
          "curation is about which molecules belong in the list at all, and the "
          "two are independent. The enrichment values are NOT recomputed on the "
          "filtered subset: each is a per-molecule measurement against an "
          "isotropic baseline, not a rank within the displayed set, so it means "
          "the same thing whoever else is on screen."),
    Scope("Candidate dossier", True,
          "A view of one shortlisted candidate. A molecule the chemist has "
          "ruled out should not be offered for inspection as though it were "
          "still in contention — though it stays reachable behind an explicit "
          "toggle, because 'why did you drop this one' is a fair question."),
    Scope("Pose clusters", False,
          "This panel describes the GEOMETRY of the modes a docking run "
          "produced for a candidate — how many distinct binding modes there "
          "are, and which real mode is their medoid. A curation filter says "
          "which molecules a chemist would consider; it says nothing about how "
          "many ways a molecule was posed. Filtering here would change "
          "'7 clusters of 9 modes' into a statement about a subset of modes "
          "that were never separately docked, which is not a curated version "
          "of the fact but a false one — the same reason a rank denominator is "
          "left alone. The CANDIDATE SELECTOR is a different question and may "
          "be scoped later if the shortlist filter should drive it; that would "
          "be a decision, not a default."),
    Scope("Convergence", True,
          "Cross-approach pairs are drawn FROM the shortlists, so a pair whose "
          "members are both excluded is not evidence about anything the "
          "chemist is still considering."),
    Scope("Shared axes", True,
          "Medians over the shortlists. Filtered — but a numeric constraint on "
          "an axis truncates that axis by construction, so the panel says so "
          "rather than reporting a median the filter created."),
    Scope("Within-stratum", True,
          "Two leaderboards over shortlisted candidates; same argument as the "
          "shortlists."),
    Scope("Decisions", False,
          "The choreography's decision log. It records why the pipeline is "
          "built the way it is and has no candidates in it."),
    Scope("Why this file?", False,
          "Decision records matching a path fragment. No candidates."),
    Scope("Provenance", False,
          "Run manifests — what code produced what output. Curating these "
          "would misrepresent what actually ran."),
    Scope("Open questions", False,
          "Unresolved decisions and unverified sources. Not a candidate view."),
)

_SCOPE_BY_PANEL = {s.panel: s for s in PANEL_SCOPE}


#: Facts about the FULL population that must survive curation untouched, quoted
#: in the GUI wherever one of them is displayed next to a curated table. These
#: are the answers to "where would filtering be wrong", and they are wrong for
#: the same reason in each case: they are counts over what the pipeline
#: produced, not over what the chemist is willing to consider.
UNFILTERED_FACTS: tuple[tuple[str, str], ...] = (
    ("synthesizability delta",
     "`shortlist_delta` counts how many candidates the synthesizable rebuild "
     "dropped and promoted across the approach's whole quota. Curation happens "
     "downstream of that rebuild and cannot change what it did."),
    ("rank denominators",
     "`rank` and `group_n_docked` — 'ranked 3rd of 1,204 docked' — describe "
     "the docked population. Re-deriving them after a filter would produce a "
     "number that was never true."),
    ("cross-approach convergence lookup",
     "Whether another approach ALSO ranked this molecule is a fact about the "
     "pipeline's output, not about the chemist's preferences, so the lookup "
     "runs over every ranked row regardless of the filter."),
    ("enrichment gate verdicts",
     "ROC-AUC, its CI and EF1% are measured on known actives against "
     "property-matched decoys. No candidate of ours is in that calculation."),
)


def scope_for(panel: str) -> Scope:
    """The declared scope for one panel. Raises on an undeclared panel."""
    try:
        return _SCOPE_BY_PANEL[panel]
    except KeyError:
        raise KeyError(
            f"panel {panel!r} has no curation scope declared. Add it to "
            "curate.PANEL_SCOPE with a reason — a panel that shows candidates "
            "must not default to unfiltered by omission.") from None


#: curate property -> the shared physicochemical axis it constrains. Every pair
#: here is the SAME RDKit call on both sides (shared/descriptors.py computes MW
#: with Descriptors.MolWt, HAC with GetNumHeavyAtoms, cLogP with Crippen.MolLogP
#: and TPSA with CalcTPSA; PROPERTIES above calls exactly those), so a bound on
#: the left truncates the column on the right exactly, not approximately.
PROPERTY_AXIS: dict[str, str] = {
    "mw": "MW",
    "heavy_atoms": "HAC",
    "logp": "cLogP",
    "tpsa": "TPSA",
}


def bounded_axes(rules: Iterable[Rule]) -> dict[str, str]:
    """{axis column: the rule text that bounds it} for numeric constraints.

    Used by the shared-axes panel to refuse to present a median as a property
    of an approach when the filter is what set its limit.
    """
    out: dict[str, str] = {}
    for rule in rules:
        if rule.kind != "bound":
            continue
        m = re.match(r"^(\w+)", rule.detail)
        axis = PROPERTY_AXIS.get(m.group(1)) if m else None
        if axis:
            out[axis] = rule.text
    return out


def describe(rules: Iterable[Rule]) -> str:
    """The active constraints as one compact line, each with its own count."""
    rules = list(rules)
    if not rules:
        return "no constraints"
    return " · ".join(f"`{r.text}` −{r.removed}" for r in rules)


def banner(rules: Iterable[Rule], n_before: int, n_after: int,
           panel: str | None = None) -> str:
    """The persistent indicator: what is filtered, and how much it removed.

    Deliberately states the KEPT count and the REMOVED count both. A filter
    that removes 24 of 25 and one that removes nothing produce tables that look
    equally plausible, and only the numbers separate them.
    """
    rules = list(rules)
    removed = n_before - n_after
    where = f" in {panel}" if panel else ""
    if not rules:
        return f"No curation filter active — showing all {n_before}{where}."
    return (f"**Curation active{where}: showing {n_after} of {n_before}** "
            f"({removed} hidden) — {describe(rules)}")
