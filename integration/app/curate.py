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
"""

from __future__ import annotations

import re
from dataclasses import dataclass

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
