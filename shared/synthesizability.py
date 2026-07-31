"""
Purpose: Reject molecules that cannot be made, using structural rules only.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: SMILES
Output: a list of violated rules, empty when nothing is obviously wrong

Issue #1, T_1 note: "We cannot rely on computed synthesizability metrics as the
current state of these tools are not accurate. Rather we need to consider
basic/simple rules or parameters to quickly eliminate impossible to synthesize
molecules (for example, a few molecules have joined carbon rings with multiple
groups attached that cannot be synthesized together and joined)."

WHY NOT SA SCORE. It is already computed for every candidate and consumed by
nothing, deliberately. An SA score is a statistical resemblance to known
synthetic chemistry, not a claim about whether a route exists, and it moves
smoothly -- so any threshold on it removes borderline-unusual molecules along
with impossible ones. These rules are the opposite: each is a specific,
nameable structural impossibility, and a molecule either has it or does not.

THE RULES ARE A STARTING SET AND ARE MEANT TO BE ARGUED WITH. Each carries a
one-line rationale so a chemist can disagree with a specific rule rather than
with "the filter". Adding a rule is adding a row here.

CONSERVATIVE BY CONSTRUCTION. Every rule is checked against the known Pin1
binders in tests/test_synthesizability.py, and a rule that rejects a molecule
somebody has actually made is wrong by definition -- the compound exists. A
filter is only worth having if its false-positive rate is near zero, because
the cost of dropping a good candidate is invisible and permanent.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


@dataclass(frozen=True)
class Rule:
    name: str
    smarts: str
    why: str
    max_allowed: int = 0


#: Each rule fires when its SMARTS matches more than `max_allowed` times.
RULES: tuple[Rule, ...] = (
    Rule("anti_bredt_bridgehead_alkene",
         "[R2]=[R2]",
         "A double bond at a bridgehead of a fused/bridged system violates "
         "Bredt's rule below ~8 ring atoms; the alkene cannot be planar."),
    Rule("adjacent_quaternary_ring_carbons",
         # `!@` = the bond to the substituent must be NON-RING. Without it the
         # pattern counts a carbon's own ring bonds as substituents, so every
         # ordinary fused bicyclic matches -- it flagged 47% of T_1 and rejected
         # EGCG, a molecule found in tea.
         "[CX4;R](!@[!#1])(!@[!#1])[CX4;R](!@[!#1])(!@[!#1])",
         "Two neighbouring ring carbons each carrying two EXOCYCLIC heavy "
         "substituents. This is the case flagged in issue #1: joined carbon "
         "rings with multiple groups that cannot be installed together."),
    Rule("peroxide_or_higher",
         "[OX2][OX2]",
         "O-O bonds are explosive hazards and are not installed casually."),
    Rule("geminal_diol_or_hemiketal",
         "[CX4]([OX2H])([OX2H])",
         "A geminal diol is the hydrate of a carbonyl and is not isolable in "
         "the form drawn."),
    Rule("nitrogen_nitrogen_nitrogen_chain",
         "[NX3][NX3][NX3]",
         "A non-aromatic N-N-N chain outside a tetrazole/azide is unstable."),
    Rule("more_than_four_contiguous_heteroatoms",
         "[!#6;!#1][!#6;!#1][!#6;!#1][!#6;!#1]",
         "Four heteroatoms in a row is almost always a generator artefact."),
    Rule("strained_small_ring_fusion",
         "[C;r3;R2]",
         "A cyclopropane carbon shared with a second ring. Bicyclobutane-like "
         "strain; makeable only by dedicated routes, not incidentally."),
)


def violations(smiles: str) -> list[Rule]:
    """Which rules this molecule breaks. Empty means nothing obvious is wrong.

    Empty is NOT a claim that the molecule is easy to make. It says only that
    none of the named impossibilities is present, which is the only claim these
    rules can support.
    """
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else smiles
    if m is None:
        return [Rule("unparseable", "", "not a valid molecule")]
    out = []
    for rule in RULES:
        patt = Chem.MolFromSmarts(rule.smarts)
        if patt is None:
            continue
        if len(m.GetSubstructMatches(patt)) > rule.max_allowed:
            out.append(rule)
    return out


def is_plausible(smiles: str) -> bool:
    """True when no rule fires. See `violations` for what that does not mean."""
    return not violations(smiles)


def explain(smiles: str) -> str:
    v = violations(smiles)
    if not v:
        return "no structural red flags"
    return "; ".join(f"{r.name}: {r.why}" for r in v)
