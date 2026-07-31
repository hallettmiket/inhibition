"""
Purpose: PROPOSED additions to shared/synthesizability.py, arising from the
         full workup of t1_db179d172dda (T_1 rank #10, SAscore 3.88).
Author:  blacksmith (with Mike Hallett)
Date:    2026-07-31
Input:   SMILES
Output:  rule definitions + the false-positive evidence that qualifies them

WHY THESE TWO, AND WHY NOT THE OBVIOUS ONES
-------------------------------------------
t1_db179d172dda passes all seven current rules and scores SAscore 3.88 -- the
LOWEST (easiest) of T_1's top 10. AiZynthFinder cannot solve it in 3,000
iterations / 30,438 nodes / 7,968 routes, and across 197 distinct molecules in
the retained trees it NEVER disconnects the cyclic acyl phosphate: 183 of 197
still carry it, and the 14 that do not are commodity reagents (POCl3, Ac2O,
mCPBA, BnBr, ...). So SAscore and the seven rules are both wrong on this
molecule, and they are wrong for a reason neither instrument can see.

THE FALSE-POSITIVE TEST (rules/D-rule: a rule that rejects a molecule somebody
has actually made is wrong by definition). Every candidate rule below was run
against all 19 parseable binders in data/reference/pin1_reference_binders_3.csv:

  rule                          fires_on_target   binders_killed   verdict
  acyl_phosphate                     yes                0          ADOPTABLE
  stereogenic_phosphorus             yes                0          ADOPTABLE
  -- rejected candidates, kept here so nobody re-proposes them --
  any phosphorus [P]                 yes                4          KILLED
  alkyl phosphate [CX4][OX2][PX4]    yes                4          KILLED
  catechol                           yes                1 (EGCG)   KILLED
  any 7-membered carbocycle          yes                0          NOT PROPOSED (see note)

The four binders killed by the naive phosphorus rules are
Wildemann-macrocyclic-peptide, Guo-Pfizer-benzothiophene-phosphonate,
Liu-Pei-cyclic-peptide and Jiang-Pei-bicyclic-CPP-peptide -- all pSer/pThr-Pro
mimetics that legitimately carry OP(=O)(O)O. Pin1 is a phosphate-binding
enzyme; a rule that bans phosphorus bans the target's own pharmacophore. The
acyl_phosphate SMARTS is safe precisely because it demands a CARBONYL carbon on
the ester oxygen, which no reference phosphate has.

MEASURED IMPACT (latest frame per approach):
  approach   n      acyl_phosphate   P_stereocentre   currently rejected
  T_1     4803           24              222                1284
  T_2     1882            0                0                  15
  T_3     5396            0                0                  44
  T_4     1782            0                0                   0
Both rules are T_1-only. T_2/T_3/T_4 are seeded from real molecules and
generate by reaction or enumeration, so neither motif ever arises. This is a
DiffSBDD-specific failure and the rules are cheap for everyone else.

SAscore does not see it: median SAscore of the 24 acyl-phosphate T_1 molecules
is 4.38 against 4.14 for T_1 as a whole. The motif is invisible to SAscore
because SAscore scores ECFP4 fragment RARITY, and every radius-2 environment
here (aryl-OMe, aryl-O-aryl, C-O-P, P=O, phenol) is individually common. It is
the COMBINATION that has no precedent, and a fragment-additive score cannot
represent a combination.
"""

from __future__ import annotations

from rdkit import Chem

# ---------------------------------------------------------------------------
# RULE 1 -- add to RULES in shared/synthesizability.py
# ---------------------------------------------------------------------------
ACYL_PHOSPHATE_RULE = dict(
    name="acyl_phosphate",
    smarts="[CX3](=[OX1])[OX2][PX4]",
    why=(
        "A mixed carboxylic-phosphoric anhydride (acyl phosphate). This is a "
        "high-energy acyl-transfer group -- biology uses it precisely because "
        "it is transient (acetyl phosphate, aminoacyl-adenylate, the E2-P "
        "aspartyl phosphate of the Ca2+-ATPase). It cannot be carried through "
        "a synthesis, chromatographed, stored, or held intact in an assay "
        "buffer. NOTE the CARBONYL in the pattern: ordinary phosphate esters "
        "and phosphonates are NOT matched, so the pSer/pThr-Pro mimetics that "
        "make up a third of the Pin1 reference set survive this rule."
    ),
)


# ---------------------------------------------------------------------------
# RULE 2 -- needs a function, not a SMARTS: stereogenicity is not substructural
# ---------------------------------------------------------------------------
def n_stereogenic_phosphorus(mol: Chem.Mol) -> int:
    """Count stereogenic phosphorus atoms (assigned or not).

    Uses FindPotentialStereo rather than FindMolChiralCenters: CIP labelling
    raises `Digraph generation failed: more than 100000 nodes` on the large
    peptidic candidates in the reference set.

    WHY THIS IS A SYNTHESIZABILITY RULE AND NOT A NICE-TO-HAVE. A generator
    that emits a P-stereocentre is writing down a configuration it has no
    means to deliver. Stereocontrolled P is hard even where the industry has
    spent the most money on it -- the ProTide/oligonucleotide phosphoramidates
    (sofosbuvir, remdesivir) needed chiral auxiliaries, asymmetric catalysis or
    engineered phosphotriesterase biocatalysts to reach a single diastereomer.
    Unlike carbon, there is no general asymmetric method. Treat an unrequested
    P-stereocentre as a route the pipeline cannot cost.
    """
    try:
        elements = Chem.FindPotentialStereo(mol)
    except Exception:
        return 0
    return sum(
        1
        for e in elements
        if e.type == Chem.StereoType.Atom_Tetrahedral
        and mol.GetAtomWithIdx(e.centeredOn).GetSymbol() == "P"
    )


STEREOGENIC_PHOSPHORUS_RULE = dict(
    name="stereogenic_phosphorus",
    predicate=n_stereogenic_phosphorus,
    max_allowed=0,
    why=n_stereogenic_phosphorus.__doc__.split("\n\n")[-1].strip(),
)


# ---------------------------------------------------------------------------
# NOT PROPOSED, and the reason recorded so it is not re-litigated
# ---------------------------------------------------------------------------
NOT_PROPOSED = {
    "seven_membered_carbocycle": (
        "[#6;r7]",
        "523 of 4803 T_1 candidates (10.9%) contain one, and 4 of T_1's top 10 "
        "do -- a real and suspicious signal that DiffSBDD's 3D->graph bond "
        "inference closes 7-rings where 6-rings belong. But it must NOT become "
        "a rule: benzodiazepines, tropanes and colchicine all have 7-rings, and "
        "this molecule's own 6-5-7 core is cyclohepta[b]indole (PubChem CID "
        "251959, CAS 246-06-0, 1-azabenz[b]azulene), a real named ring system "
        "present in the ZINC purchasable stock. The right response is to "
        "INSTRUMENT the 7-ring rate as a generator-health metric, not to filter "
        "on it."
    ),
    "catechol": (
        "[OX2H]c1ccccc1[OX2H]",
        "Kills EGCG. Already caught and reported by shared/alerts.py "
        "('catechol_A(92)|catechol|ortho_hydroquinone'). The real gap is that "
        "for T_1 the alert GATE never fires -- alert_gate_pass is True for all "
        "4803 rows because the gate only applies to approaches passing "
        "core_smarts (D0026). Fix the gate's scope, do not duplicate the alert "
        "as a synthesizability rule."
    ),
    "any_phosphorus": ("[P]", "Kills 4 of 19 reference binders. See header."),
}
