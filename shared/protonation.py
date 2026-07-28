"""
Purpose: Put every ligand in its dominant protonation state at physiological pH.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: a candidate SMILES
Output: the SMILES of the dominant microspecies at pH 7.4, and its formal charge

WHY THIS EXISTS. Generators emit neutral SMILES. `CC(=O)O` is a carboxylic
ACID, and at pH 7.4 a carboxylic acid (pKa ~4-5) is more than 99% deprotonated —
it is a carboxylATE. Scoring the neutral form is scoring a species that is not
there, which is the same class of error as D0022 (docking the pre-reaction
ligand) and D0030 (docking an unsaturated Michael adduct): the question is
always which species is actually being modelled.

WHY IT MATTERS HERE SPECIFICALLY, AND A LOT. 14 of T_2's 25 shortlisted
candidates carry a -COOH, and Pin1's pocket is strongly cationic — measured from
Cys113 SG in our own receptor, Arg69 sits at 6.5 A, Lys63 at 7.3 A and Arg68 at
11.6 A. Two errors compound if the acid is left neutral, and they do NOT cancel:

  1. the electrostatic attraction to that basic cluster is omitted, and
  2. so is the desolvation penalty the real anion pays — the hydration free
     energy of acetate is roughly -80 kcal/mol against about -7 for acetic acid.

MM-GBSA is where this bites hardest. antechamber assigns AM1-BCC charges under
the total charge it is GIVEN, so passing a carboxylate as neutral does not merely
omit a charge; it redistributes a whole electron across the molecule.

METHOD, AND ITS LIMITS. Open Babel's `-p` transformation model applies pKa rules
to the standard ionizable groups. It is a rules engine, not a pKa predictor: it
does not know local environment, it will not shift a pKa because a neighbouring
group perturbs it, and it takes no view on tautomers. For the groups that
dominate here — carboxylic acids, aliphatic amines, guanidines — the dominant
state at 7.4 is not in doubt, and a rules engine gets it right. For anything
with a genuinely borderline pKa the answer is a measurement, not a better
default, so `is_confident()` reports which is which rather than hiding it.

THE POSE IS UNAFFECTED. Changing protonation adds or removes hydrogens; the
heavy-atom skeleton is identical. A pose docked in the neutral form therefore
remains usable as geometry. What it does NOT license is claiming the DOCKING saw
the right species — Vina scored the neutral form, and that is recorded rather
than quietly corrected.
"""

from __future__ import annotations

import logging
import subprocess

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

OBABEL = "/data/lab_vm/envs/dwi_cheminf/bin/obabel"
PHYSIOLOGICAL_PH = 7.4

# Groups whose dominant state at pH 7.4 is unambiguous. Anything ionizable that
# is NOT on this list gets flagged rather than trusted.
CONFIDENT_GROUPS = {
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "carboxylate": "[CX3](=O)[OX1-]",
    "aliphatic_amine": "[NX3;H2,H1;!$(NC=O);!$(Nc)]",
    "ammonium": "[NX4+]",
    "guanidine": "[NX3][CX3](=[NX2])[NX3]",
    "guanidinium": "[NX3][CX3](=[NX2+])[NX3]",
    "sulfonic_acid": "[SX4](=O)(=O)[OX2H1]",
    "phosphate": "[PX4](=O)([OX2H1])",
    "tetrazole": "c1nnn[nH]1",
}

# Ionizable in principle, but with a pKa near enough to 7.4 that the dominant
# state is a real question rather than a lookup.
BORDERLINE_GROUPS = {
    "imidazole": "c1cnc[nH]1",
    "aromatic_amine": "[NX3;H2,H1][c]",
    "thiol": "[SX2H1]",
    "phenol": "[OX2H1][c]",
}


def _match(smiles: str, patterns: dict[str, str]) -> list[str]:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return []
    out = []
    for name, sma in patterns.items():
        p = Chem.MolFromSmarts(sma)
        if p is not None and m.HasSubstructMatch(p):
            out.append(name)
    return out


def is_confident(smiles: str) -> tuple[bool, list[str]]:
    """Whether the dominant state at pH 7.4 is a lookup or a judgement call."""
    borderline = _match(smiles, BORDERLINE_GROUPS)
    return (not borderline), borderline


def dominant_state(smiles: str, ph: float = PHYSIOLOGICAL_PH) -> dict:
    """The dominant microspecies at `ph`, with its formal charge.

    Returns a dict rather than a bare SMILES because the CHANGE is the
    interesting part: a stage that silently swaps a neutral ligand for an anion
    has changed what it is scoring, and that has to be visible in the frame.
    """
    neutral = Chem.MolFromSmiles(smiles)
    if neutral is None:
        raise ValueError(f"unparseable SMILES {smiles!r}")
    q0 = Chem.GetFormalCharge(neutral)

    proc = subprocess.run(
        [OBABEL, f"-:{smiles}", "-osmi", "-p", str(ph)],
        capture_output=True, text=True, timeout=120)
    out = (proc.stdout or "").strip().split("\t")[0].strip()
    m = Chem.MolFromSmiles(out) if out else None
    if m is None:
        log.warning("obabel returned no usable state for %s at pH %s; keeping "
                    "the input unchanged", smiles, ph)
        return {"protonated_smiles": Chem.MolToSmiles(neutral),
                "protonated_charge": q0, "charge_changed": False,
                "protonation_confident": False,
                "protonation_note": "obabel produced no state; input kept"}

    q = Chem.GetFormalCharge(m)
    confident, borderline = is_confident(smiles)
    if not confident:
        log.info("%s carries borderline group(s) %s; its state at pH %s is a "
                 "judgement call, not a lookup", smiles, borderline, ph)
    return {
        "protonated_smiles": Chem.MolToSmiles(m),
        "protonated_charge": q,
        "charge_changed": bool(q != q0),
        "protonation_confident": bool(confident),
        "protonation_note": ("; ".join(borderline) if borderline else ""),
    }


def protonate_frame(df, smiles_col: str = "canonical_smiles", ph: float = PHYSIOLOGICAL_PH):
    """Add protonation columns to a frame. Row count is unchanged."""
    rows = []
    for s in df[smiles_col]:
        try:
            rows.append(dominant_state(s, ph))
        except Exception as exc:  # noqa: BLE001 - one bad SMILES must not stop the frame
            log.warning("protonation failed for %r: %s", s, exc)
            rows.append({"protonated_smiles": s, "protonated_charge": 0,
                         "charge_changed": False, "protonation_confident": False,
                         "protonation_note": f"failed: {exc}"[:120]})
    import pandas as pd

    add = pd.DataFrame(rows, index=df.index)
    out = df.copy()
    for c in add.columns:
        out[c] = add[c]
    n = int(add["charge_changed"].sum())
    log.info("protonation at pH %s: %d/%d candidates change charge", ph, n, len(out))
    return out
