"""
Purpose: the charge a molecule actually carries at pH 7.4, and whether it carries a phosphate.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: canonical SMILES
Output: `charge_ph74` (int) and `has_phosphate` (bool) per molecule

#6 items 5 and 7 -- "T_2 phosphate: label, not protect" and charge-stratified
ranking. Both were decided and neither was implemented, and they land together
because the phosphate label is only honest once phosphate-free molecules are
actually evaluated rather than filtered out by charge stratification.

## Why `formal_charge` could not be used, even though it exists

`shared/descriptors.py` already computes `formal_charge` -- and it is **0 for
essentially every molecule in the project**: all 4,803 T_1 rows, all 1,882 T_2
rows, 5,378 of 5,396 T_3. That is not a bug in the descriptor. It is
`Chem.GetFormalCharge()` on the NEUTRAL canonical SMILES, which is exactly what
it says it is.

But docking protonated for pH 7.4 (`obabel -p 7.4`), so the molecule that was
SCORED is not the molecule that column describes. Stratifying on it would put
every candidate in one stratum and look like it had worked -- a populated,
plausible column that does not mean what its name suggests to a reader who
wants "the charge". Measured on the five T_2 seeds, all of which are
`formal_charge = 0`:

    atra          -1   carboxylate
    du_xu         -1   carboxylate
    guo_pfizer    -2   phosphate dianion
    potter_astex  +1   ammonium
    sulfopin       0   neutral

Five seeds, four distinct charge states, and the existing column reports one.

## Why obabel and not a pKa model

The charge must be the one the DOCKED structure carried, or the stratification
describes a different molecule from the one being ranked. `obabel -p 7.4` is
what `noncovalent_dock_run` used, so calling the same tool makes the two agree
by construction rather than by a model of mine that would need its own
validation. It is also fast enough to be uninteresting: 1,882 molecules per
second in batch.

## Identity, not position

obabel is invoked with a `SMILES<TAB>id` file and the result is matched back
**by the id it echoes**, never by line order. A molecule obabel drops -- and it
does drop some -- would otherwise shift every subsequent charge onto the wrong
candidate, which is the failure this project has catalogued more than any
other. Missing ids come back as None so the caller can see them, rather than
being silently filled with a plausible zero.
"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from pathlib import Path

log = logging.getLogger(__name__)

OBABEL = "/data/lab_vm/envs/dwi_cheminf/bin/obabel"

#: The pH docking used. Imported by callers rather than restated.
DEFAULT_PH = 7.4

# Phosphate / phosphonate ester or acid: a tetravalent P with a double-bonded O
# and at least two further oxygens. Deliberately covers the acid, the mono- and
# di-anion, and the ester, because #6 item 5 asks whether the MOLECULE carries
# a phosphate, not what its protonation state happens to be in one drawing.
PHOSPHATE_SMARTS = "[PX4](=O)([OX2,OX1-])([OX2,OX1-])"

_BATCH = 5000


def protonate(smiles_by_id: dict[str, str], ph: float = DEFAULT_PH,
              ) -> dict[str, str]:
    """{id: protonated SMILES} at `ph`, via the same obabel call docking used.

    Ids obabel does not return are ABSENT from the result rather than defaulted.
    """
    out: dict[str, str] = {}
    todo = list(smiles_by_id.items())
    # ONE BAD MOLECULE HALTS THE WHOLE STREAM, SILENTLY. Measured 2026-08-05 on
    # T_1: obabel reported "143 molecules converted" from a 4,803-line file,
    # wrote NOTHING to stderr, and exited 0. The 144th molecule is
    # `OC[P@TB14](O)(O)(O)CO` -- pentavalent phosphorus with trigonal-
    # bipyramidal stereochemistry, a DiffSBDD artefact -- and everything after
    # it was lost. 97% of the arm, with no error of any kind.
    #
    # Per-molecule conversion (what `noncovalent_dock_run._prepare_one` does)
    # is immune, which is why docking never hit this. Batch conversion is not.
    #
    # So: shrink and retry. A chunk that comes back short is split until the
    # failures are isolated to single molecules, which are then reported as
    # missing rather than dragging their neighbours down with them. Identity
    # matching is what makes this detectable at all -- by position, a short
    # return is invisible and every subsequent charge lands on the wrong
    # candidate.
    out.update(_convert_recursive(todo, ph, depth=0))
    missing = len(smiles_by_id) - len(out)
    if missing:
        log.warning("obabel returned no structure for %d of %d molecules; "
                    "their charge is None, not 0", missing, len(smiles_by_id))
    return out


def _convert_recursive(items: list[tuple[str, str]], ph: float,
                       depth: int) -> dict[str, str]:
    """Convert a chunk; on a short return, split it and retry the halves."""
    out = _convert_chunk(items, ph)
    if len(out) == len(items) or len(items) == 1 or depth > 24:
        return out
    mid = len(items) // 2
    merged = _convert_recursive(items[:mid], ph, depth + 1)
    merged.update(_convert_recursive(items[mid:], ph, depth + 1))
    return merged


def _convert_chunk(items: list[tuple[str, str]], ph: float) -> dict[str, str]:
    out: dict[str, str] = {}
    for start in range(0, len(items), _BATCH):
        chunk = items[start:start + _BATCH]
        with tempfile.NamedTemporaryFile("w", suffix=".smi", delete=False) as fh:
            for cid, smi in chunk:
                # The id is the TITLE field; obabel echoes it, and it is how the
                # result is matched back. Tabs and newlines in an id would break
                # that, so they are refused rather than silently mangled.
                if "\t" in cid or "\n" in cid:
                    raise ValueError(f"id contains a tab or newline: {cid!r}")
                fh.write(f"{smi}\t{cid}\n")
            path = Path(fh.name)
        try:
            proc = subprocess.run(
                [OBABEL, str(path), "-osmi", "-p", str(ph)],
                capture_output=True, text=True, timeout=1800)
            for line in proc.stdout.splitlines():
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                out[parts[1].strip()] = parts[0].strip()
        finally:
            path.unlink(missing_ok=True)
    return out


def charge_from_smiles(smiles: str) -> int | None:
    """Net formal charge of an already-protonated SMILES."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else Chem.GetFormalCharge(mol)


def charge_at_ph(smiles_by_id: dict[str, str], ph: float = DEFAULT_PH,
                 ) -> dict[str, int | None]:
    """{id: net charge at `ph`}, or None where obabel or RDKit could not say.

    None is not 0. A molecule whose charge could not be determined is a
    different fact from a neutral one, and collapsing them would put failures
    into the largest stratum.
    """
    prot = protonate(smiles_by_id, ph)
    return {cid: (charge_from_smiles(prot[cid]) if cid in prot else None)
            for cid in smiles_by_id}


def has_phosphate(smiles: str) -> bool | None:
    """Does this molecule carry a phosphate or phosphonate group?

    LABEL, NOT A FILTER (#6 item 5, decided). Pin1 BINDS phosphate -- its most
    potent non-covalent binders carry one -- and those same compounds are
    famously cell-impermeable, which is why Pfizer replaced theirs with a
    carboxylate. Both facts are true at once, so the pipeline records the
    property and lets a chemist weigh it. A rule rejecting phosphorus was
    proposed and DISCARDED because it rejected four known binders; see #12.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    patt = Chem.MolFromSmarts(PHOSPHATE_SMARTS)
    return bool(mol.GetSubstructMatches(patt))


def charge_class(charge: int | None) -> str:
    """The stratum label. Coarse on purpose.

    Vina carries no electrostatic term, so a charge stratum is not a claim
    about an energy -- it is a statement that comparing a dianion's score with
    a cation's is comparing two different physical situations. Three classes
    (anion / neutral / cation) is as fine as that argument supports; splitting
    -1 from -2 would imply a resolution the reasoning does not have.

    `unknown` is its own class rather than being folded into `neutral`.
    """
    if charge is None:
        return "unknown"
    if charge < 0:
        return "anion"
    if charge > 0:
        return "cation"
    return "neutral"
