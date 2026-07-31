"""
Purpose: The pocket-derived size ceiling that generation is allowed to exceed by.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: none at runtime — the constants are derived and recorded here
Output: MAX_HEAVY_ATOMS and a predicate for pruning oversized candidates

WHERE THE NUMBER COMES FROM. Two independent estimates that agree, which is the
only reason to trust either.

STRUCTURAL. A grid cavity calculation on 6VAJ, restricted to the region within
6 A of the crystallographic ligand QT7 and requiring a point to be enclosed on
at least 4 of 6 directions, gives a pocket volume of **1018 A^3**. At the
typical 55% ligand packing fraction that is ~27 heavy atoms; at a tight 70%,
~35. (An earlier calculation over the whole 20 A docking box returned 2185 A^3
and 59-75 heavy atoms -- it was measuring surface grooves and open solvent, not
the pocket. Pin1's PPIase site is shallow and solvent-exposed, so a box-shaped
region centred on the ligand runs out of the pocket quickly.)

EMPIRICAL. Every non-peptidic Pin1 binder in our reference set:

    Juglone 13, Sulfopin 17, Pu-benzylguanine 21, ATRA 22, Tian 23,
    Guo-Pfizer 27, Du-Xu 27, Reddi-4d 27, Reddi-4g 29, Ieda 29, KPT-6566 30,
    PiB 32, Potter-Astex 32, EGCG 33, Liu-2024-C3 41, BJP-06-005-3 52

QT7, the 6VAJ co-crystal ligand, is 16. **Liu-2024-C3 at 41 has its own
co-crystal (9INR)**, so the pocket demonstrably accommodates 41 heavy atoms --
above the 70% packing estimate, which is what you expect for a ligand that
extends out of the core pocket into an adjoining groove.

THE CEILING IS 55, AND IT IS DELIBERATELY LOOSE. The brief was a generous limit
that removes molecules which can never fit, not a filter that shapes the
chemistry. 55 is ~1.6x the tight-packing estimate and clears every known
non-peptidic binder including BJP-06-005-3 at 52. What it excludes is the
peptidic macrocycles (65, 95, 169 heavy atoms), which are a different modality
that none of these approaches generates.

A FILTER MUST NOT REJECT A KNOWN BINDER. That is the test that sets the number,
and it is pinned in tests/test_pocket_size.py. Note T_1's existing
`max_heavy_atoms: 45` is TIGHTER than BJP-06-005-3 and would reject it.
"""

from __future__ import annotations

from rdkit import Chem

#: Pocket cavity volume, A^3. Grid calculation on 6VAJ within 6 A of QT7.
POCKET_VOLUME_A3 = 1018.0

#: Mean volume per heavy atom in an organic ligand, A^3.
VOLUME_PER_HEAVY_ATOM_A3 = 20.4

#: Ligand packing fractions used to bracket the estimate.
PACKING_TYPICAL = 0.55
PACKING_TIGHT = 0.70

#: What the structure implies, before any generosity is added.
HEAVY_ATOMS_TYPICAL_PACKING = int(
    POCKET_VOLUME_A3 * PACKING_TYPICAL / VOLUME_PER_HEAVY_ATOM_A3)   # 27
HEAVY_ATOMS_TIGHT_PACKING = int(
    POCKET_VOLUME_A3 * PACKING_TIGHT / VOLUME_PER_HEAVY_ATOM_A3)     # 34

#: The generation ceiling. See the module docstring for why it is this loose.
MAX_HEAVY_ATOMS = 55

#: The largest structurally-confirmed small-molecule binder (Liu-2024-C3, 9INR).
LARGEST_COCRYSTAL_BINDER_HEAVY_ATOMS = 41

#: The largest known non-peptidic binder of any kind (BJP-06-005-3).
LARGEST_KNOWN_NONPEPTIDIC_HEAVY_ATOMS = 52


def fits_pocket(smiles: str, max_heavy: int = MAX_HEAVY_ATOMS) -> bool:
    """Could this molecule conceivably occupy the Pin1 catalytic site?

    A deliberately weak claim. False means "too large to fit under any pose",
    which is safe to prune. True means only "not obviously too large" -- it is
    not evidence of binding and must never be read as such.
    """
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else smiles
    if m is None:
        return False
    return m.GetNumHeavyAtoms() <= max_heavy


def why(smiles: str, max_heavy: int = MAX_HEAVY_ATOMS) -> str | None:
    """A reason string when a molecule is pruned, or None when it is kept."""
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else smiles
    if m is None:
        return "unparseable"
    n = m.GetNumHeavyAtoms()
    if n > max_heavy:
        return (f"{n} heavy atoms > {max_heavy} ceiling "
                f"(pocket {POCKET_VOLUME_A3:.0f} A^3 admits "
                f"~{HEAVY_ATOMS_TIGHT_PACKING} at tight packing)")
    return None
