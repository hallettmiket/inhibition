# The shared receptor

**PDB 6VAJ** — human Pin1 + sulfopin (ligand QT7), covalent at Cys113, 1.42 Å.

Verified first-hand rather than on the spec's assertion:

```
HEADER    ISOMERASE/ISOMERASE INHIBITOR           17-DEC-19   6VAJ
TITLE     CRYSTAL STRUCTURE ANALYSIS OF HUMAN PIN1
REMARK   2 RESOLUTION.    1.42 ANGSTROMS.
LINK         SG  CYS A 113                 C10 QT7 A 201     1555   1555  1.78
```

The `LINK` record confirms a genuine covalent bond at **1.78 Å** and hands over
the exact ligand attachment atom.

## Preparation

`python -m shared.receptor_prep` produces, under
`/data/lab_vm/immutable/inhibition/receptor/`:

| Artifact | Contents |
|---|---|
| `6VAJ_prepared.pdb` | protonated at pH 7.4, ligand and solvent stripped |
| `6VAJ_prepared.pdbqt` | rigid receptor for Vina/smina |
| `box.json` | 20 Å covalent box (T_3, T_4) |
| `box_expanded.json` | 26 Å expanded box (T_1, T_2) |
| `prep_log.json` | input/output SHA-256s, atom counts, box derivation |

Result on 6VAJ: **1,215 protein atoms kept**; 16 ligand and 173 solvent atoms
removed; **0 unrecognized heteroatoms retained**. Cys113 SG sits 4.26 Å from the
box centre, inside both boxes.

## Two boxes, not one

The reference ligand is covalent at Cys113, so a tight box around it is centred
on the **warhead sub-pocket**. Correct for T_3/T_4, which attack that atom;
wrong for T_1/T_2, which would be biased toward a sub-pocket they have no reason
to prefer. See [D0002](../decisions/index.md).

## Heteroatoms: retain and report

Anything not on the strip list is **kept and counted** as
`other_het_atoms_retained` in `prep_log.json`.

!!! danger "A non-zero count means go look"
    On the first run this retained 31 atoms of **PG4** (tetraethylene glycol, a
    cryoprotectant). It sat 22.65 Å from the box centre so no result was
    affected, but it had no business in a docking receptor. Stripping blindly
    would have removed it silently — and would remove a structural metal just
    as silently. See [D0003](../decisions/index.md) and the
    [receptor runbook](../runbooks/receptor_selection.md).
