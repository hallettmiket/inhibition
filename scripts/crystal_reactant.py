"""
Purpose: rebuild the PRE-reaction molecule in its crystallographic binding pose.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-09
Input: the covalent crystal poses from crystal_controls.py + the free SMILES
Output: 00_outputs/blacksmith/crystal_control_poses/rx_<PDB>.sdf + a sidecar

WHY THE FIRST ATTEMPT COULD NOT WORK. Cleaving the bond to Cys113 leaves the
PRODUCT, and the product has no leaving group -- Sulfopin's deposited ligand is
16 heavy atoms against the free molecule's 17, because the chlorine left when the
bond formed. Our reactive SMARTS for a chloroacetamide is `C(=O)CCl`, so with no
chlorine there is no reactive atom to find:

    WARNING xtal_6VAJ: no reactive SMARTS match on the MD ligand

Both controls ran 10 ns cleanly and then had nothing to measure. That is not a
defect in the criterion; it is what a covalent crystal structure IS.

WHAT THIS BUILDS INSTEAD. The reactant: the intact molecule, in the pose the
crystallographer determined, with the leaving group put back where the chemistry
says it was.

THE ONE ASSUMPTION, STATED. The halogen's position is not in the deposition -- it
had already left. For an SN2 displacement the geometry is not free: the
nucleophile attacks anti to the leaving group, so in the product the departed
halide lay on the OPPOSITE side of the reactive carbon from the sulfur. The
halogen is therefore placed along the S->C vector, extended past C by a standard
C-X bond length. This is constructed, not measured, and the built distance is
reported per case so the assumption is visible rather than buried.

WHAT IT DOES NOT ASSUME. Nothing else moves. Every other atom keeps its
crystallographic coordinate, so the comparison against a docked pose is a
comparison of binding geometry and not of two different conformer generators.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                     # noqa: E402
from shared import run_paths as rp                     # noqa: E402

log = logging.getLogger("crystal-reactant")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
POSES = B / "crystal_control_poses"
SIDE = B / "pose_sidecars"
RECEPTOR = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb")

#: Standard bond lengths for the leaving group, Angstrom.
CX = {"Cl": 1.79, "Br": 1.94, "I": 2.14}

#: The molecules we hold trustworthy free SMILES for, and the deposited atom that
#: carried the bond to Cys113 (from each structure's own LINK record).
CASES = {
    "6VAJ": {"name": "Sulfopin", "linked": "C10",
             "free": "CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl"},
    "7F0M": {"name": "Liu-2022-ZL-Pin13", "linked": "C12",
             "free": "O=C(CCl)N1CCC2(CC1)SCC(=O)N2Cc1ccc(-c2cccc3ccccc23)o1"},
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--out-list", default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    import gemmi
    import pandas as pd
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")

    st = gemmi.read_structure(str(RECEPTOR)); st.remove_hydrogens()
    sg = None
    for ch in st[0]:
        for r in ch:
            if r.seqid.num == 113:
                a = r.find_atom("SG", "*")
                if a:
                    sg = np.array([a.pos.x, a.pos.y, a.pos.z])
    if sg is None:
        raise SystemExit("no Cys113 SG in the production receptor")

    SIDE.mkdir(parents=True, exist_ok=True)
    rows, idents = [], []
    for pdb, spec in CASES.items():
        rec = {"pdb": pdb, "name": spec["name"], "linked_atom": spec["linked"]}
        try:
            src = POSES / f"xtal_{pdb}.sdf"
            pose = next(m for m in Chem.SDMolSupplier(str(src), removeHs=False,
                                                      sanitize=False) if m)
            conf = pose.GetConformer()
            xyz = np.array([list(conf.GetAtomPosition(i))
                            for i in range(pose.GetNumAtoms())])

            free = Chem.MolFromSmiles(spec["free"])
            hal = [a for a in free.GetAtoms() if a.GetSymbol() in CX]
            if len(hal) != 1:
                raise ValueError(f"expected one halogen in the free SMILES, got {len(hal)}")
            sym = hal[0].GetSymbol()

            # the deposited carbon that held the bond to the sulfur
            i_c = int(np.argmin(np.linalg.norm(xyz - sg, axis=1)))
            d_sc = float(np.linalg.norm(xyz[i_c] - sg))
            v = xyz[i_c] - sg
            v /= np.linalg.norm(v)
            x_pos = xyz[i_c] + v * CX[sym]          # anti to the sulfur

            rw = Chem.RWMol(pose)
            new = rw.AddAtom(Chem.Atom(sym))
            rw.AddBond(i_c, new, Chem.BondType.SINGLE)
            c2 = Chem.Conformer(rw.GetNumAtoms())
            for i in range(pose.GetNumAtoms()):
                c2.SetAtomPosition(i, xyz[i].tolist())
            c2.SetAtomPosition(new, x_pos.tolist())
            built = rw.GetMol(); built.RemoveAllConformers(); built.AddConformer(c2)
            for a in built.GetAtoms():
                a.SetNoImplicit(False); a.SetNumExplicitHs(0)
            built.UpdatePropertyCache(strict=False)

            fixed = AllChem.AssignBondOrdersFromTemplate(free, built)
            if fixed.GetNumHeavyAtoms() != free.GetNumHeavyAtoms():
                raise ValueError(f"rebuilt {fixed.GetNumHeavyAtoms()} heavy atoms "
                                 f"against the free molecule's {free.GetNumHeavyAtoms()}")
            fixed.SetProp("_Name", f"rx_{pdb}")
            fixed.SetProp("pose_rank", "1"); fixed.SetProp("mode", "0")
            fixed.SetProp("source", f"crystal {pdb} reactant: deposited adduct with "
                                    f"{sym} rebuilt anti to Cys113 SG")
            ident = f"rx_{pdb}"
            w = Chem.SDWriter(str(POSES / f"{ident}.sdf")); w.write(fixed); w.close()
            (SIDE / f"{ident}.json").write_text(json.dumps(
                {"canonical_smiles": Chem.MolToSmiles(free), "charge_ph74": 0,
                 "source": f"crystal reactant {pdb} ({spec['name']}); leaving group "
                           "rebuilt anti to the sulfur"}))

            rec.update(status="ok", ident=ident, halogen=sym,
                       carbon_to_sg_a=round(d_sc, 3),
                       built_x_to_sg_a=round(float(np.linalg.norm(x_pos - sg)), 3),
                       n_heavy=fixed.GetNumHeavyAtoms())
            idents.append(ident)
        except Exception as exc:                       # noqa: BLE001
            rec["status"] = f"failed: {type(exc).__name__}: {exc}"
        rows.append(rec)

    t = pd.DataFrame(rows)
    # Run-scoped, for the same reason as crystal_controls: written flat, these
    # rows surfaced on the report rail of every later, unrelated run.
    dest = sout.Topic("blacksmith", rp.controls_topic()).write("crystal_reactant", ".csv")
    t.to_csv(dest, index=False)
    print("\n" + "=" * 74)
    print("  CRYSTAL REACTANTS — the intact molecule in its crystallographic pose")
    print("=" * 74 + "\n")
    print(t.to_string(index=False))
    print("\n  The halogen is CONSTRUCTED anti to the sulfur, not measured — it had")
    print("  already left when the structure was solved. Everything else keeps its")
    print("  crystallographic coordinate.")
    if args.out_list and idents:
        Path(args.out_list).write_text("\n".join(f"{i} 1" for i in idents) + "\n")
        print(f"\n  -> {args.out_list}")
    print(f"  -> {dest}\n")


if __name__ == "__main__":
    main()
