"""
Purpose: run the experimentally-determined covalent poses through our own criterion.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-09
Input: the covalent Pin1 crystal structures in 05_redock_benchmark/cases_1/pdb
Output: 00_outputs/blacksmith/crystal_controls/crystal_controls_<N>.csv + poses as SDF

@tt8804: *"use the xray crystal prepared poses (just cleave the bonds for non-cov),
I want to compare the controls run through our pipeline vs in the known xray
crystals."*

WHY THIS IS THE DECISIVE EXPERIMENT. #47 measured that the warhead classes with
crystal structures and measured kinetics score LAST on our near-attack criterion,
while a class with no measured Pin1 activity scores first. Two explanations fit
that equally well:

  (a) our DOCKING puts those molecules in the wrong place, so the criterion never
      sees the real geometry; or
  (b) our CRITERION is wrong, and would reject the real geometry too.

Only one experiment separates them: take the pose the crystallographer determined,
put it in our receptor frame, and score it. If the true geometry scores well, the
criterion is fine and docking is the defect. If the true geometry also scores
zero, the criterion is measuring the wrong thing and no amount of better docking
rescues it.

WHAT "CLEAVE THE BOND" MEANS HERE, AND THE TRAP IN IT. These ligands are bonded
to Cys113 SG -- the `LINK` record makes it explicit -- so the deposited geometry
is the PRODUCT, not the pre-reaction complex. The reactive carbon sits ~1.8 A from
the sulfur. Our near-attack window is 2.8-4.2 A. **A bonded crystal pose is
therefore too CLOSE to score, and would read as zero for a reason that has nothing
to do with whether the criterion works.**

So cleaving is necessary but not sufficient: the cleaved pose has to be allowed to
relax to a non-covalent equilibrium before the criterion means anything. That is
what the 10 ns sweep does, and it is why this script emits a pose for the sweep
rather than a score. The starting distance is recorded for every case so the
relaxation can be seen rather than assumed.

FRAME TRANSFER IS ON C-ALPHA, AND IS CHECKED. Each crystal is superposed onto the
production 3IKD receptor on shared backbone atoms; the RMSD of that fit is
reported per case. A pose is refused if the fit exceeds FIT_MAX_A, because a bad
superposition moves the ligand relative to the pocket and would be indistinguishable
from a bad pose.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                     # noqa: E402
from shared import run_paths as rp                     # noqa: E402

log = logging.getLogger("crystal-controls")
CASES = Path("/data/lab_vm/append_only/inhibition/05_redock_benchmark/cases_1/pdb")
RECEPTOR_3IKD = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb")
POSE_DIR = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/crystal_control_poses")

#: Backbone fit above this is refused. 1.5 A is generous for Pin1 across
#: crystal forms and still far below the scale that would move a ligand out of
#: the pocket.
FIT_MAX_A = 1.5

#: Residues that are never the ligand.
SKIP = {"HOH", "SO4", "PG4", "EDO", "GOL", "DMS", "ACT", "CL", "NA", "MG", "PO4",
        "TRS", "MPD", "PEG", "IOD", "ZN", "CA"}


def covalent_cases() -> list[dict]:
    """Every crystal whose LINK record bonds a ligand to Cys113 SG."""
    out = []
    for f in sorted(CASES.glob("*.pdb")):
        for ln in f.read_text(errors="replace").splitlines():
            if not ln.startswith("LINK"):
                continue
            if "SG  CYS" not in ln or "113" not in ln[:40]:
                continue
            # the partner residue name sits in the second half of the record
            comp = ln[47:50].strip()
            atom = ln[42:46].strip()
            if comp and comp not in SKIP:
                out.append({"pdb": f.stem, "path": f, "comp_id": comp,
                            "linked_atom": atom})
                break
    return out


def _ca(struct) -> dict:
    d = {}
    for ch in struct[0]:
        for r in ch:
            a = r.find_atom("CA", "*")
            if a:
                d[r.seqid.num] = np.array([a.pos.x, a.pos.y, a.pos.z])
    return d


def superpose(mobile, target):
    """Kabsch on shared C-alpha; returns (R, mobile_centroid, target_centroid, rmsd)."""
    m, t = _ca(mobile), _ca(target)
    common = sorted(set(m) & set(t))
    if len(common) < 40:
        raise ValueError(f"only {len(common)} shared CA atoms")
    P = np.array([m[i] for i in common])
    Q = np.array([t[i] for i in common])
    Pm, Qm = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - Pm).T @ (Q - Qm))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    rmsd = float(np.sqrt((((P - Pm) @ R.T - (Q - Qm)) ** 2).sum(1).mean()))
    return R, Pm, Qm, rmsd, len(common)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--out-list", default=None,
                    help="write the resulting idents for the sweep worklist")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    import gemmi
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, rdDetermineBonds
    RDLogger.DisableLog("rdApp.*")
    import pandas as pd

    cases = covalent_cases()
    log.info("%d crystals carry a covalent LINK to Cys113", len(cases))
    if not cases:
        raise SystemExit("no covalent cases found")

    target = gemmi.read_structure(str(RECEPTOR_3IKD))
    target.remove_alternative_conformations(); target.remove_hydrogens()
    sg = None
    for ch in target[0]:
        for r in ch:
            if r.seqid.num == 113:
                a = r.find_atom("SG", "*")
                if a:
                    sg = np.array([a.pos.x, a.pos.y, a.pos.z])
    if sg is None:
        raise SystemExit("no Cys113 SG in the production receptor")

    POSE_DIR.mkdir(parents=True, exist_ok=True)
    rows, idents = [], []
    for c in cases:
        rec = {"pdb": c["pdb"], "comp_id": c["comp_id"], "linked_atom": c["linked_atom"]}
        try:
            st = gemmi.read_structure(str(c["path"]))
            st.remove_alternative_conformations(); st.remove_hydrogens()
            R, Pm, Qm, fit, n_ca = superpose(st, target)
            rec.update(fit_rmsd_a=round(fit, 3), n_ca=n_ca)
            if fit > FIT_MAX_A:
                rec["status"] = f"refused: backbone fit {fit:.2f} A > {FIT_MAX_A}"
                rows.append(rec); continue

            # ONE COPY, THE ONE AT THIS CYSTEINE. 7F0M deposits the ligand in
            # several chains; pooling them gave 120 atoms and a 32.7 A "bond".
            copies = []
            for ch in st[0]:
                for r in ch:
                    if r.name != c["comp_id"]:
                        continue
                    cx, cn, ce = [], [], []
                    for a in r:
                        cx.append([a.pos.x, a.pos.y, a.pos.z])
                        cn.append(a.name); ce.append(a.element.name)
                    copies.append((np.array(cx), cn, ce))
            if not copies:
                rec["status"] = "ligand not found"; rows.append(rec); continue
            def _near(cp):
                q = (cp[0] - Pm) @ R.T + Qm
                return float(np.linalg.norm(q - sg, axis=1).min())
            xyz, names, elems = min(copies, key=_near)
            rec["n_copies"] = len(copies)
            xyz = (xyz - Pm) @ R.T + Qm                  # into the 3IKD frame

            # THE BOND IS CLEAVED BY OMISSION: the ligand is written on its own,
            # with no link to the cysteine. Nothing moves -- which is the point,
            # and why the starting distance below is the product geometry.
            d_sg = np.linalg.norm(xyz - sg, axis=1)
            i_link = names.index(c["linked_atom"]) if c["linked_atom"] in names else int(d_sg.argmin())
            rec.update(n_atoms=len(xyz),
                       linked_atom_to_sg_a=round(float(d_sg[i_link]), 3),
                       closest_atom_to_sg_a=round(float(d_sg.min()), 3))

            mol = Chem.RWMol()
            conf = Chem.Conformer(len(xyz))
            for k, (e, p) in enumerate(zip(elems, xyz)):
                mol.AddAtom(Chem.Atom(e.capitalize()))
                conf.SetAtomPosition(k, p.tolist())
            m = mol.GetMol(); m.AddConformer(conf)
            rdDetermineBonds.DetermineConnectivity(m)
            m.SetProp("_Name", f"xtal_{c['pdb']}_{c['comp_id']}")
            m.SetProp("pose_rank", "1"); m.SetProp("mode", "0")
            m.SetProp("source", f"crystal {c['pdb']}, covalent bond cleaved")
            ident = f"xtal_{c['pdb']}"
            w = Chem.SDWriter(str(POSE_DIR / f"{ident}.sdf")); w.write(m); w.close()
            rec["status"] = "ok"; rec["ident"] = ident
            idents.append(ident)
        except Exception as exc:                       # noqa: BLE001
            rec["status"] = f"failed: {type(exc).__name__}: {exc}"
        rows.append(rec)

    t = pd.DataFrame(rows)
    # Run-scoped: a control is evidence about ONE screen against ONE receptor.
    # Written flat, it was picked up by every later run's report rail.
    dest = sout.Topic("blacksmith", rp.controls_topic()).write("crystal_controls", ".csv")
    t.to_csv(dest, index=False)

    print("\n" + "=" * 76)
    print("  CRYSTAL CONTROLS — the experimentally determined pose, in our frame")
    print("=" * 76 + "\n")
    cols = [c for c in ("pdb", "comp_id", "status", "fit_rmsd_a", "n_atoms",
                        "linked_atom_to_sg_a", "closest_atom_to_sg_a") if c in t.columns]
    print(t[cols].to_string(index=False))
    okd = t[t.status == "ok"] if "status" in t else t
    if len(okd) and "linked_atom_to_sg_a" in okd:
        print(f"\n  The bonded carbon sits a median {okd.linked_atom_to_sg_a.median():.2f} A "
              f"from Cys113 SG.")
        print("  Our near-attack window is 2.8-4.2 A, so the PRODUCT geometry is too")
        print("  CLOSE to score. These poses must relax under MD before the criterion")
        print("  means anything -- that is what the sweep is for.")
    if args.out_list and idents:
        Path(args.out_list).write_text("\n".join(f"{i} 1" for i in idents) + "\n")
        print(f"\n  -> {args.out_list} ({len(idents)} poses)")
    print(f"  -> {dest}\n")


if __name__ == "__main__":
    main()
