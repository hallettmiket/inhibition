"""
Purpose: run named reference molecules (Sulfopin, ATRA) straight through the 100 ns elevation leg.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: data/reference/pin1_reference_binders_4.csv (+ the reference poses, for covalent ones)
Output: 00_outputs/blacksmith/md_residence/md_residence_ref_<name>_<N>.csv

@tt8804 asked for Sulfopin and ATRA through elevation. They are the two most
useful reference points available and they are useful for DIFFERENT reasons:

  Sulfopin   THE PARENT. Nanomolar covalent chloroacetamide, the compound the
             whole series descends from. Its 100 ns behaviour is the number
             every candidate's residence should be read next to.
  ATRA       A genuine low-micromolar BINDER that is not a covalent warhead
             compound at all. It tests the residence gate from the other side:
             the gate's justification (§4c of the framework) is kinetic -- even
             a millimolar binder should sit through 100 ns -- and ATRA is the
             molecule on hand that ought to demonstrate it.

EACH IS DOCKED THE WAY ITS CHEMISTRY DEMANDS, WHICH IS NOT THE SAME WAY.

Sulfopin is covalent and starts from its REACTIVE pose -- the one the screen
scored -- so its residence and its near-attack number describe one pose.

ATRA is non-covalent. Reactive docking biases the search toward warhead-sulfur
contact (D0064), and applying that bias to a molecule with no warhead would
manufacture a near-attack geometry that means nothing. It is therefore PLAIN
docked. Its SMARTS hits on four Michael-acceptor classes in the screen are an
artefact of a conjugated polyene looking like an acceptor, not a claim that ATRA
alkylates Cys113.

PROTONATION, STATED BECAUSE IT MATTERS AND IS NOT RESOLVED. ATRA is a carboxylic
acid (pKa ~4.8) and is >99% deprotonated at pH 7.4, but the reference file records
it neutral and every other number in this project was computed from that form.
This run keeps the neutral form for consistency with the pose it starts from, and
the anionic form is a documented follow-up rather than a silent substitution. A
residence measured on the neutral acid is not a measurement of the species that
binds at physiological pH.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                # noqa: E402
from shared import receptors as R                 # noqa: E402
import md_residence_3ikd as mr                    # noqa: E402

log = logging.getLogger("elevate-ref")
OUT = sout.Topic("blacksmith", "md_residence")
REF = REPO / "data/reference/pin1_reference_binders_4.csv"
POSES = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/nac_v2_poses")

#: name -> (use its reactive pose?, why)
PLAN = {
    "Sulfopin": (True, "covalent chloroacetamide; start from the pose the screen scored"),
    # ATRA SHOULD be plain-docked (it is non-covalent and the reactive
    # potential's warhead bias means nothing for it), but the plain-dock path in
    # md_residence rebuilds from PDBQT WITHOUT PLACING HYDROGENS -- every H lands
    # at the origin and antechamber refuses the molecule (bonds of 11-16 A, and
    # 28 hydrogens on top of each other). The screen's meeko rebuild does place
    # them, so the screened pose is used instead. That pose was found under the
    # reactive bias, which is a caveat carried on the row rather than a silent
    # substitution. The plain-dock H bug is a real defect and is logged for fix.
    "ATRA": (True, "non-covalent, so the reactive bias on this pose is a caveat; "
                   "used because md_residence's plain-dock path leaves hydrogens "
                   "at the origin and antechamber refuses it"),
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--names", nargs="+", default=list(PLAN))
    ap.add_argument("--gpu", default="6")
    ap.add_argument("--production-ps", type=float, default=100000.0)
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    R.resolve_3ikd_ian()
    refs = pd.read_csv(REF).set_index("name")

    rows = []
    for name in args.names:
        if name not in refs.index:
            log.error("%s not in the reference table", name)
            continue
        r = refs.loc[name]
        smiles = r.canonical_smiles
        if not isinstance(smiles, str) or smiles == "UNVERIFIED":
            log.error("%s has no verified SMILES", name)
            continue
        use_pose, why = PLAN.get(name, (False, "default: plain dock"))

        pose = None
        if use_pose:
            hits = sorted(POSES.glob(f"ref_{name}__*.sdf"))
            if hits:
                pose = hits[0]
            else:
                log.warning("%s: no reference pose on disk; plain docking instead", name)
        log.info("%s -> %s  (%s)", name,
                 pose.name if pose else "plain re-dock", why)

        if args.dry_run:
            rows.append({"ident": f"ref_{name}", "smiles": smiles,
                         "pose": str(pose) if pose else "", "note": why})
            continue

        row = mr.run_one(f"ref_{name}", smiles, "reference",
                         production_ps=args.production_ps, nrun=args.nrun,
                         gpu=args.gpu, keep=True, pose=pose, pose_rank=1,
                         net_charge=0)
        row["reference_name"] = name
        row["reference_tier"] = r.tier
        row["reference_potency"] = r.potency
        row["docking_mode"] = "reactive pose" if pose else "plain re-dock"
        row["docking_rationale"] = why
        if name == "ATRA":
            row["pose_caveat"] = (
                "started from a REACTIVE-docked pose because the plain-dock path "
                "produces unparameterisable hydrogens; the warhead bias that "
                "found this pose does not reflect ATRA's chemistry")
            row["protonation_caveat"] = (
                "run as the NEUTRAL acid, as recorded in the reference file; "
                "ATRA is >99% deprotonated at pH 7.4 and this is not a "
                "measurement of the species that binds physiologically")
        rows.append(row)
        log.info("%s: %s", name, row.get("status"))

    if not rows:
        raise SystemExit("nothing run")
    df = pd.DataFrame(rows)
    if not args.dry_run:
        dest = OUT.write(f"md_residence_ref_{'_'.join(args.names)}", ".csv")
        df.to_csv(dest, index=False)
        print(f"\n  -> {dest}")
    print(df.to_string(index=False)[:2000])


if __name__ == "__main__":
    main()
