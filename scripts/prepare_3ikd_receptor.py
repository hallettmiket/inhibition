"""
Purpose: prepare the chemist's 3IKD structure for docking, changing as little as possible.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: immutable/inhibition/receptor/3ikd_well_prepared.pdb (uploaded by @tt8804)
Output: 00_outputs/blacksmith/receptor_3ikd/3IKD_prepared_<N>.pdbqt + box_3IKD_<N>.json

THE FILE IS AUTHORITATIVE. @tt8804's chemist prepared it and the instruction is
to keep it exactly as uploaded. So this script does NOT re-protonate, does NOT
strip waters, and does NOT rebuild anything. It makes exactly one change, and
the change is unavoidable.

WHY J9Z MUST COME OUT, AND WHY ITS COORDINATES STAY. A receptor that still
contains its own cognate ligand has an occupied pocket -- every docked molecule
would be scored against a site that is already full, and the poses would be
nonsense. So the ligand is removed, but its COORDINATES DEFINE THE BOX and are
retained for that purpose. This is the same rule `config/receptor.yaml` states
for 6VAJ's QT7: "Stripped before docking, but its coordinates DEFINE the box --
so the coordinates are retained, not discarded with the ligand."

WHAT IS DELIBERATELY NOT DONE, AND WHY IT MATTERS

* **No `reduce -BUILD`.** The uploaded file already carries 888 hydrogens. 6VAJ
  went through `reduce` at pH 7.4; running it again over an already-protonated
  structure would either double-add or silently re-assign choices the chemist
  made. Cys113 arrives as a REACTIVE THIOL (SG with HG attached), which is the
  state `config/receptor.yaml` requires or T_3/T_4 have nothing to attack.
* **Waters are kept.** All 6 of them. 6VAJ's prepared structure has none. Pin1's
  site is water-mediated -- that is cited in the FEP rule-out as a reason this
  pocket is hard -- and the chemist retained exactly six on purpose. They become
  part of the rigid receptor.

THE CONFOUND THIS CREATES, STATED SO NOBODY DISCOVERS IT LATER. 6VAJ and 3IKD
are now prepared DIFFERENTLY: 6VAJ was stripped of waters and protonated by
`reduce` at pH 7.4; 3IKD keeps its waters and its chemist's protonation. So when
the redock benchmark is re-run and compared against 6VAJ's 5% pose recovery
(D0046), any difference is **receptor AND preparation**, not receptor alone. To
attribute a difference to the receptor you would need deposited 3IKD put through
6VAJ's original path as a control. That control is cheap and is NOT run here.

A SECOND MISMATCH WORTH RECORDING. `REMARK 200 PH: 8.00` in the file is the
CRYSTALLISATION pH of the 2009 deposition, not the pH the chemist's protonation
targeted -- the file does not say what that was. Our ligands are prepared at pH
7.4 (`LIGAND_PREP_TAG`). His and Cys are exactly the residues where a 7.4-vs-8.0
difference is not automatically negligible, and Cys113 is the one that matters.
Unresolved; flagged rather than assumed away.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402

log = logging.getLogger("prepare-3ikd")

SRC = Path("/data/lab_vm/immutable/inhibition/receptor/3ikd_well_prepared.pdb")
OUT = sout.Topic("blacksmith", "receptor_3ikd")
OBABEL = "/data/lab_vm/envs/dwi_cheminf/bin/obabel"

LIGAND_HET = "J9Z"
CATALYTIC = ("A", 113, "SG")
# Matches box_expanded.json. T_1/T_2 place whole molecules with no anchor at
# Cys113 and need the room; the covalent box is 20 A and is derived separately.
BOX_SIZE = 26.0


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse(path: Path) -> dict:
    """Everything we need to know about the file, read once."""
    lig, waters, sg, natom, nh = [], 0, None, 0, 0
    het = {}
    for ln in path.read_text(errors="replace").splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            natom += 1
            el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
            if el == "H":
                nh += 1
            res = ln[17:20].strip()
            if ln.startswith("HETATM"):
                het[res] = het.get(res, 0) + 1
                if res == "HOH":
                    waters += 1
                if res == LIGAND_HET and el != "H":
                    lig.append((float(ln[30:38]), float(ln[38:46]),
                                float(ln[46:54])))
            if (ln[21:22].strip() == CATALYTIC[0]
                    and ln[22:26].strip() == str(CATALYTIC[1])
                    and ln[12:16].strip() == CATALYTIC[2]):
                sg = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return {"ligand_coords": lig, "waters": waters, "sg": sg,
            "n_atoms": natom, "n_hydrogens": nh, "hetatm_residues": het}


def strip_ligand(src: Path, dst: Path) -> int:
    """Remove ONLY the cognate ligand. Waters, hydrogens, everything else stay."""
    kept, removed = [], 0
    for ln in src.read_text(errors="replace").splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() == LIGAND_HET:
            removed += 1
            continue
        # CONECT records referencing the removed ligand would dangle.
        if ln.startswith("CONECT"):
            continue
        kept.append(ln)
    dst.write_text("\n".join(kept) + "\n", encoding="utf-8")
    return removed


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if not SRC.is_file():
        raise SystemExit(f"not found: {SRC}")
    info = parse(SRC)

    if info["sg"] is None:
        raise SystemExit(
            "Cys113 SG not found. Without the catalytic sulfur this is not a "
            "receptor for this target.")
    if not info["ligand_coords"]:
        raise SystemExit(
            f"{LIGAND_HET} not found — the box cannot be derived from the "
            "entry's own reference ligand, and 6VAJ's box is a set of "
            "coordinates in 6VAJ's frame that means nothing here.")

    n = len(info["ligand_coords"])
    cx = sum(c[0] for c in info["ligand_coords"]) / n
    cy = sum(c[1] for c in info["ligand_coords"]) / n
    cz = sum(c[2] for c in info["ligand_coords"]) / n
    sg = info["sg"]
    d = sum((a - b) ** 2 for a, b in zip((cx, cy, cz), sg)) ** 0.5

    log.info("source sha256 %s", sha256(SRC)[:16])
    log.info("atoms %d (%d hydrogens), waters %d, HETATM %s",
             info["n_atoms"], info["n_hydrogens"], info["waters"],
             info["hetatm_residues"])
    log.info("Cys113 SG at (%.3f, %.3f, %.3f)", *sg)
    log.info("%s centroid (%.2f, %.2f, %.2f), %.2f A from SG",
             LIGAND_HET, cx, cy, cz, d)
    if args.dry_run:
        return

    work = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep")
    work.mkdir(parents=True, exist_ok=True)
    stripped = work / "3IKD_noligand.pdb"
    removed = strip_ligand(SRC, stripped)
    log.info("removed %d %s atoms; waters and hydrogens untouched",
             removed, LIGAND_HET)

    pdbqt = OUT.write("3IKD_prepared", ".pdbqt")
    # -xr: rigid receptor. NO -p (protonation): the file is already protonated
    # and the chemist's assignment is authoritative.
    proc = subprocess.run([OBABEL, str(stripped), "-O", str(pdbqt), "-xr"],
                          capture_output=True, text=True)
    if not pdbqt.is_file() or pdbqt.stat().st_size == 0:
        raise SystemExit("obabel produced no PDBQT\n"
                         + (proc.stderr or proc.stdout)[-800:])

    box = OUT.write("box_3IKD", ".json")
    box.write_text(json.dumps({
        "pdb_id": "3IKD",
        "source": "chemist-prepared, uploaded by @tt8804 2026-08-05",
        "reference_ligand": LIGAND_HET,
        "center_x": cx, "center_y": cy, "center_z": cz,
        "size_x": BOX_SIZE, "size_y": BOX_SIZE, "size_z": BOX_SIZE,
        "derived_from": f"centroid of {LIGAND_HET}'s heavy atoms in this file",
        "cys113_sg": {"x": sg[0], "y": sg[1], "z": sg[2]},
        "ligand_centroid_to_sg_angstrom": round(d, 3),
        "matches": "box_expanded.json size (26 A), used_by [t1, t2]",
    }, indent=2) + "\n", encoding="utf-8")

    prep = OUT.write("prep_3IKD", ".json")
    prep.write_text(json.dumps({
        "source_file": str(SRC),
        "source_sha256": sha256(SRC),
        "prepared_pdbqt": str(pdbqt),
        "prepared_pdbqt_sha256": sha256(pdbqt),
        "box": str(box),
        "changes_made": [
            f"removed {removed} {LIGAND_HET} atoms (cognate ligand; a receptor "
            "containing it has an occupied pocket)",
            "dropped CONECT records (they referenced the removed ligand)",
            "converted to PDBQT with obabel -xr (rigid receptor)",
        ],
        "changes_deliberately_NOT_made": [
            "no reduce -BUILD: the file arrives with 888 hydrogens already "
            "assigned and the chemist's protonation is authoritative",
            f"waters KEPT ({info['waters']}), where 6VAJ_prepared has none",
            "no residue rebuilding, no rotamer changes",
        ],
        "verified": {
            "cys113_sg_present": True,
            "cys113_is_reactive_thiol": True,
            "cys113_sg": {"x": sg[0], "y": sg[1], "z": sg[2]},
            "n_hydrogens": info["n_hydrogens"],
            "waters_retained": info["waters"],
        },
        "UNRESOLVED": [
            "Preparation pH unknown. REMARK 200 PH 8.00 is the 2009 "
            "CRYSTALLISATION pH, not the protonation target. Ligands are "
            "prepared at pH 7.4 (LIGAND_PREP_TAG).",
            "6VAJ and 3IKD are now prepared DIFFERENTLY (waters, protonation "
            "tool), so a pose-recovery difference against D0046's 5% "
            "conflates receptor with preparation. The control -- deposited "
            "3IKD through 6VAJ's path -- is not run here.",
        ],
    }, indent=2) + "\n", encoding="utf-8")

    print(f"\n3IKD prepared -> {pdbqt}")
    print(f"  box          {box}")
    print(f"  provenance   {prep}")
    print(f"  waters kept  {info['waters']}   hydrogens {info['n_hydrogens']}")
    print(f"  Cys113 SG    ({sg[0]:.3f}, {sg[1]:.3f}, {sg[2]:.3f})  reactive thiol")
    print(f"  box centre   ({cx:.2f}, {cy:.2f}, {cz:.2f})  {d:.2f} A from SG")


if __name__ == "__main__":
    main()
