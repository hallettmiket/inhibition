"""
Purpose: prepare 3IKG, 3IKD and 9INR for the receptor ensemble, each with its own box.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: RCSB (downloaded), shared/receptor_prep.py's exact preparation path
Output: 00_outputs/blacksmith/ensemble_receptors/<PDB>_prepared_<N>.pdbqt + box_<PDB>_<N>.json

#6 item 6 / D0052. 6VAJ is already prepared; these are the other three members
of the ensemble, and each one's box is derived from ITS OWN reference ligand.

WHY THE SAME PATH AS 6VAJ, NOT A NEW ONE. `receptor_prep`'s strip -> reduce
-BUILD -> obabel -xr sequence is reused function by function rather than
reimplemented. An ensemble exists to vary the RECEPTOR while holding everything
else fixed; a second preparation path would vary the preparation too, and the
resulting spread would be uninterpretable -- we could not say whether a score
difference came from the conformer or from how it was protonated.

WHY THE BOX IS PER RECEPTOR. A box is a set of coordinates in one structure's
frame. 6VAJ's box centred on QT7 means nothing in 3IKG's frame; using it would
dock into whatever happens to sit at those coordinates in the other crystal --
plausibly empty space beside the site -- and return affinities that look
entirely ordinary. Each box is centred on that entry's own cognate ligand, at
the same 26 A used by `box_expanded.json`, because T_1/T_2 place whole
molecules with no anchor at Cys113 and need the room.

WHY OUTPUTS GO TO append_only AND NOT NEXT TO 6VAJ. A prepared receptor is a
DERIVED artefact. `immutable/` is for original source data that no code may
modify; 6VAJ_prepared.pdbqt living there is a pre-existing inconsistency, not a
precedent to copy. These are versioned by `shared/outputs.py` like every other
derived artefact, so a re-preparation never overwrites one a frame was docked
against.

NOT RUN BLIND. The ligand het code for each entry is asserted against the
downloaded file before anything is prepared: a typo'd code would silently give
`parse_ligand_coords` nothing, `centroid` would raise or return garbage, and a
box centred on garbage is exactly the populated-and-plausible failure this
project keeps writing.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import receptor_prep as rp              # noqa: E402

log = logging.getLogger("prepare-ensemble")

OUT = sout.Topic("blacksmith", "ensemble_receptors")
RCSB = "https://files.rcsb.org/download/{pdb}.pdb"

# (pdb_id, reference ligand het code, why it is in the ensemble).
# The het codes come from #6 item 6 and are VERIFIED against the download
# before use -- see the module docstring.
MEMBERS = [
    ("3IKG", "J8Z", "cognate for the Potter-Astex seed"),
    ("3IKD", "J9Z", "cognate for the Du-Xu seed"),
    ("9INR", "A1D9K", "cognate for the Liu-2024-C3 seed"),
]

BOX_SIZE = 26.0          # matches box_expanded.json; T_1/T_2 need the room
CATALYTIC = ("A", 113, "SG")


def download(pdb_id: str, dest: Path) -> Path:
    if dest.is_file():
        log.info("%s already downloaded", pdb_id)
        return dest
    url = RCSB.format(pdb=pdb_id)
    log.info("downloading %s", url)
    with urllib.request.urlopen(url, timeout=60) as r:
        dest.write_bytes(r.read())
    return dest


def prepare_one(pdb_id: str, het: str, note: str, *, work: Path) -> dict:
    raw = download(pdb_id, work / f"{pdb_id}.pdb")

    coords = rp.parse_ligand_coords(raw, het)
    if not coords:
        raise SystemExit(
            f"{pdb_id}: no atoms found for het code {het!r}. The box would be "
            "centred on nothing. Check the code against the entry before "
            "re-running -- a wrong centre docks into empty space and returns "
            "ordinary-looking affinities.")
    centre = rp.centroid(coords)
    log.info("%s: %s has %d atoms, centroid (%.2f, %.2f, %.2f)",
             pdb_id, het, len(coords), *centre)

    if rp.find_atom(raw, *CATALYTIC) is None:
        raise SystemExit(
            f"{pdb_id}: no {CATALYTIC[0]}:{CATALYTIC[1]}:{CATALYTIC[2]} in the "
            "structure. Every ensemble member must carry the catalytic "
            "cysteine or it is not a receptor for this target.")

    stripped = work / f"{pdb_id}_stripped.pdb"
    counts = rp.strip_structure(raw, stripped, ligand_het=het)
    protonated = work / f"{pdb_id}_prepared.pdb"
    rp.protonate(stripped, protonated, 7.4)
    rp.assert_preserved(stripped, protonated, catalytic=CATALYTIC)

    pdbqt = OUT.write(f"{pdb_id}_prepared", ".pdbqt")
    rp.to_pdbqt(protonated, pdbqt)

    box = OUT.write(f"box_{pdb_id}", ".json")
    box.write_text(json.dumps({
        "pdb_id": pdb_id,
        "reference_ligand": het,
        "note": note,
        "center_x": centre[0], "center_y": centre[1], "center_z": centre[2],
        "size_x": BOX_SIZE, "size_y": BOX_SIZE, "size_z": BOX_SIZE,
        "derived_from": "centroid of this entry's own cognate ligand",
        "matches": "box_expanded.json size (26 A), used_by [t1, t2]",
    }, indent=2) + "\n", encoding="utf-8")

    log.info("%s -> %s, %s", pdb_id, pdbqt.name, box.name)
    return {"pdb_id": pdb_id, "het": het, "pdbqt": str(pdbqt),
            "box": str(box), "centre": centre, **counts}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--only", nargs="*", help="prepare only these PDB ids")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    work = Path("/data/lab_vm/modifiable/inhibition/ensemble_receptor_prep")
    work.mkdir(parents=True, exist_ok=True)

    todo = [m for m in MEMBERS
            if not args.only or m[0] in {s.upper() for s in args.only}]
    results = [prepare_one(*m, work=work) for m in todo]

    summary = OUT.write("ensemble_receptor_prep", ".json")
    summary.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"\nprepared {len(results)} receptor(s) -> {OUT.dir}")
    for r in results:
        print(f"  {r['pdb_id']}  ligand {r['het']}  "
              f"centre ({r['centre'][0]:.2f}, {r['centre'][1]:.2f}, "
              f"{r['centre'][2]:.2f})")
    print(f"\nsummary: {summary}")


if __name__ == "__main__":
    main()
