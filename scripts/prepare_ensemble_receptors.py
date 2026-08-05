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
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import receptor_prep as rp              # noqa: E402

log = logging.getLogger("prepare-ensemble")

OUT = sout.Topic("blacksmith", "ensemble_receptors")
RCSB = "https://files.rcsb.org/download/{pdb}.pdb"
RCSB_CIF = "https://files.rcsb.org/download/{pdb}.cif"

# (pdb_id, reference ligand het code, why it is in the ensemble).
# The het codes come from #6 item 6 and are VERIFIED against the download
# before use -- see the module docstring.
MEMBERS = [
    ("3IKG", "J8Z", "cognate for the Potter-Astex seed"),
    ("3IKD", "J9Z", "cognate for the Du-Xu seed"),
    ("9INR", "A1D9K", "cognate for the Liu-2024-C3 seed"),
]

BOX_SIZE = 26.0          # matches box_expanded.json; T_1/T_2 need the room

# TWO SIMILAR-LOOKING TUPLES THAT MEAN DIFFERENT THINGS, so they are named
# rather than reused. `find_atom` takes an ATOM name; `assert_preserved` takes
# a RESIDUE name. Passing ("A", 113, "SG") to both -- which is the obvious
# thing to do, and what this file did first -- makes the post-condition
# compare an atom name against a residue name and fail with "catalytic SG113
# missing or renamed (found 'CYS')" on a receptor whose Cys113 SG is present
# and untouched. The guard was right; the constant was wrong.
CATALYTIC_ATOM = ("A", 113, "SG")     # chain, resid, ATOM name
CATALYTIC_RES = ("A", 113, "CYS")     # chain, resid, RESIDUE name


def download(pdb_id: str, dest: Path) -> Path:
    """Fetch the entry, falling back to mmCIF where no legacy PDB exists.

    NEWER ENTRIES HAVE NO .pdb FILE AT ALL. 9INR -- the cognate structure for
    the Liu-2024-C3 seed, and so the one this ensemble most needs -- returns
    **404** for `.pdb` and 200 for `.cif`. The legacy format cannot represent
    every deposition and RCSB simply does not publish one, so a hard-coded
    `.pdb` URL fails on exactly the recent structures an ensemble most wants,
    with an HTTP error that says nothing about why.

    The CIF is converted with obabel rather than parsed here, so everything
    downstream (strip -> reduce -> obabel -xr) is identical to the path 6VAJ,
    3IKG and 3IKD took. An ensemble exists to vary the RECEPTOR while holding
    everything else fixed; a second ingestion path would vary the ingestion.
    """
    if dest.is_file():
        log.info("%s already downloaded", pdb_id)
        return dest
    try:
        url = RCSB.format(pdb=pdb_id)
        log.info("downloading %s", url)
        with urllib.request.urlopen(url, timeout=60) as r:
            dest.write_bytes(r.read())
        return dest
    except urllib.error.HTTPError as exc:
        if exc.code != 404:
            raise
        log.warning("%s has no legacy .pdb (HTTP 404) — falling back to mmCIF",
                    pdb_id)

    cif = dest.with_suffix(".cif")
    log.info("downloading %s", RCSB_CIF.format(pdb=pdb_id))
    with urllib.request.urlopen(RCSB_CIF.format(pdb=pdb_id), timeout=180) as r:
        cif.write_bytes(r.read())
    proc = subprocess.run(["obabel", str(cif), "-O", str(dest)],
                          capture_output=True, text=True)
    if not dest.is_file() or dest.stat().st_size == 0:
        raise SystemExit(f"{pdb_id}: mmCIF -> PDB conversion produced nothing.\n"
                         + (proc.stderr or proc.stdout)[-600:])
    log.info("%s: converted mmCIF -> %s", pdb_id, dest.name)
    return dest



def ligand_coords_from_cif(cif: Path, het: str) -> list[tuple[float, float, float]]:
    """Heavy-atom coordinates of `het` read from the mmCIF `_atom_site` loop.

    NEEDED BECAUSE THE LEGACY PDB CANNOT HOLD THE LIGAND'S NAME. `A1D9K` is a
    5-character CCD code; the PDB format caps HET codes at 3, which is *why*
    9INR publishes no `.pdb`. obabel's CIF -> PDB conversion drops every
    HETATM, so the converted file has no ligand at all -- and asking it for
    `A1D9K` returns nothing, which reads as "wrong structure" rather than
    "this format cannot express the question".

    Columns are located BY NAME from the loop header, never by position:
    `_atom_site` column order is not fixed across depositions, and indexing by
    position would read the wrong column and return coordinates that are
    perfectly plausible numbers from the wrong field.
    """
    lines = cif.read_text(errors="replace").splitlines()
    cols: dict[str, int] = {}
    i, n = 0, len(lines)
    while i < n:
        if lines[i].strip() == "loop_":
            j, names = i + 1, []
            while j < n and lines[j].lstrip().startswith("_atom_site."):
                names.append(lines[j].strip().split(".", 1)[1])
                j += 1
            if names:
                cols = {name: k for k, name in enumerate(names)}
                i = j
                break
        i += 1
    if not cols:
        raise SystemExit(f"{cif.name}: no _atom_site loop found")

    def col(*candidates: str) -> int:
        for c in candidates:
            if c in cols:
                return cols[c]
        raise SystemExit(f"{cif.name}: none of {candidates} in _atom_site")

    c_comp = col("label_comp_id", "auth_comp_id")
    c_sym = col("type_symbol")
    cx, cy, cz = col("Cartn_x"), col("Cartn_y"), col("Cartn_z")

    out = []
    for line in lines[i:]:
        st = line.strip()
        if not st or st.startswith(("#", "loop_", "_")):
            break
        f = st.split()
        if len(f) <= max(c_comp, c_sym, cx, cy, cz):
            continue
        if f[c_comp] != het or f[c_sym].upper() == "H":
            continue
        try:
            out.append((float(f[cx]), float(f[cy]), float(f[cz])))
        except ValueError:
            continue
    return out


def prepare_one(pdb_id: str, het: str, note: str, *, work: Path) -> dict:
    raw = download(pdb_id, work / f"{pdb_id}.pdb")

    try:
        coords = rp.parse_ligand_coords(raw, het)
    except Exception:  # noqa: BLE001 - the PDB may not be able to name it
        coords = []
    if not coords:
        cif = raw.with_suffix(".cif")
        if cif.is_file():
            log.warning("%s: %r absent from the PDB (5-char CCD codes cannot "
                        "be written to it) — reading coordinates from the CIF",
                        pdb_id, het)
            coords = ligand_coords_from_cif(cif, het)
    if not coords:
        raise SystemExit(
            f"{pdb_id}: no atoms found for het code {het!r}. The box would be "
            "centred on nothing. Check the code against the entry before "
            "re-running -- a wrong centre docks into empty space and returns "
            "ordinary-looking affinities.")
    centre = rp.centroid(coords)
    log.info("%s: %s has %d atoms, centroid (%.2f, %.2f, %.2f)",
             pdb_id, het, len(coords), *centre)

    if rp.find_atom(raw, *CATALYTIC_ATOM) is None:
        raise SystemExit(
            f"{pdb_id}: no {CATALYTIC_ATOM[0]}:{CATALYTIC_ATOM[1]}:"
            f"{CATALYTIC_ATOM[2]} in the "
            "structure. Every ensemble member must carry the catalytic "
            "cysteine or it is not a receptor for this target.")

    stripped = work / f"{pdb_id}_stripped.pdb"
    counts = rp.strip_structure(raw, stripped, ligand_het=het)
    protonated = work / f"{pdb_id}_prepared.pdb"
    rp.protonate(stripped, protonated, 7.4)
    rp.assert_preserved(stripped, protonated, catalytic=CATALYTIC_RES,
                        label=pdb_id)

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
