#!/usr/bin/env python3
"""
Purpose: Render one docked pose as an interaction map -- contacting residues as
         sticks, polar contacts dashed and measured, centred on the ligand.
Author:  Timothy Wu (with Claude Code)
Date:    2026-08-25
Input:   a candidate id and a pose column
Output:  a self-contained HTML page (3Dmol.js) for screenshotting

WHY THIS IS NOT IN pose3d.py YET. `pose_html` draws a pose on a labelled
sub-pocket surface and offers three canned framings; it has no notion of which
residues a given pose actually touches. This adds that. It is written against
pose3d's own constants rather than beside them -- CATALYTIC_RESI, CATALYTIC_ATOM
and SUBPOCKETS are imported, never retyped -- because a residue label is a claim
about the structure and the cheapest way to get it wrong is to transcribe a
number.

THE NUMBERING GUARD. `docs/state_of_the_project.md` records a renderer that drew
a GLUTAMATE labelled as the target cysteine, because it bypassed the guard that
knows the MD system renumbers from 1 (PIN1_OFFSET = 50). This reads the
crystallographic receptor, where Cys113 really is residue 113 -- so it ASSERTS
that residue 113 is a CYS carrying an SG before drawing anything. If that fails,
the numbering is not what this code assumes and no figure should be produced.
"""
from __future__ import annotations

import collections
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "integration" / "app"))
import pose3d as p3d  # noqa: E402

#: A contact. Beyond this the residue is scenery, not interaction.
CONTACT_A = 4.5
#: Both partners N/O and within this: drawn as a dashed polar contact. Wider
#: than a textbook H-bond because the pose is docked, not refined, and a
#: hard 3.2 A cut would hide contacts the geometry plainly intends.
POLAR_A = 3.6
POLAR_ELEMENTS = {"N", "O"}


def receptor_residues(pdb: Path) -> dict:
    res = collections.defaultdict(list)
    for ln in pdb.read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        name = ln[12:16].strip()
        el = (ln[76:78].strip() or name[0]).upper()
        if el == "H":
            continue
        res[(ln[21], int(ln[22:26]), ln[17:20].strip())].append(
            (name, el, float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return res


def assert_numbering(res: dict) -> tuple:
    """Residue 113 must be a CYS with an SG, or this receptor is not what we think."""
    hits = [(k, v) for k, v in res.items() if k[1] == p3d.CATALYTIC_RESI]
    if not hits:
        raise SystemExit(f"no residue {p3d.CATALYTIC_RESI} in the receptor")
    (ch, num, rn), atoms = hits[0]
    sg = [a for a in atoms if a[0] == p3d.CATALYTIC_ATOM]
    if rn != "CYS" or not sg:
        raise SystemExit(
            f"residue {num} is {rn} with no {p3d.CATALYTIC_ATOM}; the numbering "
            "is not crystallographic and this renderer must not draw it")
    return sg[0]


def subpocket_of(num: int) -> tuple[str, str] | None:
    for sp in p3d.SUBPOCKETS:
        if num in sp.resi:
            return sp.label, sp.colour
    return None


def analyse(candidate: str, pose_column: str = "nac3_pose_path"):
    import data as D
    f, _ = D.load_frame(candidate.split("_")[0])
    row = f[f["candidate_id"] == candidate].iloc[0]
    pose = p3d.read_poses(Path(row[pose_column]))[0]
    rec = p3d.receptor_for(pose_column)
    res = receptor_residues(rec)
    sg = assert_numbering(res)

    lig = [(e.upper(), x, y, z)
           for e, x, y, z in p3d.pose_atoms(pose) if e.upper() != "H"]
    d_sg = min(math.dist(sg[2:], a[1:]) for a in lig)

    contacts = []
    for key, atoms in res.items():
        ch, num, rn = key
        best, polar = 1e9, None
        for an, ae, ax, ay, az in atoms:
            for le, lx, ly, lz in lig:
                d = math.dist((ax, ay, az), (lx, ly, lz))
                if d < best:
                    best = d
                if (ae in POLAR_ELEMENTS and le in POLAR_ELEMENTS
                        and d <= POLAR_A
                        and (polar is None or d < polar[0])):
                    polar = (d, (ax, ay, az), (lx, ly, lz))
        if best <= CONTACT_A:
            contacts.append({"chain": ch, "num": num, "resn": rn,
                             "min": best, "polar": polar,
                             "sub": subpocket_of(num)})
    contacts.sort(key=lambda c: c["min"])
    return {"row": row, "pose": pose, "receptor": rec, "sg": sg,
            "d_sg": d_sg, "contacts": contacts, "n_heavy": len(lig)}


def render(a: dict, *, width: int = 1200, height: int = 1000,
           zoom: float = 0.80, rotate_y: float = 0.0,
           rotate_x: float = 0.0, backdrop: bool = False,
           label_distances: bool = False,
           residue_sticks: bool = False) -> str:
    """Draw the pose with its contacts.

    `backdrop` adds a grey surface over the pocket residues that no sub-pocket
    claims. It is OFF by default: a FOURTH transparent surface makes headless
    swiftshader composite an all-white frame, every attempt, deterministically
    (three surfaces render; adding the grey one never does). Nothing in the
    scene is wrong -- the software GL simply will not draw it -- so the figure
    is built from three and the backdrop stays available for interactive use,
    where a real GPU handles it.
    """
    import py3Dmol
    v = py3Dmol.view(width=width, height=height)
    v.addModel(a["receptor"].read_text(), "pdb")
    v.setStyle({"model": 0}, {})

    lining = list(p3d.pocket_resi())
    claimed = {i for sp in p3d.SUBPOCKETS for i in sp.resi}
    rest = [i for i in lining if i not in claimed]
    # A THIN surface, not the opaque one pose_html draws. This panel is about
    # which residues touch the ligand, and an 0.78-opacity surface hides every
    # side chain that answers that.
    if backdrop and rest:
        v.addSurface("VDW", {"opacity": 0.22, "color": "lightgrey"},
                     {"model": 0, "resi": rest})
    # OPACITY FOLLOWS THE REPRESENTATION. With side chains drawn, the surface
    # has to be thin enough to see through; without them it is the only thing
    # representing the protein and a 0.26 wash reads as fog.
    # 0.85, NOT 1.0. A fully opaque surface makes headless swiftshader
    # composite an all-white frame every attempt; 0.85 renders and still reads
    # as solid. Measured, not chosen: 1.0 blank, 0.85 fine, 0.26 fine.
    opacity = 0.26 if residue_sticks else 0.85
    for sp in p3d.SUBPOCKETS:
        v.addSurface("VDW", {"opacity": opacity, "color": sp.colour},
                     {"model": 0, "resi": sp.resi})

    # Receptor side chains are OFF by default: the protein is the shape the
    # ligand sits in, and 18 sets of sticks bury that shape in line work. The
    # contacts are still carried -- by the labels, by the dashed polar
    # contacts, and by the sub-pocket colour under each one.
    if residue_sticks:
        for c in a["contacts"]:
            colour = c["sub"][1] if c["sub"] else "white"
            v.addStyle({"model": 0, "resi": c["num"], "chain": c["chain"]},
                       {"stick": {"radius": 0.13, "color": colour}})
        v.addStyle({"model": 0, "resi": p3d.CATALYTIC_RESI},
                   {"stick": {"colorscheme": "yellowCarbon", "radius": 0.20}})

    v.addModel(a["pose"].text, a["pose"].fmt)
    v.setStyle({"model": 1},
               {"stick": {"colorscheme": "cyanCarbon", "radius": 0.17}})

    # polar contacts, dashed and measured
    for c in a["contacts"]:
        if not c["polar"]:
            continue
        d, p_rec, p_lig = c["polar"]
        v.addCylinder({"start": {"x": p_rec[0], "y": p_rec[1], "z": p_rec[2]},
                       "end": {"x": p_lig[0], "y": p_lig[1], "z": p_lig[2]},
                       "radius": 0.045, "color": "magenta", "dashed": True,
                       "fromCap": 1, "toCap": 1})
        # The numeric distances are OFF by default. At figure size a dashed
        # line already says "polar contact", and five floating numbers collide
        # with the residue labels that identify what is contacting what --
        # "2.8" landed on top of "His59" and the pair read as "2.8 59". The
        # distances are reported in the panel's text instead, where they can be
        # read without a magnifier.
        if label_distances:
            mid = [(x + y) / 2 for x, y in zip(p_rec, p_lig)]
            v.addLabel(f"{d:.1f}",
                       {"position": {"x": mid[0], "y": mid[1], "z": mid[2]},
                        "backgroundColor": "white", "backgroundOpacity": 0.78,
                        "fontColor": "black", "fontSize": 10, "inFront": True,
                        "borderThickness": 0.3, "borderColor": "magenta"})

    # residue labels -- ONLY where the label carries information: a residue a
    # sub-pocket claims, a residue making a polar contact, or Cys113 itself.
    # Labelling all 18 contacts turns the panel into a word cloud, and the
    # unclaimed ones at 4 A are scenery.
    labelled = [c for c in a["contacts"]
                if c["sub"] or c["polar"] or c["num"] == p3d.CATALYTIC_RESI]
    res_all = receptor_residues(a["receptor"])
    for c in labelled:
        atoms = res_all[(c["chain"], c["num"], c["resn"])]
        cx = sum(x for _, _, x, _, _ in atoms) / len(atoms)
        cy = sum(y for _, _, _, y, _ in atoms) / len(atoms)
        cz = sum(z for _, _, _, _, z in atoms) / len(atoms)
        colour = c["sub"][1] if c["sub"] else "#f2f2f2"
        v.addLabel(f"{c['resn'].title()}{c['num']}",
                   {"position": {"x": cx, "y": cy, "z": cz},
                    "backgroundColor": colour, "backgroundOpacity": 0.85,
                    "fontColor": "black", "fontSize": 11, "inFront": True,
                    "borderThickness": 0.4, "borderColor": "#555555"})

    # the warhead-to-SG distance, whatever it is
    sg = a["sg"]
    lig_atoms = [(e.upper(), x, y, z)
                 for e, x, y, z in p3d.pose_atoms(a["pose"]) if e.upper() != "H"]
    near = min(lig_atoms, key=lambda t: math.dist(sg[2:], t[1:]))
    v.addCylinder({"start": {"x": sg[2], "y": sg[3], "z": sg[4]},
                   "end": {"x": near[1], "y": near[2], "z": near[3]},
                   "radius": 0.075, "color": "gold", "dashed": True,
                   "fromCap": 1, "toCap": 1})
    mid = [(sg[2] + near[1]) / 2, (sg[3] + near[2]) / 2, (sg[4] + near[3]) / 2]
    v.addLabel(f"Cys113 SG · {a['d_sg']:.2f} Å",
               {"position": {"x": mid[0], "y": mid[1], "z": mid[2]},
                "backgroundColor": "gold", "backgroundOpacity": 0.9,
                "fontColor": "black", "fontSize": 12, "inFront": True})

    # CENTRE ON THE LIGAND, then back off enough to hold its contacts.
    # pose_html's "pocket" framing zooms to every pocket-lining residue, which
    # puts the ligand off to one side of the frame.
    v.zoomTo({"model": 1})
    v.zoom(zoom)
    if rotate_y:
        v.rotate(rotate_y, "y")
    if rotate_x:
        v.rotate(rotate_x, "x")
    return v._make_html()


if __name__ == "__main__":
    cid = sys.argv[1]
    out = Path(sys.argv[2])
    zoom = float(sys.argv[3]) if len(sys.argv) > 3 else 0.80
    ry = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    rx = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    a = analyse(cid)
    out.write_text(render(a, zoom=zoom, rotate_y=ry, rotate_x=rx))
    print(f"{cid}  heavy={a['n_heavy']}  SG={a['d_sg']:.2f} A  "
          f"contacts={len(a['contacts'])}  polar="
          f"{sum(1 for c in a['contacts'] if c['polar'])}")
    for c in a["contacts"]:
        tag = f" polar {c['polar'][0]:.2f}" if c["polar"] else ""
        sub = c["sub"][0] if c["sub"] else "-"
        print(f"   {c['resn'].title()}{c['num']:<4} {c['min']:.2f} A  {sub}{tag}")
    print("WROTE", out)
