"""
Purpose: A labelled sub-pocket must be the residues it claims to be.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: integration/app/pose3d.py + the prepared receptor
Output: pass/fail

WHY THIS EXISTS. Issue #3 asked for the receptor drawn as a surface with Pin1's
three sub-pockets labelled: the proline-binding pocket, the phosphate-binding
Arg loop and the Cys113 pocket. A residue number in that labelling is a claim,
and the cheapest way to get one wrong is to copy it from a paper that numbered a
different construct. A mislabelled pocket is WORSE than an unlabelled one: it is
confidently wrong, it looks authoritative, and a chemist reasons from it.

So the numbers are not trusted from the source they were read out of. Every
declared residue is checked against `6VAJ_prepared.pdb` itself, by name. If the
receptor is ever re-prepared, re-numbered or swapped, this fails here rather
than silently relabelling the pocket in front of a user.

The residue sets themselves come from the primary literature, not from the
structure alone — the structure can say what is near what, but not what the
region is conventionally called:
  Ranganathan et al., Cell 89:875 (1997), PMID 9200606 — basic triad
    Lys63/Arg68/Arg69; proline pocket Leu122/Met130/Phe134.
  Behrsin et al., J Mol Biol 365:1143 (2007), PMID 17113106 — corroborates both.
  Dubiella et al., Nat Chem Biol 17:954 (2021), PMID 33972797 — sulfopin/6VAJ,
    i.e. this project's own receptor.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "integration" / "app"))
import pose3d as p3d  # noqa: E402

# READABLE, NOT PRESENT. `is_file()` is True for a file this process cannot
# open -- the data root's ACL is invisible to the client and denies per user.
# With a presence check these skips did not fire, and thirteen tests reported
# assertion failures about Arg-loop residues and pocket size when the real
# cause was that the receptor could not be read at all.
RECEPTOR_READABLE = p3d.receptor_readable()
needs_receptor = pytest.mark.skipif(
    not RECEPTOR_READABLE,
    reason="prepared receptor is not readable by this user (absent, or denied "
           "by the data root's ACL)")


def test_the_three_sub_pockets_issue_3_asked_for_are_all_present():
    keys = {s.key for s in p3d.SUBPOCKETS}
    assert keys == {"proline", "phosphate", "cys113"}, (
        "issue #3 names exactly three sub-pockets: the proline-binding pocket, "
        "the phosphate-binding Arg loop and the Cys113 pocket")


@needs_receptor
def test_every_declared_residue_matches_the_actual_receptor():
    """The guard. A wrong residue number fails here, not in front of a user."""
    problems = p3d.verify_subpockets()
    assert not problems, (
        "sub-pocket residues disagree with 6VAJ_prepared.pdb:\n  "
        + "\n  ".join(problems))


@needs_receptor
def test_the_basic_triad_is_lys63_arg68_arg69():
    """Named explicitly in issue #3; verified against the file, not assumed."""
    res = p3d._receptor_residues()
    triad = p3d.SUBPOCKETS_BY_KEY["phosphate"]
    assert triad.resi == [63, 68, 69]
    assert [res[i][0] for i in triad.resi] == ["LYS", "ARG", "ARG"]


@needs_receptor
def test_the_catalytic_tetrad_is_cys113_his59_his157_ser154():
    res = p3d._receptor_residues()
    cys = p3d.SUBPOCKETS_BY_KEY["cys113"]
    assert set(cys.resi) == {59, 113, 154, 157}
    assert {res[i][0] for i in cys.resi} == {"HIS", "CYS", "SER"}


@needs_receptor
def test_the_proline_pocket_keeps_the_canonical_trio():
    """Leu122/Met130/Phe134 are the 1997/2007 definition and must not drift."""
    res = p3d._receptor_residues()
    pro = p3d.SUBPOCKETS_BY_KEY["proline"]
    assert {122, 130, 134} <= set(pro.resi)
    assert res[122][0] == "LEU" and res[130][0] == "MET" and res[134][0] == "PHE"


def test_the_regions_are_disjoint():
    """A residue in two coloured regions makes the surface unreadable.

    Dubiella's inhibitor-era proline-pocket description also sweeps in Thr152
    and His157, which the older definition assigns to the catalytic site. The
    overlap is resolved by claiming each residue once and SAYING so in `why`,
    not by colouring it twice.
    """
    seen: dict[int, str] = {}
    for sp in p3d.SUBPOCKETS:
        for resi in sp.resi:
            assert resi not in seen, (
                f"residue {resi} is claimed by both {seen[resi]!r} and "
                f"{sp.key!r}; 3Dmol draws the last surface over the first, so "
                "one of the two labels would be silently wrong")
            seen[resi] = sp.key


def test_each_region_has_its_own_colour():
    colours = [s.colour for s in p3d.SUBPOCKETS]
    assert len(set(colours)) == len(colours), (
        "two sub-pockets sharing a colour cannot be told apart on the surface")


def test_every_region_explains_itself():
    for sp in p3d.SUBPOCKETS:
        assert len(sp.why) > 60, (
            f"{sp.key} needs a real explanation — the label alone does not "
            "tell a reader what the region is for")


@needs_receptor
def test_each_centroid_is_near_its_own_residues():
    """A label placed at the wrong centroid points at the wrong pocket."""
    import math

    res = p3d._receptor_residues()
    for sp in p3d.SUBPOCKETS:
        c = p3d.subpocket_centroid(sp.key)
        assert c is not None, f"no centroid for {sp.key}"
        nearest = min(
            math.dist((x, y, z), c)
            for resi in sp.resi for _, x, y, z in res[resi][1])
        assert nearest < 6.0, (
            f"{sp.key}'s label sits {nearest:.1f} A from its nearest own atom")


@needs_receptor
def test_the_catalytic_sg_is_read_from_the_receptor_not_hardcoded():
    """It must agree with what receptor_prep recorded, to 3 decimal places."""
    import json

    sg = p3d.catalytic_sg()
    assert sg is not None
    log = json.loads(
        (p3d.RECEPTOR.parent / "prep_log.json").read_text(encoding="utf-8"))
    recorded = tuple(log["cys113_sg"])
    assert all(abs(a - b) < 1e-3 for a, b in zip(sg, recorded)), (
        f"pose3d reads SG at {sg}, prep_log recorded {recorded}")


@needs_receptor
def test_the_grey_surface_covers_the_measured_pocket_not_a_sequence_window():
    """`resi 101..125` is a sequence window and is not the binding site.

    The earlier surface used exactly that, which included residues pointing
    away from the site and excluded the entire Arg loop.
    """
    lining = set(p3d.pocket_resi())
    assert {63, 68, 69} <= lining, "the Arg loop must be on the surface"
    assert len(lining) > 25
    assert lining - set(range(101, 126)), (
        "the lining must not be a contiguous sequence window around Cys113")


def test_the_legend_names_every_region():
    legend = p3d.subpocket_legend()
    for sp in p3d.SUBPOCKETS:
        assert sp.label in legend and sp.colour in legend


def test_the_pymol_script_is_generated_from_the_same_definitions():
    """The export exists to be trusted MORE than the embedded viewer.

    If the .pml named its own residues, the session and the GUI could disagree
    about where the proline pocket is — and the user would be looking at the
    session.
    """
    pml = p3d.pymol_script("poses.sdf", covalent=True)
    for sp in p3d.SUBPOCKETS:
        assert f"select {sp.key}, receptor and resi " in pml
        for resi in sp.resi:
            assert str(resi) in pml
        assert f"color {sp.colour}, {sp.key}" in pml
    assert "show surface, receptor" in pml
    assert "set all_states, on" in pml, "every docked mode must load, not just one"
    assert str(p3d.CATALYTIC_RESI) in pml and p3d.CATALYTIC_ATOM in pml
