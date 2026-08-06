"""
Purpose: guards for the single-candidate workup — the ones that would have caught today's defects.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06

Every test here corresponds to a defect found while running the full simulation
workup on `t4_72f5671e89cb`, and each is written so that it FAILS on the code as
it stood this morning. A test that cannot fail on the broken version documents a
belief rather than checking a behaviour.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem                                    # noqa: E402
from rdkit.Chem import AllChem                            # noqa: E402


# --------------------------------------------------------------------------
# The junction frcmod covers the BDHI attachment carbon (GAFF2 c2)
# --------------------------------------------------------------------------

def test_frcmod_parser_reads_fixed_width_type_fields():
    """`S -c2` is one term, not the two whitespace tokens it looks like."""
    import covalent_workup_one as cw
    from shared import mmgbsa as mg

    terms = cw._frcmod_terms(mg.JUNCTION_FRCMOD)
    assert ("S", "c2") in terms["BOND"]
    assert ("2C", "S", "c2") in terms["ANGLE"]
    assert ("X", "c2", "S", "X") in terms["DIHE"]


def test_frcmod_parser_reports_an_absent_type_as_absent():
    """The check must be able to say no, or it is not a check.

    A looser parser (substring matching over the whole file) reported every
    term as present because some other line happened to contain its letters.
    """
    import covalent_workup_one as cw
    from shared import mmgbsa as mg

    terms = cw._frcmod_terms(mg.JUNCTION_FRCMOD)
    assert ("S", "cz") not in terms["BOND"]
    assert ("zz", "zz", "S") not in terms["ANGLE"]


def test_junction_covers_the_bdhi_sp2_attachment_carbon():
    """bdhi_c5 attaches through an sp2 C=N carbon, which GAFF2 types `c2`.

    D0067 corrected the GEOMETRY this warhead is scored with; this is the
    force-field half of the same fact. The terms needed once Cys113's SG
    replaces the cap hydrogen are S-c2, 2C-S-c2, the angles to the attachment
    carbon's heavy neighbours (c5 and n2 in this ring), and X-c2-S-X.
    """
    import covalent_workup_one as cw
    from shared import mmgbsa as mg

    terms = cw._frcmod_terms(mg.JUNCTION_FRCMOD)
    required = [("BOND", ("S", "c2")),
                ("ANGLE", ("2C", "S", "c2")),
                ("ANGLE", ("c5", "c2", "S")),
                ("ANGLE", ("n2", "c2", "S")),
                ("DIHE", ("X", "c2", "S", "X"))]
    missing = [(s, t) for s, t in required
               if t not in terms[s] and tuple(reversed(t)) not in terms[s]]
    assert not missing, f"{mg.JUNCTION_FRCMOD.name} is missing {missing}"


# --------------------------------------------------------------------------
# The receptor is checked, not inherited
# --------------------------------------------------------------------------

def _cys_pdb(tmp_path: Path, xyz, resnum: int = 63) -> Path:
    p = tmp_path / "receptor_cyx.pdb"
    x, y, z = xyz
    p.write_text(
        f"ATOM      1  SG  CYX A{resnum:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           S\n")
    return p


def test_receptor_identity_accepts_the_receptor_it_claims(tmp_path):
    import covalent_workup_one as cw
    out = cw.assert_receptor_identity(
        _cys_pdb(tmp_path, cw.RECEPTORS["3IKD"]["sg_xyz"]), "3IKD", 63)
    assert out["receptor_verified"] == "3IKD"
    assert out["sg_offset_a"] < cw.SG_TOLERANCE_A
    # The two receptors really are ~48.6 A apart — the number D0059 quotes.
    assert 48.0 < out["sg_offset_to_other_receptors_a"]["6VAJ"] < 49.5


def test_receptor_identity_rejects_the_default_6vaj_receptor(tmp_path):
    """The defect this guard exists for.

    `mmgbsa.prepare_receptor(workdir, receptor_pdb=None)` defaults to 6VAJ. A
    leg that forgets the override builds against a pocket 48.6 A away and
    returns a plausible dG. Claiming 3IKD while holding 6VAJ's coordinates
    must raise.
    """
    import covalent_workup_one as cw
    with pytest.raises(cw.WorkupError, match="6VAJ, not 3IKD"):
        cw.assert_receptor_identity(
            _cys_pdb(tmp_path, cw.RECEPTORS["6VAJ"]["sg_xyz"]), "3IKD", 63)


def test_receptor_identity_looks_up_the_renumbered_residue(tmp_path):
    """`prepare_receptor` renumbers, so Cys113 is not residue 113.

    Searching for the literal number 113 finds nothing in either prepared
    receptor, and a guard that cannot locate its own atom fails on correct
    input.
    """
    import covalent_workup_one as cw
    p = _cys_pdb(tmp_path, cw.RECEPTORS["3IKD"]["sg_xyz"], resnum=63)
    assert cw.assert_receptor_identity(p, "3IKD", 63)["sg_offset_a"] == 0.0
    with pytest.raises(cw.WorkupError, match="no CYS/CYX SG"):
        cw.assert_receptor_identity(p, "3IKD", 113)


# --------------------------------------------------------------------------
# A multi-pose SDF is selected by rank, never by position
# --------------------------------------------------------------------------

def _multipose_sdf(tmp_path: Path, ranks=(1, 2, 3), reverse=False) -> Path:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    p = tmp_path / "poses.sdf"
    w = Chem.SDWriter(str(p))
    for r in (reversed(ranks) if reverse else ranks):
        m = Chem.Mol(mol)
        m.SetProp("pose_rank", str(r))
        m.SetProp("nac_distance_A", f"{3.0 + r / 10:.3f}")
        w.write(m)
    w.close()
    return p


def test_read_pose_accepts_the_multi_pose_export(tmp_path):
    """The export writes the top N poses; the reader used to refuse all of them."""
    import bpmd_run
    p = _multipose_sdf(tmp_path)
    mol, props = bpmd_run.read_pose(p, pose_rank=2)
    assert props["pose_rank"] == 2


def test_read_pose_selects_by_rank_not_by_file_order(tmp_path):
    """Taking mols[0] passes today and silently biases the wrong pose later."""
    import bpmd_run
    p = _multipose_sdf(tmp_path, reverse=True)      # file order 3, 2, 1
    _, props = bpmd_run.read_pose(p, pose_rank=1)
    assert props["pose_rank"] == 1


def test_read_pose_refuses_a_rank_that_is_not_there(tmp_path):
    import bpmd_run
    p = _multipose_sdf(tmp_path)
    with pytest.raises(bpmd_run.BPMDRunError, match="pose_rank=9"):
        bpmd_run.read_pose(p, pose_rank=9)


# --------------------------------------------------------------------------
# BPMD's CV is bounded by a force, not by a stop flag GROMACS ignores
# --------------------------------------------------------------------------

def test_plumed_input_walls_the_cv_inside_the_metad_grid():
    """Every replica of the convergence run died indexing METAD off its grid.

    COMMITTOR fired at 681 ps and GROMACS ignored the stop flag, so the CV kept
    going to 1.97 nm against GRID_MAX=2.0. The wall is a force and cannot be
    ignored; the grid must have headroom above it.
    """
    from shared import bpmd
    txt = bpmd.plumed_input(warhead_idx=100, sg_idx=200)
    assert "UPPER_WALLS" in txt
    assert bpmd.WALL_NM > bpmd.UNBOUND_NM, \
        "a wall at or below the unbound threshold would bias the escape barrier"
    assert bpmd.GRID_MAX_NM > bpmd.WALL_NM, \
        "the wall must sit inside the grid, or it cannot prevent the crash"


def test_plumed_grid_is_fine_enough_for_the_hill_width():
    """Bins coarser than SIGMA blur the barrier the run exists to measure."""
    from shared import bpmd
    assert bpmd.GRID_SPACING_NM <= bpmd.HILL_SIGMA_NM / 2
    txt = bpmd.plumed_input(warhead_idx=1, sg_idx=2)
    nbin = int([l for l in txt.splitlines() if "GRID_BIN" in l][0].split("=")[1])
    span = bpmd.GRID_MAX_NM - bpmd.GRID_MIN_NM
    assert abs(span / nbin - bpmd.GRID_SPACING_NM) < 1e-9


def test_plumed_input_still_refuses_zero_based_indices():
    """The pre-existing guard must survive the rewrite."""
    from shared import bpmd
    with pytest.raises(bpmd.BPMDError):
        bpmd.plumed_input(warhead_idx=0, sg_idx=5)


# --------------------------------------------------------------------------
# MM-GBSA can type a receptor that kept its crystallographic water
# --------------------------------------------------------------------------

def test_mmgbsa_sources_a_water_force_field():
    """3IKD keeps 2 ordered waters and ff19SB+gaff2 cannot type them.

    Without this, tleap exits 31 and the failure reads as "no usable complex
    topology" — which points at the junction parameters, which were fine.
    """
    from shared import mmgbsa as mg
    assert "water" in mg.WATER_FF
    src = Path(mg.__file__).read_text()
    assert "source {WATER_FF}" in src, \
        "WATER_FF is defined but the tleap script does not source it"


# --------------------------------------------------------------------------
# The pose handed to antechamber is the whole molecule
# --------------------------------------------------------------------------

def test_saved_pose_keeps_its_hydrogens(tmp_path):
    """Stripping Hs from a sanitized SDF leaves noImplicit set on the heavies.

    `Chem.AddHs` then adds nothing, RDKit sanitizes the heavy-atom skeleton
    happily and reports the correct SMILES for it, and antechamber is the first
    thing to object ("Weird atomic valence (2)"). The formula check catches it
    where the heavy-atom count does not.
    """
    import md_residence_3ikd as mr

    smiles = "CCO"
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    mol.SetProp("pose_rank", "1")
    src = tmp_path / "pose.sdf"
    w = Chem.SDWriter(str(src))
    w.write(mol)
    w.close()

    out = mr.saved_pose_to_sdf(src, tmp_path, smiles, pose_rank=1)
    built = [m for m in Chem.SDMolSupplier(str(out), removeHs=False) if m][0]
    assert sum(1 for a in built.GetAtoms() if a.GetAtomicNum() == 1) == 6


def test_saved_pose_rejects_a_pose_that_is_not_the_candidate(tmp_path):
    import md_residence_3ikd as mr

    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    mol.SetProp("pose_rank", "1")
    src = tmp_path / "pose.sdf"
    w = Chem.SDWriter(str(src))
    w.write(mol)
    w.close()
    with pytest.raises(mr.ResidenceError):
        mr.saved_pose_to_sdf(src, tmp_path, "CCCCN", pose_rank=1)


# --------------------------------------------------------------------------
# The undefined warhead stereocentre is inherited, not re-drawn
# --------------------------------------------------------------------------

def test_the_candidates_warhead_stereocentre_is_undefined():
    """The premise of the whole stereo problem, asserted rather than assumed."""
    free = "O=S1(=O)CC[C@@H](N(Cc2cnoc2)C2CC(Br)=NO2)C1"
    centres = Chem.FindMolChiralCenters(
        Chem.MolFromSmiles(free), includeUnassigned=True,
        useLegacyImplementation=False)
    assert dict(centres)[13] == "?", \
        "if this centre becomes defined, the inheritance machinery is moot"


def test_adduct_inherits_the_poses_configuration(tmp_path):
    """Re-embedding the adduct picked the OPPOSITE diastereomer.

    Measured on this candidate: the NAC pose is (5R, 13R) and the adduct
    embedded from SMILES with randomSeed=42 is (5R, 13S). The covalent leg was
    scoring a different molecule from the non-covalent legs and nothing raised.
    """
    import covalent_workup_one as cw

    free = "O=S1(=O)CC[C@@H](N(Cc2cnoc2)[C@H]2CC(Br)=NO2)C1"   # the pose: 13R
    adduct = "O=S1(=O)CC[C@@H](N(Cc2cnoc2)C2CC=NO2)C1"          # unassigned
    mol = Chem.AddHs(Chem.MolFromSmiles(free))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    mol.SetProp("pose_rank", "1")
    src = tmp_path / "pose.sdf"
    w = Chem.SDWriter(str(src))
    w.write(mol)
    w.close()

    _, stereo = cw.adduct_from_pose(src, adduct, "Br", tmp_path / "add.sdf")
    assert dict((i, l) for i, l in stereo["pose_free_centres"])[13] == "R"
    assert dict((i, l) for i, l in stereo["adduct_centres"])[13] == "R"
    assert 13 in stereo["unassigned_in_smiles"]


def test_adduct_from_pose_refuses_a_different_molecule(tmp_path):
    """The pose minus its leaving group must BE the library's adduct."""
    import covalent_workup_one as cw

    mol = Chem.AddHs(Chem.MolFromSmiles("BrCCO"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    mol.SetProp("pose_rank", "1")
    src = tmp_path / "pose.sdf"
    w = Chem.SDWriter(str(src))
    w.write(mol)
    w.close()
    with pytest.raises(cw.WorkupError, match="not the same molecule"):
        cw.adduct_from_pose(src, "CCCCN", "Br", tmp_path / "add.sdf")
