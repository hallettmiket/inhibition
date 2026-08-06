"""
Purpose: the BPMD driver's index resolution and convergence arithmetic, on cases with known answers.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06

`tests/test_bpmd.py` covers the scoring layer. This covers the DRIVER, and it
concentrates on the two places where a wrong answer would still look right.

1. LOCATING THE WARHEAD IN THE SOLVATED SYSTEM. `shared/bpmd.py` refuses a
   0-based index, but nothing can refuse the *wrong atom's* index. The driver
   maps the pose's warhead onto the built system by internal geometry, so the
   tests here rotate and translate the molecule, permute its atom order, and
   check the same atom comes back. A test that only ran the identity case would
   pass against an implementation that ignored the geometry entirely.

2. THE FIRST SMOKE RUN FAILED HERE, FOR REAL. `combine {rec LIG}` implies the
   ligand sits at residue `n_residues`, and it does not: `solvatebox` moves every
   water to the end, including the receptor's own crystallographic ones, so the
   index landed on a WAT. On 3IKD that raised. Had the two crystallographic
   waters been ordered differently it would have landed on a protein residue and
   produced a plausible number, so the regression is worth pinning.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import bpmd_run as br                              # noqa: E402
from shared import bpmd                            # noqa: E402


# --------------------------------------------------------------------------
# mapping the pose's atoms onto the built system
# --------------------------------------------------------------------------

def _molecule():
    """Six heavy atoms with no internal symmetry, so every atom is distinguishable."""
    xyz = np.array([[0.0, 0.0, 0.0],
                    [1.5, 0.0, 0.0],
                    [2.4, 1.2, 0.0],
                    [1.9, 2.7, 0.3],
                    [0.2, 3.4, 1.1],
                    [-1.1, 1.3, 0.7]])
    z = [6, 6, 6, 8, 7, 17]
    return xyz, z


def _rigid(xyz, angle=0.7, shift=(11.0, -4.0, 3.0)):
    c, s = np.cos(angle), np.sin(angle)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return xyz @ rot.T + np.array(shift)


def test_identity_order_is_recognised_and_reported_as_verified():
    xyz, z = _molecule()
    out = br._locate_ligand_atom(xyz, z, 2, xyz.copy(), list(z), np.array([9.0, 9.0, 9.0]))
    assert out["sys_heavy_idx"] == 2
    assert "order preserved" in out["atom_map"]


def test_a_rigid_motion_does_not_move_the_warhead():
    """solvatebox places the solute in a box; the CV atom must survive that."""
    xyz, z = _molecule()
    out = br._locate_ligand_atom(xyz, z, 5, _rigid(xyz), list(z),
                                 np.array([0.0, 0.0, 0.0]))
    assert out["sys_heavy_idx"] == 5


def test_a_permuted_atom_order_still_finds_the_same_atom():
    """THE TEST THAT MATTERS. antechamber and tleap may renumber, and taking the
    pose's index straight across would then bias a different atom in a run that
    completes normally."""
    xyz, z = _molecule()
    perm = [3, 0, 5, 2, 4, 1]
    sys_xyz = _rigid(xyz)[perm]
    sys_z = [z[i] for i in perm]
    for target in range(len(z)):
        out = br._locate_ligand_atom(xyz, z, target, sys_xyz, sys_z,
                                     np.array([0.0, 0.0, 0.0]))
        assert perm[out["sys_heavy_idx"]] == target, f"atom {target} mis-mapped"
        assert "internal geometry" in out["atom_map"]


def test_a_different_molecule_is_refused_rather_than_matched():
    xyz, z = _molecule()
    other = xyz.copy()
    other[4] += 2.5                      # same atom count, different molecule
    with pytest.raises(br.BPMDRunError, match="not this pose"):
        br._locate_ligand_atom(xyz, z, 4, other[[1, 0, 2, 3, 4, 5]],
                               [z[i] for i in [1, 0, 2, 3, 4, 5]],
                               np.array([0.0, 0.0, 0.0]))


def test_a_changed_atom_count_is_refused():
    xyz, z = _molecule()
    with pytest.raises(br.BPMDRunError, match="heavy-atom count changed"):
        br._locate_ligand_atom(xyz, z, 0, xyz[:4], z[:4], np.array([0.0, 0.0, 0.0]))


def test_a_symmetric_molecule_resolves_toward_the_sulfur_and_says_so():
    """A genuine internal symmetry is a tie no rigid transform can break.

    Four corners of a rectangle all carry the same distance profile, so the
    fingerprint cannot name one. They are chemically equivalent, so the one
    facing Cys113 is taken — the same rule the SMARTS-level choice uses — and the
    tie is REPORTED rather than absorbed.
    """
    xyz = np.array([[-2.0, -1.0, 0.0], [2.0, -1.0, 0.0],
                    [2.0, 1.0, 0.0], [-2.0, 1.0, 0.0]])
    z = [6, 6, 6, 6]
    perm = [1, 0, 2, 3]                     # not one of the rectangle's symmetries
    sys_xyz = xyz[perm]
    sg = np.array([2.0, 1.0, 5.0])          # directly above the corner at (2, 1)
    out = br._locate_ligand_atom(xyz, z, 0, sys_xyz, [z[i] for i in perm], sg)
    assert out["map_ties"] == 4, "every corner is equivalent; the tie must surface"
    chosen = sys_xyz[out["sys_heavy_idx"]]
    assert chosen == pytest.approx(np.array([2.0, 1.0, 0.0])), \
        "the tie must break toward the sulfur"


def test_a_close_but_distinguishable_atom_is_not_swept_into_a_tie():
    """The tie rule must fire on SYMMETRY, not on similarity.

    A threshold loose enough to catch merely-similar atoms hands the
    nearest-the-sulfur tie-break two atoms that are simply different, and it then
    picks one — which is how the first version of this mapped a warhead onto the
    wrong carbon while reporting a clean match.
    """
    xyz, z = _molecule()
    perm = [3, 0, 5, 2, 4, 1]
    sys_xyz = _rigid(xyz)[perm]
    out = br._locate_ligand_atom(xyz, z, 2, sys_xyz, [z[i] for i in perm],
                                 np.array([0.0, 0.0, 0.0]))
    assert out["map_ties"] == 1
    assert perm[out["sys_heavy_idx"]] == 2


# --------------------------------------------------------------------------
# the residue-layout regression
# --------------------------------------------------------------------------

def test_pdb_residue_names_counts_residues_not_atoms(tmp_path):
    pdb = tmp_path / "receptor_cys.pdb"
    pdb.write_text(
        "ATOM      1  N   GLU A   1      0.000   0.000   0.000\n"
        "ATOM      2  CA  GLU A   1      1.000   0.000   0.000\n"
        "ATOM      3  N   CYS A   2      2.000   0.000   0.000\n"
        "TER\n"
        "HETATM    4  O   HOH A   3      3.000   0.000   0.000\n"
        "END\n")
    assert br._pdb_residue_names(pdb) == ["GLU", "CYS", "HOH"]


def test_sg_in_pose_frame_reads_the_renumbered_cysteine(tmp_path):
    """The pose-frame sulfur, taken by the index prepare_receptor renumbered to.

    Pin1 has two cysteines, so matching on "a CYS with an SG" is not enough — the
    wrong one would give a well-formed coordinate and a wrong reactive centre.
    """
    pdb = tmp_path / "receptor_cys.pdb"
    pdb.write_text(
        "ATOM      1  SG  CYS A   7      1.000   2.000   3.000\n"
        "ATOM      2  SG  CYS A  63      9.000   8.000   7.000\n"
        "END\n")
    assert br.sg_in_pose_frame(pdb, 63) == pytest.approx(np.array([9.0, 8.0, 7.0]))
    assert br.sg_in_pose_frame(pdb, 7) == pytest.approx(np.array([1.0, 2.0, 3.0]))
    with pytest.raises(br.BPMDRunError, match="no CYS99"):
        br.sg_in_pose_frame(pdb, 99)


def test_the_ligand_is_not_at_the_receptor_residue_count():
    """The bug the first smoke run caught, stated as arithmetic.

    receptor_cys.pdb for 3IKD holds 113 protein residues + 2 crystallographic
    waters = 115, so `combine {rec LIG}` reads as putting the ligand at 0-based
    index 115. tleap moves the waters to the end, so it is actually at 113. The
    protein-only conversion the driver does is what gets both the ligand and
    Cys113 right.
    """
    pdb_names = ["ALA"] * 113 + ["HOH", "HOH"]
    n_res = len(pdb_names)
    protein_only = sum(1 for n in pdb_names if n in br.PROTEIN_NAMES)
    assert n_res == 115 and protein_only == 113
    assert protein_only != n_res, "indexing the ligand at n_residues hits a water"


def test_cys_index_converts_into_the_protein_only_sequence():
    """A water sorted BEFORE Cys113 would shift the index; the conversion absorbs it."""
    pdb_names = ["ALA"] * 20 + ["HOH"] + ["ALA"] * 41 + ["CYS"] + ["ALA"] * 51
    cyx_index = pdb_names.index("CYS") + 1                    # 1-based, all residues
    pos = sum(1 for n in pdb_names[:cyx_index - 1] if n in br.PROTEIN_NAMES)
    protein = [n for n in pdb_names if n in br.PROTEIN_NAMES]
    assert protein[pos] == "CYS"
    assert pos != cyx_index - 1, "the water must shift the index, or this proves nothing"


# --------------------------------------------------------------------------
# convergence arithmetic
# --------------------------------------------------------------------------

def _colvar(path, distances, biases, dt=10.0):
    lines = ["#! FIELDS time d metad.bias"]
    lines += [f"{i*dt:.1f} {d:.4f} {b:.4f}"
              for i, (d, b) in enumerate(zip(distances, biases))]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_truncation_keeps_only_the_frames_within_the_window(tmp_path):
    src = _colvar(tmp_path / "COLVAR", [0.35] * 101, list(np.linspace(0, 100, 101)))
    dest = tmp_path / "cut"
    last = br.truncate_colvar(src, 500.0, dest)
    assert last == pytest.approx(500.0)
    kept = [l for l in dest.read_text().splitlines() if not l.startswith("#")]
    assert len(kept) == 51                      # t = 0 .. 500 ps at 10 ps
    assert dest.read_text().startswith("#! FIELDS"), "the header must survive"


def test_truncating_a_replica_that_stopped_early_is_a_no_op(tmp_path):
    """A COMMITTOR stop means there are no frames to cut; the answer must not change."""
    src = _colvar(tmp_path / "COLVAR", [0.35] * 31, [0.0] * 31)
    short = br.truncate_colvar(src, 10_000.0, tmp_path / "cut")
    assert short == pytest.approx(300.0)
    assert bpmd.analyse_replica(tmp_path / "cut", 0).frac_in_window == \
        pytest.approx(bpmd.analyse_replica(src, 0).frac_in_window)


def test_truncation_below_the_first_frame_is_an_error(tmp_path):
    src = _colvar(tmp_path / "COLVAR", [0.35] * 10, [0.0] * 10, dt=100.0)
    with pytest.raises(br.BPMDRunError, match="no frames"):
        br.truncate_colvar(src, -1.0, tmp_path / "cut")


def test_subsets_are_enumerated_exhaustively_when_there_are_few():
    rng = np.random.default_rng(0)
    assert len(br.draw_subsets(10, 2, 25, rng)) == 45 or \
        len(br.draw_subsets(10, 2, 100, rng)) == 45
    assert len(br.draw_subsets(4, 2, 100, rng)) == 6      # C(4,2), all of them
    assert len(br.draw_subsets(10, 5, 25, rng)) == 25     # C(10,5)=252, sampled


def test_subset_draws_are_reproducible():
    """A convergence verdict must not move because the analysis was re-run."""
    a = br.draw_subsets(10, 5, 25, np.random.default_rng(br.SUBSET_SEED))
    b = br.draw_subsets(10, 5, 25, np.random.default_rng(br.SUBSET_SEED))
    assert a == b


def test_convergence_report_recovers_a_protocol_that_does_not_matter(tmp_path):
    """Ten identical replicas: the score cannot depend on how many you use."""
    rows = []
    for i in range(10):
        p = _colvar(tmp_path / f"COLVAR{i}", [0.35] * 1001, [0.0] * 1001)
        rows.append({"status": "ok", "colvar": str(p)})
    df = br.convergence_report(rows, (1_000.0, 10_000.0), (2, 4, 10))
    assert df.score_sd.max() == pytest.approx(0.0, abs=1e-12)
    assert df.score_mean.nunique() == 1, "identical replicas, one answer"


def test_convergence_report_exposes_a_protocol_that_does(tmp_path):
    """Half the replicas hold and half escape immediately.

    A 2-replica protocol can draw two holders or two escapers, so its min and max
    must straddle the 10-replica answer. That spread is the finding the
    convergence test exists to produce.
    """
    rows = []
    for i in range(10):
        if i % 2:
            d, b = [0.35] * 501, [0.0] * 501
        else:
            d, b = [1.4] * 501, list(np.linspace(0, 50, 501))
        rows.append({"status": "ok",
                     "colvar": str(_colvar(tmp_path / f"C{i}", d, b))})
    df = br.convergence_report(rows, (5_000.0,), (2, 10))
    two = df[df.n_replicas == 2].iloc[0]
    ten = df[df.n_replicas == 10].iloc[0]
    assert two.score_min < ten.score_mean < two.score_max
    assert two.score_sd > 0


def test_convergence_refuses_a_single_replica(tmp_path):
    p = _colvar(tmp_path / "COLVAR", [0.35] * 10, [0.0] * 10)
    with pytest.raises(br.BPMDRunError, match="nothing to converge"):
        br.convergence_report([{"status": "ok", "colvar": str(p)}], (100.0,), (2,))


def test_a_replicate_on_disk_declares_the_length_it_was_run_at(tmp_path):
    """The 100 ps verification run and the 10 ns production run look identical on
    disk — both leave a COLVAR and a prod.gro. `nsteps` is what tells them apart,
    and pooling them would let a replica that was never given time to escape
    lower the escape count of a set that was."""
    mdp = tmp_path / "prod.mdp"
    mdp.write_text("integrator = md\ndt = 0.002\nnsteps = 50000\nref-t = 300\n")
    assert br._mdp_production_ps(mdp) == pytest.approx(100.0)
    mdp.write_text("nsteps = 5000000\n")
    assert br._mdp_production_ps(mdp) == pytest.approx(10_000.0)
    assert br._mdp_production_ps(tmp_path / "absent.mdp") is None


def test_failed_replicas_are_excluded_from_convergence_not_counted_as_zero(tmp_path):
    """A replicate that failed has not been measured; it must not score 0."""
    good = [{"status": "ok",
             "colvar": str(_colvar(tmp_path / f"C{i}", [0.35] * 101, [0.0] * 101))}
            for i in range(4)]
    df = br.convergence_report(good + [{"status": "failed: boom", "colvar": None}],
                               (1_000.0,), (2, 4))
    assert int(df.n_replicas.max()) == 4
