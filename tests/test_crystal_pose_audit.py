"""
Purpose: prove the crystal-pose audit refuses the two ways it could report a
         confidently wrong recovery number.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-12
Input: none — synthetic molecules and frames with a known answer
Output: pass/fail

WHY THIS EXISTS. #64 asked whether sulfopin's bottom-7% rank is a pose
GENERATION failure or a pose SELECTION failure, and both readings had evidence.
The measurement that settles it is a recovery rate against the crystal pose, and
a recovery rate has exactly two ways to be wrong while looking completely
normal:

  1. THE WRONG FRAME. #64 was blocked because the only sulfopin SDF the author
     found was origin-framed. Measured against it: best RMSD 12.07 A, 0 of 455
     poses within 2.5 A. That reads as a total docking failure and is a
     coordinate-system mismatch. `_check_frame` exists to raise instead.

  2. THE WRONG PAIRING. The criterion cannot be recomputed from the persisted
     cloud -- the screen measures the approach to the FLEXIBLE Cys113 sulfur,
     whose per-pose position is not in the archive -- so the screen's own rows
     are joined in BY ORDER. An order-based join is this project's defining
     defect shape, so it is checked two ways: the per-mode counts must fall out,
     and the geometry must agree pose by pose.

Both guards are tested by making them FAIL, per `how_this_project_breaks.md`:
a guard whose failure path is never exercised is a guard that passes for free.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import crystal_pose_audit as cpa                      # noqa: E402


# --------------------------------------------------------------------------
# the reactive atoms are found by CONNECTIVITY, not by the depositor's names
# --------------------------------------------------------------------------
def _adduct() -> Chem.Mol:
    """Sulfopin's deposited form: the free ligand minus its chlorine."""
    m = Chem.AddHs(Chem.MolFromSmiles("CC(=O)N(CC(C)(C)C)C1CCS(=O)(=O)C1"))
    AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
    return Chem.RemoveHs(m)


def test_reactive_atoms_are_derived_from_connectivity():
    """The reacting carbon is the terminal one on the amide carbonyl.

    6VAJ calls it C10; the other five covalent Pin1 entries call the equivalent
    atom C19, C14, C24, C12 and C3. Anything that selected it by name would
    pass on 6VAJ and quietly pick a different atom on the next structure.
    """
    m = _adduct()
    r, region = cpa.reactive_atoms(m)
    assert m.GetAtomWithIdx(r).GetSymbol() == "C"
    # its one heavy neighbour is the carbonyl carbon, which bears the amide N
    nbrs = m.GetAtomWithIdx(r).GetNeighbors()
    assert len(nbrs) == 1
    carbonyl = nbrs[0]
    assert {a.GetSymbol() for a in carbonyl.GetNeighbors()} == {"C", "O", "N"}
    # the region is reactive C, carbonyl C, carbonyl O, amide N -- four atoms
    assert len(region) == 4 and region[0] == r
    assert sorted(m.GetAtomWithIdx(i).GetSymbol() for i in region) == \
        ["C", "C", "N", "O"]


def test_reactive_atoms_refuses_an_ambiguous_ligand():
    """Two terminal carbons on one carbonyl is not something to resolve by luck."""
    m = Chem.AddHs(Chem.MolFromSmiles("CC(=O)NC"))     # N-methylacetamide
    AllChem.EmbedMolecule(m, randomSeed=1)
    m = Chem.RemoveHs(m)
    # the amide N here carries a terminal methyl too, but the carbonyl still has
    # exactly one -- so this must SUCCEED, and the refusal is exercised below.
    cpa.reactive_atoms(m)

    two = Chem.AddHs(Chem.MolFromSmiles("CC(=O)N(C)C(=O)C"))   # an imide
    AllChem.EmbedMolecule(two, randomSeed=2)
    with pytest.raises(ValueError, match="expected one amide carbonyl"):
        cpa.reactive_atoms(Chem.RemoveHs(two))


# --------------------------------------------------------------------------
# the frame guard
# --------------------------------------------------------------------------
def test_frame_guard_rejects_an_origin_framed_reference():
    """THE GUARD #64 NEEDED, exercised on the exact failure it exists for.

    A ligand whose reactive carbon is 13.5 A from the anchor sulfur is not a bad
    pose; it is a different coordinate system. Reporting it as recovery is how
    "0 of 455 poses within 2.5 A" got written down as though it were a docking
    result.
    """
    m = _adduct()
    r, _ = cpa.reactive_atoms(m)
    conf = m.GetConformer()
    # place the molecule at the origin and the anchor 13.5 A away
    xyz = np.array(conf.GetPositions())
    xyz -= xyz.mean(axis=0)
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, p.tolist())
    far = np.array(xyz[r]) + np.array([13.5, 0.0, 0.0])
    with pytest.raises(cpa.FrameError, match="not in the receptor's frame"):
        cpa._check_frame(m, r, far)


def test_frame_guard_accepts_a_bond_length():
    """And it must PASS at the distance a covalent adduct actually sits at."""
    m = _adduct()
    r, _ = cpa.reactive_atoms(m)
    sg = np.array(m.GetConformer().GetPositions()[r]) + np.array([1.78, 0, 0])
    assert cpa._check_frame(m, r, sg) == pytest.approx(1.78, abs=1e-6)


# --------------------------------------------------------------------------
# the join
# --------------------------------------------------------------------------
def _cloud_and_table(n_per_mode=(3, 4), scramble=False):
    """A cloud and the screen rows it came from, with a known pairing.

    The cloud is written the way `nac_screen_v2` writes it -- grouped by mode,
    conformer order preserved within a mode -- and the table carries the
    noise mode the cloud omits.
    """
    rows, cloud = [], []
    idx = 0
    for mode, n in enumerate(n_per_mode):
        for k in range(n):
            rows.append({"ident": "x", "mode": mode, "pose_idx": idx,
                         "energy": -5.0 - k, "energy_rank": idx + 1,
                         "distance": 3.0 + 0.1 * idx, "angle": 100.0 + idx,
                         "viable": bool(k == 0), "in_range": True})
            cloud.append({"mode": mode, "pose_rank": len(cloud) + 1,
                          "distance_rigid_sg": 3.0 + 0.1 * idx + 0.05})
            idx += 1
    rows.append({"ident": "x", "mode": -1, "pose_idx": idx, "energy": -1.0,
                 "energy_rank": idx + 1, "distance": 9.9, "angle": 10.0,
                 "viable": False, "in_range": False})
    tbl = pd.DataFrame(rows)
    cl = pd.DataFrame(cloud)
    if scramble:
        cl["distance_rigid_sg"] = cl["distance_rigid_sg"].values[::-1]
    return cl, tbl


def test_join_drops_the_noise_mode_and_pairs_in_order():
    cl, tbl = _cloud_and_table()
    out = cpa._join_screen_rows(cl, tbl, "x")
    assert len(out) == len(cl)
    assert (out["mode"].values == cl["mode"].values).all()
    # the screen's distance arrives, offset by the sidechain shift only
    assert np.allclose(out["distance"] + 0.05, out["distance_rigid_sg"])


def test_join_refuses_a_cloud_that_is_a_different_run():
    """A cloud of a different size is two runs, not one -- and must not join."""
    cl, tbl = _cloud_and_table()
    with pytest.raises(ValueError, match="different runs"):
        cpa._join_screen_rows(cl.iloc[:-1], tbl, "x")


def test_join_check_catches_a_scrambled_pairing():
    """THE TEST THAT MATTERS: counts alone cannot see a within-mode scramble.

    A pairing that stays inside its mode reproduces every per-mode count
    exactly, so a count-based check passes it. Only the pose-by-pose geometry
    notices, which is why the correlation test exists.
    """
    cl, tbl = _cloud_and_table(n_per_mode=(6, 6), scramble=True)
    out = cpa._join_screen_rows(cl, tbl, "x")
    agg = pd.DataFrame([
        {"mode": m, "n_poses_mode": int((out["mode"] == m).sum()),
         "n_in_range": int(out.loc[out["mode"] == m, "in_range"].sum()),
         "n_viable": int(out.loc[out["mode"] == m, "viable"].sum())}
        for m in sorted(out["mode"].unique())])
    with pytest.raises(ValueError, match="not paired correctly"):
        cpa._check_join(out, agg)


def test_join_check_passes_the_correct_pairing():
    cl, tbl = _cloud_and_table(n_per_mode=(6, 6))
    out = cpa._join_screen_rows(cl, tbl, "x")
    agg = pd.DataFrame([
        {"mode": m, "n_poses_mode": int((out["mode"] == m).sum()),
         "n_in_range": int(out.loc[out["mode"] == m, "in_range"].sum()),
         "n_viable": int(out.loc[out["mode"] == m, "viable"].sum())}
        for m in sorted(out["mode"].unique())])
    got = cpa._check_join(out, agg)
    assert got["join_counts_reproduced"] is True
    assert got["join_distance_r"] > 0.99
