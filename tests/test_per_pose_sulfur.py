"""Every pose is measured against ITS OWN Cys113 sulfur.

THE DEFECT (D0109). Cys113 is docked as a flexible sidechain, so each pose has
its own SG. `nac_screen.sg_position` returned the FIRST docked model's sulfur --
its docstring saying, correctly, that the position "must come from the pose being
measured" -- and `measure_poses` broadcast that one value across all 640
conformers. Poses 2..N were measured against where the sulfur sat in pose 1.

WHY IT WAS INVISIBLE FOR SO LONG. The error is small for most poses (median
0.18 A of SG movement) and one-directional, so distances were plausible and
slightly too short. It only became visible at the tail: 7 of 200 poses came out
below 1.81 A, shorter than a C-S bond. Nothing raised, because a short distance
is a legal float.

WHAT THESE TESTS DO NOT DO. They do not check that `sg_positions` parses a DLG
-- that needs a docking run. They check the CONTRACT that made the silent
broadcast possible: `measure_poses` must refuse one sulfur for many conformers.
An API that cannot be misused this way is the fix; the call-site edits are just
today's instance.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from shared import nac_criterion as nac                    # noqa: E402


def _mol(n_conformers: int, mechanism_smarts: str = "CCCl"):
    """A tiny molecule with `n_conformers` DISTINCT conformers."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(mechanism_smarts))
    AllChem.EmbedMultipleConfs(m, numConfs=n_conformers, randomSeed=42)
    if m.GetNumConformers() < n_conformers:      # embedding is not guaranteed
        pytest.skip("RDKit produced fewer conformers than requested")
    return m


def test_one_sulfur_for_many_poses_is_refused():
    """The exact call the screen used to make must now raise."""
    m = _mol(5)
    match = (0, 2)                                # C, Cl -- sn2_displacement
    with pytest.raises(ValueError, match="one SG position was given"):
        nac.measure_poses(m, match, "sn2_displacement",
                          np.array([1.0, 2.0, 3.0]))


def test_a_static_sulfur_must_be_asked_for_by_name():
    """`allow_static_sg` is the opt-in, so the choice appears at the call site."""
    m = _mol(5)
    out = nac.measure_poses(m, (0, 2), "sn2_displacement",
                            np.array([1.0, 2.0, 3.0]), allow_static_sg=True)
    assert len(out) == 5


def test_a_wrong_length_sulfur_array_is_refused():
    """Pairing pose j with pose k's sulfur is worse than the original defect."""
    m = _mol(5)
    with pytest.raises(ValueError, match="expected"):
        nac.measure_poses(m, (0, 2), "sn2_displacement", np.zeros((3, 3)))


def test_per_pose_sulfur_actually_changes_the_answer():
    """The per-pose array must be USED, not accepted and then ignored.

    An implementation that took the array and indexed `sg[0]` every time would
    pass every test above. So this feeds two sulfur sets that differ only after
    pose 0 and requires the measured distances to differ -- which they cannot if
    only the first row is read.
    """
    m = _mol(6)
    match = (0, 2)
    base = np.tile(np.array([5.0, 0.0, 0.0]), (6, 1))
    moved = base.copy()
    moved[1:] += np.array([3.0, 0.0, 0.0])        # every pose but the first

    d_base = [r.distance for r in
              nac.measure_poses(m, match, "sn2_displacement", base)]
    d_moved = [r.distance for r in
               nac.measure_poses(m, match, "sn2_displacement", moved)]

    assert d_base[0] == pytest.approx(d_moved[0]), (
        "pose 0's sulfur was identical in both runs; its distance must match")
    differing = sum(1 for a, b in zip(d_base[1:], d_moved[1:])
                    if abs(a - b) > 1e-6)
    assert differing == 5, (
        f"only {differing} of 5 later poses changed when their sulfurs moved "
        f"3 A. The per-pose array is being accepted and then not indexed per "
        f"pose -- which is the D0109 defect wearing the new signature.")


def test_singular_sg_position_is_still_reachable_but_named():
    """`sg_position` must survive as a NAME, so its misuse is legible.

    Deleting it would turn every stale caller into an AttributeError somewhere
    unrelated. Keeping it means a reader who finds it sees the deprecation and
    the record number.
    """
    from scripts import nac_screen as ns                    # noqa: F401
    assert hasattr(ns, "sg_positions"), "the plural accessor must exist"
    assert hasattr(ns, "sg_position"), (
        "the singular accessor was deleted rather than deprecated; stale "
        "callers now fail somewhere that does not name the reason")
    assert "D0109" in (ns.sg_position.__doc__ or ""), (
        "sg_position's docstring must name the record, or the next reader "
        "cannot tell why the plural exists")
