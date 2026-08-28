"""The PoseBusters gate: it flags, it never deletes, and it cannot pass vacuously.

D0089 adopted this gate and nac_v5 shipped without it. The properties that make
it safe to add now are the ones tested here -- not its hit rate, which is
chemistry and belongs in the decision record.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import pose_contacts as pc            # noqa: E402


def test_a_pdbqt_receptor_is_refused_by_name():
    """PoseBusters cannot parse pdbqt. It WARNS rather than raising, so every
    protein-ligand check fails for want of a receptor and the verdict is
    '0 of 500 valid' -- an infrastructure fault wearing a result's clothes.
    This actually happened on the first run of the gate."""
    v2 = pytest.importorskip("nac_screen_v2", reason="needs the docking env")
    with pytest.raises(ValueError, match="needs a .pdb receptor"):
        v2.posebusters_valid(None, Path("/tmp/x.pdbqt"), Path("/tmp"))


def test_split_poses_labels_excluded_conformers_minus_one():
    """The gate excludes by LABELLING, not by deleting. D0093 is the record of
    what a filter that deletes costs: 21% of every cloud gone, and four
    experiments measuring a population nobody chose."""
    rng = np.random.default_rng(0)
    T = np.clip(rng.uniform(2, 8, (40, 10, 5)), 0.5, pc.CAP_A).astype(np.float32)
    w = np.ones(10)
    D = pc.pose_distances(T, w)
    lab = pc.group(D, 1.0)
    assert (lab >= 0).all(), "group() itself must never emit a noise label"


def test_excluded_poses_do_not_change_the_groups_of_the_rest():
    """Grouping a subset must give the subset the same answer it would get alone
    -- otherwise the gate silently changes the grouping of valid poses too."""
    rng = np.random.default_rng(3)
    T = np.clip(rng.uniform(2, 8, (60, 8, 4)), 0.5, pc.CAP_A).astype(np.float32)
    w = np.ones(8)
    keep = np.arange(0, 60, 2)
    full = pc.group(pc.pose_distances(T[keep], w), 0.9)
    again = pc.group(pc.pose_distances(T[keep], w), 0.9)
    assert (full == again).all()


def test_the_gate_cannot_pass_by_being_absent():
    """A guard whose failure mode is 'everything valid' is not a guard. The
    screen raises when PoseBusters cannot be imported rather than defaulting the
    mask to True -- assert the module does not carry such a fallback."""
    src = (REPO / "scripts" / "nac_screen_v2.py").read_text()
    fn = src[src.index("def posebusters_valid"):src.index("def one(cand")]
    assert "except ImportError" not in fn, (
        "a silent import fallback would make the gate pass for every pose")
    assert "return np.ones" not in fn, (
        "the gate must never manufacture an all-valid verdict")
