"""
Purpose: Every docked mode must be readable, scored, and correctly ordered.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: integration/app/pose3d.py, synthetic pose files + the real ones
Output: pass/fail

WHY THIS EXISTS. Issue #3 asked to "show all poses for a molecule not just one".
The viewer drew model 1 and stopped. For Vina that was at least the best-scoring
mode — Vina sorts by affinity. For gnina it was NOT: the SDF comes back sorted
by CNNscore, and the candidate frame's `affinity_kcal` is read off the FIRST
record, so on this build 11 of 25 T_3 shortlist entries and 6 of 27 T_4 entries
have a pose in the SAME FILE with a better `minimizedAffinity` than the one that
was drawn — by up to 2.86 kcal/mol.

That is not a missing feature, it is a hidden disagreement between two scores
that the pipeline reports side by side. So the tests below pin three things: all
modes are read, each carries its own scores, and the direction of every score is
respected (kcal/mol lower-is-better; CNN* higher-is-better). Getting that
direction backwards would present the worst pose as the best, which is the one
failure that looks completely normal on screen.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "integration" / "app"))
import pose3d as p3d  # noqa: E402


# --- synthetic fixtures: format handling, independent of what is on disk ----

def _sdf(n: int = 3, affinities=(-5.0, -7.0, -6.0), cnn=(0.9, 0.4, 0.2)) -> str:
    """A minimal multi-record SDF shaped like gnina's, blank title lines and all."""
    out = []
    for i in range(n):
        out += [
            "", "", "",
            "  1  0  0  0  0  0  0  0  0  0999 V2000",
            f"  {-12.0 + i:8.4f}  {-35.5:8.4f}  {9.0:8.4f} C   0  0  0  0  0",
            "M  END",
            "> <minimizedAffinity>", f"{affinities[i]}", "",
            "> <CNNscore>", f"{cnn[i]}", "",
            "$$$$",
        ]
    return "\n".join(out) + "\n"


def _pdbqt(affinities=(-8.3, -8.1, -7.9)) -> str:
    out = []
    for i, a in enumerate(affinities, 1):
        out += [
            f"MODEL {i}",
            f"REMARK VINA RESULT:      {a}      0.000      0.000",
            f"ATOM      1  C   UNL     1     {-12.0 + i:7.3f} {-35.500:7.3f} "
            f"{9.000:7.3f}  0.00  0.00    +0.000 C ",
            "ENDMDL",
        ]
    return "\n".join(out) + "\n"


def test_all_sdf_records_are_read_not_just_the_first(tmp_path):
    f = tmp_path / "d_test_docked.sdf"
    f.write_text(_sdf())
    poses = p3d.read_poses(f)
    assert [p.index for p in poses] == [1, 2, 3]
    assert [p.scores["minimizedAffinity"] for p in poses] == [-5.0, -7.0, -6.0]
    assert [p.scores["CNNscore"] for p in poses] == [0.9, 0.4, 0.2]


def test_sdf_blocks_keep_their_blank_header_lines(tmp_path):
    """A split-then-strip parse eats the title line and shifts the counts row.

    3Dmol then reads the atom count out of the wrong row and draws nothing —
    which looks exactly like a molecule that failed to load.
    """
    f = tmp_path / "d_test_docked.sdf"
    f.write_text(_sdf(n=2))
    for pose in p3d.read_poses(f):
        lines = pose.text.splitlines()
        assert lines[0] == "" and lines[1] == "" and lines[2] == ""
        assert "V2000" in lines[3], "the counts line must stay on line 4"


def test_all_pdbqt_models_are_read(tmp_path):
    f = tmp_path / "t1_test_out.pdbqt"
    f.write_text(_pdbqt())
    poses = p3d.read_poses(f)
    assert [p.index for p in poses] == [1, 2, 3]
    assert [p.scores["vina_affinity"] for p in poses] == [-8.3, -8.1, -7.9]
    assert all(p.text.startswith("MODEL") and "ENDMDL" in p.text for p in poses)


def test_a_single_model_file_without_a_wrapper_is_still_one_pose(tmp_path):
    f = tmp_path / "t1_bare_out.pdbqt"
    f.write_text("ATOM      1  C   UNL     1     -12.000 -35.500   9.000\n")
    assert len(p3d.read_poses(f)) == 1


def test_score_direction_is_respected_per_score(tmp_path):
    """kcal/mol is lower-better; CNN* is higher-better. Mixing them inverts it."""
    f = tmp_path / "d_test_docked.sdf"
    f.write_text(_sdf())
    table = p3d.pose_score_table(p3d.read_poses(f))
    starred = {row["pose"]: [k for k, v in row.items()
                             if k != "pose" and v and "★" in str(v)]
               for row in table}
    assert starred[2] == ["minimizedAffinity"], "best affinity is the MOST negative"
    assert starred[1] == ["CNNscore"], "best CNNscore is the LARGEST"


def test_hidden_better_pose_reports_the_gap(tmp_path):
    f = tmp_path / "d_test_docked.sdf"
    f.write_text(_sdf())
    hit = p3d.hidden_better_pose(p3d.read_poses(f), shown=1,
                                 score="minimizedAffinity")
    assert hit is not None
    idx, gap = hit
    assert idx == 2 and abs(gap - 2.0) < 1e-9


def test_hidden_better_pose_is_none_when_pose_one_already_wins(tmp_path):
    """Vina sorts by affinity, so this must NOT cry wolf on T_1/T_2."""
    f = tmp_path / "t1_test_out.pdbqt"
    f.write_text(_pdbqt())
    assert p3d.hidden_better_pose(p3d.read_poses(f), shown=1,
                                  score="vina_affinity") is None


def test_a_non_covalent_pose_gets_no_bond_drawn(tmp_path):
    """Drawing a link into T_1/T_2 would assert a mechanism they do not have."""
    # Readable, not present -- see `pose3d.receptor_readable`.
    if not p3d.receptor_readable():
        pytest.skip("prepared receptor is not readable by this user "
                    "(absent, or denied by the data root's ACL)")
    f = tmp_path / "t1_far_out.pdbqt"
    f.write_text(
        "MODEL 1\nREMARK VINA RESULT:      -8.0      0.000      0.000\n"
        "ATOM      1  C   UNL     1      50.000  50.000  50.000  0.00  0.00\n"
        "ENDMDL\n")
    assert p3d.covalent_link(p3d.read_poses(f)[0]) is None


# --- the real files --------------------------------------------------------

def _first_pose_file(approach: str, pattern: str) -> Path | None:
    """The first pose file this process can actually OPEN, or None.

    Readability, not existence -- the same distinction as
    `pose3d.receptor_readable`. Returning a matched-but-unreadable path made
    the caller's `pytest.skip("no docked poses on disk")` miss, and the test
    then failed with a bare `PermissionError` raised from inside the parser,
    which reads as a broken parser rather than as an access problem.
    """
    d = p3d.DOCKING_DIRS.get(approach)
    if d is None or not d.is_dir():
        return None
    return next((p for p in sorted(d.glob(pattern)) if os.access(p, os.R_OK)),
                None)


@pytest.mark.parametrize("approach", ["t3", "t4"])
def test_real_covalent_poses_are_bonded_to_cys113(approach):
    """gnina was told to bond to A:113:SG, so EVERY mode must sit on it.

    A mode that drifted off the SG would mean the covalent restraint did not
    hold, which is a pipeline failure the viewer would otherwise render as a
    perfectly ordinary picture.
    """
    # Readable, not present -- see `pose3d.receptor_readable`.
    if not p3d.receptor_readable():
        pytest.skip("prepared receptor is not readable by this user "
                    "(absent, or denied by the data root's ACL)")
    f = _first_pose_file(approach, "d_*_docked.sdf")
    if f is None:
        pytest.skip(f"{approach}: no docked poses on disk")
    poses = p3d.read_poses(f)
    assert len(poses) > 1, "gnina writes several modes; one means a parse bug"
    for pose in poses:
        link = p3d.covalent_link(pose)
        assert link is not None, (
            f"{f.name} pose {pose.index} has no atom within "
            f"{p3d.COVALENT_BOND_MAX_A} A of the Cys113 SG")
        assert 1.4 < link[2] < p3d.COVALENT_BOND_MAX_A, (
            f"link length {link[2]:.2f} A is not a C-S bond")


@pytest.mark.parametrize("approach", ["t1", "t2"])
def test_real_vina_poses_are_ordered_best_first(approach):
    f = _first_pose_file(approach, "poses/*_out.pdbqt")
    if f is None:
        pytest.skip(f"{approach}: no docked poses on disk")
    poses = p3d.read_poses(f)
    aff = [p.scores["vina_affinity"] for p in poses if "vina_affinity" in p.scores]
    assert len(aff) > 1
    assert aff == sorted(aff), (
        "Vina is documented to return modes best-first; if that ever changes, "
        "the caption claiming pose 1 is the best mode becomes false")


def test_the_export_bundle_carries_everything_needed_to_reopen_the_view():
    """The hand-off is only useful if it is self-contained."""
    import io
    import zipfile

    # Readable, not present -- see `pose3d.receptor_readable`.
    if not p3d.receptor_readable():
        pytest.skip("prepared receptor is not readable by this user "
                    "(absent, or denied by the data root's ACL)")
    f = _first_pose_file("t3", "d_*_docked.sdf")
    if f is None:
        pytest.skip("no covalent poses on disk")
    z = zipfile.ZipFile(io.BytesIO(p3d.pose_bundle(f, covalent=True)))
    names = set(z.namelist())
    assert {f.name, p3d.RECEPTOR.name, "view_pin1.pml", "README.txt"} == names
    # The pose file must go over WHOLE — exporting one mode would reproduce the
    # very truncation this feature exists to undo.
    assert len(p3d.read_poses(f)) == z.read(f.name).decode().count("$$$$")
    readme = z.read("README.txt").decode()
    assert "CNNscore" in readme, "the reader must be told what the order means"


def test_pose_html_renders_the_surface_labels_and_every_requested_mode(tmp_path):
    """A Streamlit change that only compiles is not done."""
    pytest.importorskip("py3Dmol")
    # Readable, not present -- see `pose3d.receptor_readable`.
    if not p3d.receptor_readable():
        pytest.skip("prepared receptor is not readable by this user "
                    "(absent, or denied by the data root's ACL)")
    f = _first_pose_file("t3", "d_*_docked.sdf")
    if f is None:
        pytest.skip("no covalent poses on disk")

    n = len(p3d.read_poses(f))
    html = p3d.pose_html(f, show=tuple(range(1, n + 1)), height=400)

    # one receptor + every requested pose
    assert html.count("addModel(") == n + 1
    # a grey lining surface plus one per sub-pocket
    assert html.count("addSurface(") == len(p3d.SUBPOCKETS) + 1
    for sp in p3d.SUBPOCKETS:
        assert sp.label in html, f"{sp.key} is not labelled in the render"
        assert f'"color": "{sp.colour}"' in html
    assert "addCylinder(" in html, "the covalent link must be drawn"
    assert "covalent Cys113 SG" in html
    assert html.count("addCylinder(") == n, "every mode's link, not just one"

    single = p3d.pose_html(f, show=(1,), height=400)
    assert single.count("addModel(") == 2, "receptor + exactly one pose"
    assert len(single) < len(html)
