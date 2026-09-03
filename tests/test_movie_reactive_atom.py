"""The movie's "warhead->SG" readout must measure the warhead.

@twu383, 2026-09-02, looking at a mode selected under 3.0 A whose movie opened
at 5.35 A: *"audit why the movie shows a starting warhead distance higher than
3.0 and does not match the rmsd plots"*.

It did not match because the movie and the plot beneath it asked two different
questions. `mdprio_report.nac_series` located the reactive atom by the warhead
SMARTS; `elevation_report.surface_payload` took the atom NAMED `C10`. Measured
over 98 finished nac_v8 sweeps, `C10` was the reactive atom in **0** of them --
the readout was a median of 3.11 A away from the truth and up to 10.08 A.

`tests/test_crystal_pose_audit.py` had already written down why that fails:
"6VAJ calls it C10; the other five covalent Pin1 entries call the equivalent
atom C19, C14, C24, C12 and C3."
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))


def _er():
    try:
        import elevation_report as m
    except Exception as exc:                              # noqa: BLE001
        pytest.skip(f"elevation_report not importable here: {exc}")
    return m


def _movie(tmp_path: Path) -> Path:
    """Two frames: Cys113 SG, and a MOL residue whose atom names lie.

    `C10` is placed FAR from the sulfur and the real reactive atom (index 0)
    close to it, so any code that still selects by name reports a large distance
    and fails loudly rather than being off by a plausible amount.
    """
    p = tmp_path / "movie.pdb"
    L = []
    for m in (1, 2):
        L.append(f"MODEL     {m}")
        # Cys113 -> resid 63 with PIN1_OFFSET = 50
        L.append("ATOM      1  SG  CYS A  63       0.000   0.000   0.000  1.00  0.00")
        # MOL atom 0 -- the real warhead, 3.0 A away
        L.append("HETATM    2  C1  MOL A 300       3.000   0.000   0.000  1.00  0.00")
        L.append("HETATM    3  C2  MOL A 300       4.000   0.000   0.000  1.00  0.00")
        # MOL atom 2, NAMED C10, parked 20 A away
        L.append("HETATM    4  C10 MOL A 300      20.000   0.000   0.000  1.00  0.00")
        L.append("ENDMDL")
    p.write_text("\n".join(L) + "\n")
    return p


def test_selecting_the_reactive_atom_by_index_measures_the_right_atom(tmp_path,
                                                                     monkeypatch):
    er = _er()
    monkeypatch.setattr(er, "KEY_SITES", {}, raising=False)
    _t, dist, _l, _p = er.surface_payload(_movie(tmp_path), reactive_idx=0)
    assert dist[0] == pytest.approx(3.0, abs=0.01), (
        f"expected the MOL atom at index 0 (3.0 A), got {dist[0]} -- if this is "
        f"20.0 the C10 selection is back")


def test_there_is_no_implicit_default(tmp_path, monkeypatch):
    """The defect was a DEFAULT nobody chose, not a wrong argument.

    Catalogue #32/#35: a value that is legal, plausible and cannot announce that
    it is not what anyone wanted. Refusing to run without it is the fix.
    """
    er = _er()
    monkeypatch.setattr(er, "KEY_SITES", {}, raising=False)
    with pytest.raises(ValueError, match="reactive"):
        er.surface_payload(_movie(tmp_path))


def test_selecting_by_name_is_still_possible_but_must_be_asked_for(tmp_path,
                                                                   monkeypatch):
    """The old behaviour survives as an explicit request -- sulfopin needs it."""
    er = _er()
    monkeypatch.setattr(er, "KEY_SITES", {}, raising=False)
    _t, dist, _l, _p = er.surface_payload(_movie(tmp_path), reactive_name="C10")
    assert dist[0] == pytest.approx(20.0, abs=0.01)


def test_an_index_from_a_different_molecule_is_refused(tmp_path, monkeypatch):
    """Silently clamping would put a real number on the wrong atom."""
    er = _er()
    monkeypatch.setattr(er, "KEY_SITES", {}, raising=False)
    with pytest.raises(ValueError, match="outside"):
        er.surface_payload(_movie(tmp_path), reactive_idx=99)


def test_the_two_readouts_share_one_resolver():
    """The plot and the movie must not each locate the warhead their own way."""
    src = (REPO / "scripts" / "mdprio_report.py").read_text()
    assert "def reactive_atom(" in src, "the shared resolver is gone"
    # nac_series must USE it rather than re-implementing the SMARTS scan
    i = src.find("def nac_series(")
    body = src[i:i + 3000]
    assert "reactive_atom(" in body, (
        "nac_series no longer calls the shared resolver -- the plot and the "
        "movie can drift apart again")
    assert "for r in wh.itertuples()" not in body, (
        "nac_series has its own SMARTS scan again")


def test_no_caller_selects_the_warhead_by_a_hardcoded_name():
    """The whole class, not the one case.

    WALKS THE AST, NOT THE TEXT. The first version of this test read the source
    line by line and flagged the docstring in `surface_payload` that QUOTES the
    removed code as an example -- the same false positive the version-pin test
    produced when it flagged four docstrings. Prose describing a defect is not
    the defect; only executable code counts.
    """
    import ast
    offenders = []
    for f in (list((REPO / "scripts").glob("*.py"))
              + list((REPO / "shared").glob("*.py"))):
        try:
            tree = ast.parse(f.read_text(errors="replace"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            # `<something> == "C10"` and its mirror, anywhere in real code
            if not isinstance(node, ast.Compare):
                continue
            for op, comp in zip(node.ops, node.comparators):
                if not isinstance(op, (ast.Eq, ast.NotEq)):
                    continue
                if (isinstance(comp, ast.Constant)
                        and isinstance(comp.value, str)
                        and re.fullmatch(r"C\d{1,2}", comp.value)):
                    offenders.append(
                        f"{f.name}:{getattr(node, 'lineno', '?')}: "
                        f"compares against {comp.value!r}")
    assert not offenders, (
        "a ligand atom looks like it is being selected by name:\n  "
        + "\n  ".join(offenders))
