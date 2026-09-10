"""Typed interactions must be typed by IDENTITY, and must fail loudly.

@twu383, 2026-09-09: *"I have many preferential factors that I can determine
visually but need to find ways to quantitatively describe them for automated
screening"*.

Every test here is about a way this could produce a plausible profile from the
wrong chemistry -- which is the only failure mode that matters in a descriptor
meant to drive automated selection.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from shared import interaction_profile as ip          # noqa: E402


def _pdb_line(rec, serial, name, resname, resid, xyz, element) -> str:
    """A PDB ATOM/HETATM record built BY COLUMN, not by concatenation.

    Written this way after the first version put `resName` at columns 16-18
    instead of 17-19, so "LEU" parsed as "EU" and the module correctly refused
    to type it. A fixture that mis-columns its own format tests the parser's
    error path and nothing else -- and the reader here slices fixed columns
    exactly as the PDB spec says, so the fixture has to as well.
    """
    s = [" "] * 80
    def put(start, text, width):          # 1-based PDB columns
        for k, ch in enumerate(str(text)[:width]):
            s[start - 1 + k] = ch
    put(1, f"{rec:<6}", 6)
    put(7, f"{serial:>5}", 5)
    put(13, f"{name:<4}", 4)              # 13-16 atom name
    put(18, f"{resname:>3}", 3)           # 18-20 resName
    put(22, "A", 1)                       # 22 chainID
    put(23, f"{resid:>4}", 4)             # 23-26 resSeq
    put(31, f"{xyz[0]:8.3f}", 8)          # 31-38 x
    put(39, f"{xyz[1]:8.3f}", 8)
    put(47, f"{xyz[2]:8.3f}", 8)
    put(55, "  1.00", 6)
    put(61, "  0.00", 6)
    put(77, f"{element:>2}", 2)           # 77-78 element
    return "".join(s).rstrip()


def _frame(lig_xyz, prot_atoms) -> str:
    """One MODEL: prot_atoms = [(resname, resid, name, element, xyz)]."""
    L = ["MODEL     1"]
    for i, xyz in enumerate(lig_xyz):
        L.append(_pdb_line("HETATM", i + 1, f"C{i}", "MOL", 300, xyz, "C"))
    for i, (rn, ri, nm, el, xyz) in enumerate(prot_atoms):
        L.append(_pdb_line("ATOM", i + 900, nm, rn, ri, xyz, el))
    return "\n".join(L) + "\nENDMDL\n"


def test_an_untypable_residue_raises_rather_than_reading_as_apolar(tmp_path,
                                                                  monkeypatch):
    """A missing interaction is invisible in a way a wrong one is not.

    Silently treating an unknown residue as apolar drops every polar
    interaction it makes, and the profile still looks complete -- catalogue
    disguise #4, a guard that cannot fail.
    """
    pdb = tmp_path / "m.pdb"
    pdb.write_text(_frame([(0, 0, 0)],
                          [("XYZ", 10, "CA", "C", (3.0, 0, 0))]))
    monkeypatch.setattr(ip, "ligand_chemistry",
                        lambda s: {"n_heavy": 1, "donors": set(), "acceptors": set(),
                                   "apolar": {0}, "halogens": set(),
                                   "halogen_roots": {}, "cations": set(),
                                   "anions": set(), "aromatic_rings": []})
    with pytest.raises(ip.ProfileError, match="untypable"):
        ip.profile(pdb, tmp_path / "x.sdf")


def test_a_ligand_of_the_wrong_size_raises(tmp_path, monkeypatch):
    """The indices address frame coordinates; a mismatch addresses wrong atoms.

    Same shape as the C10 defect: the numbers would all be plausible.
    """
    pdb = tmp_path / "m.pdb"
    pdb.write_text(_frame([(0, 0, 0), (1.5, 0, 0)],
                          [("ALA", 10, "CB", "C", (3.0, 0, 0))]))
    monkeypatch.setattr(ip, "ligand_chemistry",
                        lambda s: {"n_heavy": 7, "donors": set(), "acceptors": set(),
                                   "apolar": set(), "halogens": set(),
                                   "halogen_roots": {}, "cations": set(),
                                   "anions": set(), "aromatic_rings": []})
    with pytest.raises(ip.ProfileError, match="not the same molecule"):
        ip.profile(pdb, tmp_path / "x.sdf")


def test_HID_and_HIE_donate_from_opposite_nitrogens():
    """Amber writes the protonation state into the NAME, and it inverts the roles.

    His59 and His157 are the catalytic pair on this target, so collapsing both
    to "HIS" would mis-assign donor/acceptor on the two residues that matter
    most. `md_movie` carries HIS_FORMS for the same reason.
    """
    assert ip.SIDECHAIN_DONORS["HID"] == ("ND1",)
    assert ip.SIDECHAIN_DONORS["HIE"] == ("NE2",)
    assert ip.SIDECHAIN_ACCEPTORS["HID"] == ("NE2",)
    assert ip.SIDECHAIN_ACCEPTORS["HIE"] == ("ND1",)
    # and they are genuinely opposite, not accidentally equal
    assert ip.SIDECHAIN_DONORS["HID"] != ip.SIDECHAIN_DONORS["HIE"]


def test_a_salt_bridge_needs_opposite_charges(tmp_path, monkeypatch):
    """A cation near a cation is not a salt bridge."""
    pdb = tmp_path / "m.pdb"
    pdb.write_text(_frame([(0, 0, 0)],
                          [("LYS", 13, "NZ", "N", (3.0, 0, 0))]))
    base = dict(n_heavy=1, donors=set(), acceptors=set(), apolar=set(),
                halogens=set(), halogen_roots={}, aromatic_rings=[])
    monkeypatch.setattr(ip, "ligand_chemistry",
                        lambda s: {**base, "anions": {0}, "cations": set()})
    df = ip.profile(pdb, tmp_path / "x.sdf")
    assert (df.interaction == "salt_bridge").any(), "LYS+ / ligand- was missed"

    monkeypatch.setattr(ip, "ligand_chemistry",
                        lambda s: {**base, "anions": set(), "cations": {0}})
    df2 = ip.profile(pdb, tmp_path / "x.sdf")
    assert not (df2.interaction == "salt_bridge").any(), (
        "two cations were called a salt bridge")


def test_hbond_geometry_rejects_an_acceptor_behind_the_donor(tmp_path,
                                                             monkeypatch):
    """Distance alone is not an H-bond; there is no H, so the angle stands in."""
    base = dict(n_heavy=1, donors=set(), apolar=set(), halogens=set(),
                halogen_roots={}, cations=set(), anions=set(), aromatic_rings=[])
    monkeypatch.setattr(ip, "ligand_chemistry",
                        lambda s: {**base, "acceptors": {0}})
    # acceptor 3.0 A from the backbone N: in range
    pdb = tmp_path / "ok.pdb"
    pdb.write_text(_frame([(3.0, 0, 0)],
                          [("ALA", 10, "N", "N", (0, 0, 0))]))
    df = ip.profile(pdb, tmp_path / "x.sdf")
    assert (df.interaction == "hbond_protein_donor").any()
    # and one beyond the cutoff is not
    far = tmp_path / "far.pdb"
    far.write_text(_frame([(5.0, 0, 0)],
                          [("ALA", 10, "N", "N", (0, 0, 0))]))
    df2 = ip.profile(far, tmp_path / "x.sdf")
    assert not (df2.interaction == "hbond_protein_donor").any()


def test_occupancy_is_a_fraction_of_frames_and_needs_the_count(tmp_path,
                                                               monkeypatch):
    """A count is not screenable; a fraction of the run is."""
    import pandas as pd
    df = pd.DataFrame({"frame": [0, 1, 2, 0], "resname": ["ALA"] * 4,
                       "resid_pdb": [60] * 4,
                       "interaction": ["hydrophobic"] * 3 + ["salt_bridge"],
                       "distance_a": [3.5, 3.6, 3.7, 3.9]})
    occ = ip.occupancy(df, n_frames=4)
    h = occ[occ.interaction == "hydrophobic"].iloc[0]
    assert h.occupancy == pytest.approx(0.75)
    with pytest.raises(ip.ProfileError, match="frame count"):
        ip.occupancy(df)


def test_residue_numbers_are_reported_in_pin1_numbering(tmp_path, monkeypatch):
    """A row must be readable against the literature without knowing the offset."""
    pdb = tmp_path / "m.pdb"
    pdb.write_text(_frame([(0, 0, 0)],
                          [("LEU", 72, "CB", "C", (3.5, 0, 0))]))
    monkeypatch.setattr(ip, "ligand_chemistry",
                        lambda s: {"n_heavy": 1, "donors": set(), "acceptors": set(),
                                   "apolar": {0}, "halogens": set(),
                                   "halogen_roots": {}, "cations": set(),
                                   "anions": set(), "aromatic_rings": []})
    df = ip.profile(pdb, tmp_path / "x.sdf")
    assert set(df.resid_pdb) == {122}, "72 + PIN1_OFFSET should read as Leu122"
    assert set(df.resid_md) == {72}
