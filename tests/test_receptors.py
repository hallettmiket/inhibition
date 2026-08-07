"""Tests for shared.receptors — the receptor-identity guard.

These tests exist because "3IKD" names two different structures on this box, and
2.0.0's defect catalogue is almost entirely values taken by name rather than by
identity. The point of the guard is that it FAILS rather than warns, so the tests
assert the failure as much as the success.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import receptors as R  # noqa: E402


def test_ian_and_deposited_are_different_files():
    """The whole reason this module exists."""
    assert R.IAN_3IKD_SRC != R.DEPOSITED_3IKD


def test_resolve_returns_the_chemists_file():
    p = R.resolve_3ikd_ian()
    assert p == R.IAN_3IKD_SRC
    assert p.exists()


def test_resolve_verifies_the_digest():
    """A pinned digest that does not match must raise, not warn."""
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(R, "IAN_3IKD_SHA", "0" * 64)
        with pytest.raises(R.ReceptorIdentityError, match="changed underneath"):
            R.resolve_3ikd_ian()


def test_verify_false_skips_the_check():
    """Escape hatch exists, but is opt-in and never the default."""
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(R, "IAN_3IKD_SHA", "0" * 64)
        assert R.resolve_3ikd_ian(verify=False) == R.IAN_3IKD_SRC


def test_missing_source_names_the_other_structure():
    """The error must say WHICH file is missing and what it is not the same as."""
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(R, "IAN_3IKD_SRC", R.DATA / "does_not_exist.pdb")
        with pytest.raises(R.ReceptorIdentityError) as e:
            R.resolve_3ikd_ian()
        assert "deposited" in str(e.value)


def test_noligand_variant_resolves():
    p = R.resolve_3ikd_ian(noligand=True)
    assert p == R.IAN_3IKD_NOLIGAND
    assert p.exists()


@pytest.mark.parametrize("attr,expected", [
    ("IAN_3IKD_SRC", "3IKD_ian"),
    ("IAN_3IKD_NOLIGAND", "3IKD_ian_noligand"),
    ("DEPOSITED_3IKD", "3IKD_deposited"),
])
def test_describe_identifies_each_receptor(attr, expected):
    path = getattr(R, attr)
    if not path.exists():
        pytest.skip(f"{attr} not present on this machine")
    assert R.describe(path)["receptor"] == expected


def test_describe_refuses_to_guess_from_the_filename():
    """A path merely CONTAINING '3IKD' must not be reported as 3IKD_ian.

    This is the failure mode being guarded: the filename is exactly what is
    ambiguous, so it can never be the evidence.
    """
    assert R.describe(Path("/tmp/some_other_3IKD.pdb"))["receptor"] == "UNRECOGNISED"


def test_the_pinned_digest_matches_what_is_on_disk():
    """Guards against the pin itself being wrong — a pin nobody checks is decoration."""
    if not R.IAN_3IKD_SRC.exists():
        pytest.skip("3IKD_ian not present on this machine")
    h = hashlib.sha256(R.IAN_3IKD_SRC.read_bytes()).hexdigest()
    assert h == R.IAN_3IKD_SHA
