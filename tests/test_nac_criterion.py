"""
Purpose: the near-attack criterion is tested on hand-built geometry with a known answer.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05

WHY HAND-BUILT AND NOT DOCKED POSES. A docked pose's "correct" NAC verdict is
whatever this module says it is, so testing against one would be circular. Every
case here places atoms at coordinates whose angle is known from the arithmetic,
so the test can disagree with the implementation.

The distance cases carry the module's central claim -- that a near-attack
conformation is a van der Waals CONTACT and not a bond. A pose at 1.55 A with a
perfect 180 deg approach must be REJECTED: that is past the transition state and
unreachable by a molecule still carrying its leaving group. Reactive docking
produces exactly such poses by the hundred, so this is the assertion that keeps
them out.
"""

from __future__ import annotations

import numpy as np
import pytest

from shared import nac_criterion as nac


# --------------------------------------------------------------------------
# SN2 — backside attack, anti to the leaving group
# --------------------------------------------------------------------------

def sn2_coords(angle_deg: float, dist: float = 3.5):
    """C at the origin, Cl on +x, S placed at `angle_deg` from the C->Cl bond.

    angle_deg = 180 puts the sulfur exactly opposite the chlorine.
    """
    c = np.zeros(3)
    cl = np.array([1.78, 0.0, 0.0])
    t = np.radians(angle_deg)
    sg = np.array([dist * np.cos(t), dist * np.sin(t), 0.0])
    return np.vstack([c, cl]), sg


@pytest.mark.parametrize("angle,expect", [
    (180.0, True),    # ideal backside
    (165.0, True),    # comfortably inside the window
    (150.0, True),    # exactly at the pre-registered edge
    (149.0, False),   # just outside it
    (97.6, False),    # the MEASURED median from reactive docking's default arm
    (64.9, False),    # its worst — the leaving group blocks the sulfur
])
def test_sn2_angle_window(angle, expect):
    coords, sg = sn2_coords(angle)
    r = nac.measure("sn2_displacement", coords, sg)
    assert r.viable is expect
    assert r.angle == pytest.approx(angle, abs=0.1)
    assert r.angle_kind == "S-C-LG"


def test_sn2_bond_distance_is_rejected_despite_perfect_angle():
    """The claim this module is built on.

    1.55 A is where reactive docking's published parameterisation puts the
    warhead. Even with a flawless 180 deg approach it is not a near-attack
    conformation of the FREE form -- it is a formed bond, and the free molecule
    cannot be there while the leaving group is attached.
    """
    coords, sg = sn2_coords(180.0, dist=1.55)
    r = nac.measure("sn2_displacement", coords, sg)
    assert r.angle == pytest.approx(180.0, abs=0.1)
    assert not r.viable, "a bonded distance must not pass as near-attack"


def test_sn2_too_far_is_rejected():
    coords, sg = sn2_coords(180.0, dist=6.0)
    assert not nac.measure("sn2_displacement", coords, sg).viable


def test_sn2_ring_opening_uses_the_same_geometry():
    coords, sg = sn2_coords(175.0)
    assert nac.measure("sn2_ring_opening", coords, sg).viable


# --------------------------------------------------------------------------
# Michael / SNAr — approach along the normal to an sp2 plane
# --------------------------------------------------------------------------

def planar_coords(off_normal_deg: float, dist: float = 3.5, flip: bool = False,
                  toward: bool = False):
    """Beta-C at the origin with its two partners in the xy-plane.

    The plane normal is +z, so `off_normal_deg` is measured from +z. The tilt is
    AWAY from `p1` by default, which is the Burgi-Dunitz direction — the
    nucleophile leans away from the adjacent substituent, giving an approach
    angle of 90 + off_normal. `toward=True` tilts the wrong way, which is the
    case the second constraint exists to reject.
    """
    centre = np.zeros(3)
    p1 = np.array([1.33, 0.0, 0.0])
    p2 = np.array([-0.75, 1.30, 0.0])
    t = np.radians(off_normal_deg)
    v = np.array([np.sin(t) * (1 if toward else -1), 0.0, np.cos(t)])
    if flip:
        v[2] *= -1
    return np.vstack([centre, p1, p2]), centre + dist * v


@pytest.mark.parametrize("off,expect", [
    (0.0, True),      # straight down the normal; approach 90 deg
    (17.0, True),     # the Burgi-Dunitz ideal: approach 107 deg
    (30.0, True),     # the pre-registered edge
    (31.0, False),    # just outside the off-normal cone
    (90.0, False),    # in the plane of the alkene — no orbital overlap
])
def test_michael_off_normal_window(off, expect):
    coords, sg = planar_coords(off)
    r = nac.measure("michael_addition", coords, sg)
    assert r.viable is expect
    assert r.angle == pytest.approx(off, abs=0.1)
    assert r.angle_kind == "off-normal"


def test_burgi_dunitz_ideal_is_measured_as_such():
    """The 107 deg approach must fall inside the window, not merely near it."""
    coords, sg = planar_coords(17.0)
    r = nac.measure("michael_addition", coords, sg)
    assert r.approach_angle == pytest.approx(107.0, abs=0.5)
    assert r.viable


def test_perpendicular_but_tilted_the_wrong_way_is_rejected():
    """The constraint that off-normal alone could not express.

    Tilting TOWARD the adjacent substituent is as far out of the plane as
    tilting away from it, so the off-normal test cannot tell them apart — but
    only one is a Burgi-Dunitz trajectory. Without this, a criterion that is
    geometrically perpendicular admits chemically wrong approaches, which is
    why HTS negatives cleared the old window 57% of the time.
    """
    away, sg_away = planar_coords(25.0)
    toward, sg_toward = planar_coords(25.0, toward=True)
    r_away = nac.measure("michael_addition", away, sg_away)
    r_toward = nac.measure("michael_addition", toward, sg_toward)
    assert r_away.angle == pytest.approx(r_toward.angle, abs=0.1), \
        "off-normal cannot distinguish these — that is the point"
    assert r_away.viable and not r_toward.viable


def test_sn2_carries_no_approach_angle():
    """None, not 0.0 — an unmeasured quantity must not read as a measured zero."""
    r = nac.measure("sn2_displacement", *sn2_coords(180.0))
    assert r.approach_angle is None


def test_either_face_of_the_plane_is_equally_valid():
    """An sp2 centre has two faces and a nucleophile may use either.

    Without folding onto [0, 90] the opposite face reads as 170 deg and is
    rejected, which would discard half of all genuinely viable Michael poses.
    """
    top, sg_top = planar_coords(10.0)
    bot, sg_bot = planar_coords(10.0, flip=True)
    rt = nac.measure("michael_addition", top, sg_top)
    rb = nac.measure("michael_addition", bot, sg_bot)
    assert rt.viable and rb.viable
    assert rt.angle == pytest.approx(rb.angle, abs=0.1)


def test_snar_is_perpendicular_not_anti():
    """SNAr and SN2 disagree about the same angle, which is why mechanism matters.

    In SNAr the nucleophile adds along the ring normal to form the Meisenheimer
    complex while the leaving group stays in the ring plane -- so S-C-LG is near
    90 deg, exactly the value that is DEAD for SN2. A single distance-plus-angle
    rule applied across mechanisms would invert the verdict for one of them.
    """
    coords, sg = planar_coords(0.0)
    assert nac.measure("snar_displacement", coords, sg).viable
    # the same arrangement judged as SN2: sulfur is 90 deg from the C-LG bond
    r_sn2 = nac.measure("sn2_displacement", coords[:2], sg)
    assert r_sn2.angle == pytest.approx(90.0, abs=0.1)
    assert not r_sn2.viable


# --------------------------------------------------------------------------
# failures are raised, never scored as "not viable"
# --------------------------------------------------------------------------

def test_unknown_mechanism_raises():
    coords, sg = sn2_coords(180.0)
    with pytest.raises(nac.NACError, match="unknown mechanism"):
        nac.measure("hydrolysis", coords, sg)


def test_missing_leaving_group_raises_rather_than_failing_the_molecule():
    with pytest.raises(nac.NACError, match="leaving-group"):
        nac.measure("sn2_displacement", np.zeros((1, 3)), np.array([3.5, 0, 0]))


def test_collinear_plane_raises():
    coords = np.array([[0., 0, 0], [1.3, 0, 0], [2.6, 0, 0]])
    with pytest.raises(nac.NACError, match="collinear"):
        nac.measure("michael_addition", coords, np.array([0., 0, 3.5]))


def test_bad_coord_shape_raises():
    with pytest.raises(nac.NACError, match=r"\(n, 3\)"):
        nac.measure("sn2_displacement", np.zeros(3), np.zeros(3))


# --------------------------------------------------------------------------

def test_viable_fraction():
    good = nac.measure("sn2_displacement", *sn2_coords(180.0))
    bad = nac.measure("sn2_displacement", *sn2_coords(90.0))
    assert nac.viable_fraction([good, bad, bad, bad]) == 0.25


def test_viable_fraction_refuses_an_empty_ensemble():
    """Unmeasured is not the same as unreactive."""
    with pytest.raises(nac.NACError):
        nac.viable_fraction([])


# --------------------------------------------------------------------------
# the isotropic baseline — what makes mechanisms comparable
# --------------------------------------------------------------------------

def test_sn2_null_matches_the_closed_form():
    """A single ideal direction, so a single cone — checked against the algebra."""
    assert nac.isotropic_null("sn2_displacement") == pytest.approx(
        (1 - np.cos(np.radians(30.0))) / 2, abs=1e-9)


def test_perpendicular_null_is_bounded_by_the_cone_it_lives_inside():
    """The quadrature has no closed form, so it is bracketed instead.

    The criterion is the INTERSECTION of an off-normal cone with a
    Burgi-Dunitz window, so its measure must be strictly less than the cone's
    alone — and non-zero, since the ideal trajectory satisfies both.
    """
    cone = 1 - np.cos(np.radians(nac.PERPENDICULAR_MAX_OFF_NORMAL))
    null = nac.isotropic_null("michael_addition")
    assert 0 < null < cone, f"{null:.4f} should sit strictly inside {cone:.4f}"


def test_the_two_mechanisms_are_now_comparably_selective():
    """The confound this baseline exists to remove, and the fix for it.

    MEASURED with the incomplete criterion: warhead-matched HTS inactives
    scored 3-19% viable under SN2 and 71-81% under Michael on the same
    receptor, because the perpendicular window was 4.4x wider by solid angle
    alone. Adding the missing Burgi-Dunitz constraint brings the two windows to
    within ~1.2x of each other.

    `enrichment()` divides the baseline out either way — this asserts that the
    criteria are now of comparable STRINGENCY, so the raw fractions are not
    wildly misleading to a reader who skips it.
    """
    ratio = nac.isotropic_null("michael_addition") / nac.isotropic_null("sn2_displacement")
    assert 0.5 < ratio < 2.0, f"windows differ in selectivity by {ratio:.1f}x"


def test_enrichment_is_one_when_the_molecule_does_no_better_than_chance():
    null = nac.isotropic_null("sn2_displacement")
    n = 1000
    good = nac.measure("sn2_displacement", *sn2_coords(180.0))
    bad = nac.measure("sn2_displacement", *sn2_coords(90.0))
    results = [good] * round(null * n) + [bad] * (n - round(null * n))
    assert nac.enrichment(results) == pytest.approx(1.0, abs=0.02)


def test_enrichment_refuses_to_mix_mechanisms():
    """Pooling mechanisms would divide by the wrong baseline and invent a result."""
    a = nac.measure("sn2_displacement", *sn2_coords(180.0))
    b = nac.measure("michael_addition", *planar_coords(0.0))
    with pytest.raises(nac.NACError, match="one mechanism"):
        nac.enrichment([a, b])
