"""
Purpose: is a docked pose a near-attack conformation? Mechanism-specific, not a distance cutoff.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: a pose (RDKit conformer) + the warhead's mechanism + the Cys113 SG position
Output: the measured approach geometry, and whether it clears a pre-registered window

STAGE 3 OF `docs/ranking_rationale.md`. The rationale states the requirement and
this module is where it becomes a number.

WHY THIS MODULE HAD TO EXIST — the measurement that forced it. Reactive docking
(D0063) biases sampling with a modified vdW pair potential between the ligand's
electrophilic atom and Cys113's sulfur. It works, in the narrow sense that it
does what it says: 20/20 poses placed the warhead carbon 1.55 A from the sulfur,
against a 1.80 A target.

Every single one of those poses was chemically dead. The S-C-Cl angle had median
97.6 deg and maximum 128.9 deg; SN2 at an sp3 carbon requires the nucleophile
ANTI to the leaving group, near 180 deg. In 8 of 20, the chlorine sat physically
between the sulfur and the carbon it was supposed to be attacking.

That is not a bug in reactive docking. The pair potential is a function of
DISTANCE ONLY, so it is isotropic by construction: it rewards the reactive atoms
being close and is indifferent to the direction of approach. A distance restraint
cannot encode an angle, and no tuning of it will.

WE TESTED WHETHER STERICS COULD DO IT INSTEAD, and they cannot, quite. The
published parameterisation deflates 1-3 and 1-4 neighbour radii by 0.5x so atoms
next to the reactive centre can interpenetrate -- that is what makes a 1.55 A
approach geometrically possible, and it is also what lets the leaving group
occupy the attack vector for free. Re-parameterising to a genuine near-attack
distance (3.2 A) with FULL steric radii admitted viable backside poses for the
first time -- 3/40 above 150 deg, versus 0/40 -- but the median was still 86 deg.
Real sterics help and do not suffice.

    CONCLUSION, recorded as D0064: the reactive potential is a SAMPLER, not a
    criterion. It puts poses in the right REGION. Whether a pose is
    reaction-competent is a separate question, answered here, on the poses it
    generates.

THE WINDOWS ARE PRE-REGISTERED. Per D0045 and the rationale's pre-registration
clause, every window below is fixed from textbook stereoelectronics BEFORE any
candidate is scored, never fitted to our own actives. A window tuned until the
ranking looks reasonable is not a criterion, it is a knob. Callers get the raw
measured angle alongside the verdict precisely so the window can be audited
rather than trusted.

IDENTITY, NOT POSITION. Atoms are located by SMARTS match on the molecule meeko
rebuilds from the docking output, which preserves the original atom ordering and
bonding. This module never indexes a pose by atom name or serial: a PDBQT names
every carbon "C", and taking a value by name where names are not unique is the
defect this project keeps rediscovering.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pre-registered windows. Fixed 2026-08-05, before any candidate was scored.
# ---------------------------------------------------------------------------

# Van der Waals contact between S (1.80 A) and C (1.70 A) is 3.50 A. A
# near-attack conformation in Bruice's sense is a contact, not a bond: the
# reactant state, with the bond not yet formed. Anything at 1.8 A is already
# past the transition state and is a geometry the FREE form cannot adopt while
# the leaving group is still attached.
NAC_DIST_MIN = 2.8
NAC_DIST_MAX = 4.2

# SN2 at sp3 carbon: backside attack, collinear with the C-LG bond. 180 deg is
# ideal; 150 deg is the generous end of what is usually accepted as near-attack.
# Deliberately generous -- a criterion that is too strict cannot be
# distinguished from one that is broken.
SN2_ANGLE_MIN = 150.0

# Conjugate addition and SNAr both proceed by attack roughly along the normal to
# a planar sp2 system -- the alkene plane for Michael, the ring plane for SNAr.
# 45 deg off-normal is generous; true Burgi-Dunitz-like trajectories are nearer
# 15-30 deg.
PERPENDICULAR_MAX_OFF_NORMAL = 45.0

MECHANISMS = {
    "sn2_displacement": "anti_to_leaving_group",
    "sn2_ring_opening": "anti_to_leaving_group",
    "michael_addition": "perpendicular_to_plane",
    "snar_displacement": "perpendicular_to_plane",
}


class NACError(ValueError):
    """The pose could not be measured — never silently treated as a failure."""


@dataclass(frozen=True)
class NACResult:
    """One pose's approach geometry.

    `viable` is the verdict, but `distance` and `angle` are always populated so
    a window can be re-drawn after the fact without re-docking. `angle_kind`
    names what the angle MEANS, because 90 deg is dead for SN2 and ideal for
    SNAr — reporting a bare number across mechanisms would be meaningless.
    """

    distance: float
    angle: float
    angle_kind: str
    mechanism: str
    viable: bool

    def __str__(self) -> str:  # pragma: no cover - display only
        return (f"{self.mechanism}: d={self.distance:.2f} A, "
                f"{self.angle_kind}={self.angle:.1f} deg, "
                f"{'VIABLE' if self.viable else 'not viable'}")


# ---------------------------------------------------------------------------
# geometry primitives
# ---------------------------------------------------------------------------

def _angle(a: np.ndarray, vertex: np.ndarray, c: np.ndarray) -> float:
    """Interior angle a-vertex-c in degrees."""
    u, v = a - vertex, c - vertex
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0 or nv == 0:
        raise NACError("degenerate angle: coincident atoms")
    return float(np.degrees(np.arccos(np.clip(u @ v / (nu * nv), -1.0, 1.0))))


def _off_normal(sg: np.ndarray, centre: np.ndarray,
                p1: np.ndarray, p2: np.ndarray) -> float:
    """Angle between the approach vector and the plane normal, in degrees.

    The plane is defined by `centre`, `p1`, `p2`. Folded onto [0, 90]: a
    nucleophile approaching from either face of an sp2 centre is equally valid,
    so 170 deg off one normal is a 10 deg approach off the other.
    """
    n = np.cross(p1 - centre, p2 - centre)
    nn = np.linalg.norm(n)
    if nn < 1e-6:
        raise NACError("degenerate plane: the three atoms are collinear")
    approach = sg - centre
    na = np.linalg.norm(approach)
    if na == 0:
        raise NACError("degenerate approach: sulfur coincides with the centre")
    off = np.degrees(np.arccos(np.clip(abs((n / nn) @ (approach / na)), 0.0, 1.0)))
    return float(off)


# ---------------------------------------------------------------------------

def measure(mechanism: str, coords: np.ndarray, sg: np.ndarray) -> NACResult:
    """Approach geometry for one pose.

    `coords` holds the atoms matched by the class's `reactive_atom_smarts`, in
    SMARTS order. Across every mechanism in the library, index 0 is the atom the
    sulfur attacks:

        sn2_displacement    [CH2][Cl]                 0=electrophilic C, 1=LG
        sn2_ring_opening    [CX3]([Br])=[NX2]         0=electrophilic C, 1=LG
        michael_addition    [CX3]=[CX3][CX3]=O        0=beta C, 1=alpha C, 2=C(=O)
        snar_displacement   [c]([Cl])[n]              0=ipso c, 1=LG, 2=ring n

    Raises rather than returning a false verdict when the geometry cannot be
    measured. "Could not be measured" and "is not reaction-competent" are
    different facts, and a molecule that fails to type must not be silently
    ranked as unreactive.
    """
    if mechanism not in MECHANISMS:
        raise NACError(f"unknown mechanism {mechanism!r}; "
                       f"known: {sorted(MECHANISMS)}")
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise NACError(f"coords must be (n, 3), got {coords.shape}")

    kind = MECHANISMS[mechanism]
    target = coords[0]
    dist = float(np.linalg.norm(sg - target))
    in_range = NAC_DIST_MIN <= dist <= NAC_DIST_MAX

    if kind == "anti_to_leaving_group":
        if len(coords) < 2:
            raise NACError(f"{mechanism} needs the leaving-group atom")
        angle = _angle(sg, target, coords[1])
        viable = in_range and angle >= SN2_ANGLE_MIN
        return NACResult(dist, angle, "S-C-LG", mechanism, viable)

    if len(coords) < 3:
        raise NACError(f"{mechanism} needs three atoms to define the plane")
    angle = _off_normal(sg, target, coords[1], coords[2])
    viable = in_range and angle <= PERPENDICULAR_MAX_OFF_NORMAL
    return NACResult(dist, angle, "off-normal", mechanism, viable)


def measure_poses(mol, smarts_match: tuple[int, ...], mechanism: str,
                  sg: np.ndarray) -> list[NACResult]:
    """`measure` over every conformer on `mol`, in conformer order.

    `smarts_match` comes from a substructure search on THIS molecule, so the
    indices address the same atoms the conformers hold.
    """
    out = []
    for cid in range(mol.GetNumConformers()):
        pos = mol.GetConformer(cid).GetPositions()
        out.append(measure(mechanism, pos[list(smarts_match)], sg))
    return out


def isotropic_null(mechanism: str) -> float:
    """Fraction of approach directions that would pass the window BY CHANCE.

    The window is a cone, so its solid angle divided by the sphere's is the rate
    a sulfur arriving from a uniformly random direction would clear it. Computed
    exactly rather than sampled, for the reason `pose_selection_bench` gives:
    this is the number every result is quoted against, and sampling it would put
    noise in the denominator.

        SN2, >=150 deg           cone of half-angle 30 deg    (1-cos30)/2 = 6.7%
        perpendicular, <=45 deg  two cones of half-angle 45   1-cos45    = 29.3%

    WHY THIS EXISTS. Measured on 3IKD, warhead-matched HTS inactives scored
    3-19% viable under the SN2 criterion and 71-81% under the Michael one. That
    is not a statement about the molecules -- the perpendicular window is 4.4x
    more permissive by solid angle alone. A raw viable fraction is therefore
    comparable WITHIN a mechanism and meaningless ACROSS mechanisms, which is
    D0020's "rank within warhead class, not globally" arriving from a new
    direction.

    Dividing by this gives an enrichment over chance, which is comparable across
    mechanisms because the criterion's own permissiveness has been divided out.
    """
    if mechanism not in MECHANISMS:
        raise NACError(f"unknown mechanism {mechanism!r}")
    if MECHANISMS[mechanism] == "anti_to_leaving_group":
        half = np.radians(180.0 - SN2_ANGLE_MIN)
        return float((1.0 - np.cos(half)) / 2.0)
    # both faces of the sp2 plane are usable, so two cones
    return float(1.0 - np.cos(np.radians(PERPENDICULAR_MAX_OFF_NORMAL)))


def enrichment(results: list[NACResult]) -> float:
    """Viable fraction divided by what an isotropic approach would give.

    1.0 means the molecule reaches near-attack geometry no more often than a
    randomly-oriented approach would; above 1.0 means the pocket and the
    molecule's own shape actively steer the warhead into position, which is the
    thing worth ranking on.
    """
    if not results:
        raise NACError("no poses to summarise")
    mechs = {r.mechanism for r in results}
    if len(mechs) != 1:
        raise NACError(f"enrichment needs one mechanism, got {sorted(mechs)}")
    return viable_fraction(results) / isotropic_null(mechs.pop())


def viable_fraction(results: list[NACResult]) -> float:
    """Fraction of poses that are reaction-competent.

    THE RANKING SIGNAL this module exists to produce. Docking runs are
    independent searches, so this estimates how readily the molecule can reach a
    reaction-competent geometry at all -- which is the question
    `docs/ranking_rationale.md` says to rank on, and it is a property of the
    molecule rather than a re-reading of a score.

    Empty input is an error, not 0.0: a molecule with no poses has not been
    shown to be unreactive, it has not been measured.
    """
    if not results:
        raise NACError("no poses to summarise")
    return sum(r.viable for r in results) / len(results)
