"""
Purpose: how much do a molecule's top-scoring poses AGREE about where the reactive region sits?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: a molecule's docked poses — each an energy plus the coordinates of its reactive atoms
Output: a ConsensusResult: agreement fraction, the pairwise-RMSD distribution behind it, and the N it was measured at

COMPONENT 2 OF @tt8804's RANKING DESIGN (`docs/ranking_rationale.md`).
Component 1 is pose quality (`shared/nac_criterion.py`), 3 is BPMD and 4 is MD
residence. This is the cheap one: it re-reads poses a docking run has already
produced and starts no new simulation.

WHY IT EXISTS — D0068. The viable-NAC FRACTION does not converge. The same 300
molecules scored 2.91x median enrichment at 200 runs and 0.96x at 2,000, and the
15 crystallographic positives — never selected on score — fell identically
(2.27x -> 0.97x), so it is not winner's curse. More search finds lower-energy
poses that are less reaction-competent, which means the viable fraction is a
property of the SEARCH as much as of the molecule.

    A frequency is a ratio over the whole pose population, so it moves whenever
    the population's composition moves, and more search always moves it.

The hypothesis this module exists to test is that AGREEMENT does not inherit the
search the same way. It is measured over the top-N poses by energy — a
fixed-size window at the TOP of the ranking, not the whole population. More
search can only replace a member of that window with a lower-energy one, so if
the docking has a genuine minimum the window converges on it. Nothing in that
argument is proved here. See "How this is falsified" below.

WHAT IS COMPARED, AND WHY NOT WHOLE-MOLECULE RMSD. D0062 measured
spearman(whole-molecule RMSD, error in reactive-region placement) = +0.433 over
738 poses on 3IKD, and 52.7% of all poses place the reactive region correctly
while FAILING the 2 A RMSD test. A pose whose scaffold has swung round while the
warhead stayed put is a good answer to a covalent question, and whole-molecule
RMSD scores it as a failure. So the comparison here is over the reactive atoms
only — the same atoms `nac_criterion.measure` consumes, from the same
`reactive_atom_smarts` match, in the same SMARTS order.

NO SUPERPOSITION. The RMSD is computed in the receptor frame as the poses were
docked. Superposing two poses before comparing them would ask "is the warhead in
the same place RELATIVE TO ITS OWN MOLECULE", which is a question about the
ligand's internal conformation and is answered identically by two poses at
opposite ends of the pocket. The question here is whether the reactive region
lands in the same place IN THE SITE.

NOTHING IS EVER AVERAGED. `representative_index` is an INDEX into the poses the
caller passed in, obtained from `pose_vector.representative`, so it is always a
conformation some program actually generated. The mean of a bimodal pose set
lies in the trough where no pose ever was, and its geometry is a plausible
number computed from a structure that does not exist — the failure shape
`docs/how_this_project_breaks.md` is about.

THE SCORE SELECTS THE POSES AND DOES NOT ENTER THE NUMBER. Energy decides WHICH
poses are compared; it contributes no weight to the agreement. That keeps this
from being a re-reading of the docking score, which five levels of theory have
failed on (D0041, D0046, D0036, D0038/D0044, D0057) — though note the honest
caveat that top-N selection does import the score's ordering, so consensus is
not fully independent of it either.

    THE POSES MUST COME FROM INDEPENDENT RUNS. AutoDock-GPU reseeds every
    invocation, so a 200-run screen gives 200 independent searches and agreement
    among them means something (D0065; `n_poses == nrun` exactly, D0068). The
    nine modes of ONE Vina invocation are diversified by a minimum-RMSD floor —
    `pose_vector.cluster` carries the same warning — so a consensus computed
    over them measures Vina's output filter and will read artificially LOW.
    Nothing in this module can detect that from coordinates. It is the caller's
    responsibility.

WHY THE HEADLINE IS PAIRWISE AND NOT "IS THE BIGGEST CLUSTER BIG".
A largest-cluster fraction reports the size of the winning mode and says nothing
about how many rivals it beat: 5 poses in one tight cluster and 5 in another,
6 A away, gives a perfectly respectable 0.5 and hides that the molecule has two
incompatible answers. `agreement` counts pose PAIRS that land within tolerance
of each other, so the cross-cluster pairs — the majority in a balanced bimodal
set — count against it, and that case reads 0.44 rather than "half the poses
agree". `n_modes` and `dominant_mode_fraction` are reported alongside so the
shape is visible, never as the headline.

NEVER A BARE CENTRAL VALUE. `agreement` always travels with its jackknife range
(the same number recomputed with each pose left out in turn), and the
pairwise-RMSD distribution always travels with its IQR and its maximum. D0065's
own lesson, arriving here: four seeds of the chloroacetamide arm gave AUC 0.872,
0.881, 0.852, 0.822 and the honest report is the spread, not the single run.

N IS PART OF THE METRIC'S DEFINITION, not a tuning knob. D0068 consequence 2
says any enrichment must be quoted with the run count that produced it; the same
applies here with more force, because N is chosen rather than inherited. So
`top_n` is a REQUIRED argument with no default, every result carries the N it
was actually measured at, and `require_same_n` refuses to let results measured
at different N be compared.

HOW THIS IS FALSIFIED. Dock the same molecules at 200 and at 2,000 runs and
compute consensus at matched N and matched tolerance. If it moves the way
enrichment moved — a median shift of the same size, or a rank correlation near
D0068's rho = 0.364 — then consensus inherits the search too and the premise
above is wrong. Second falsifier: consensus at N = 5 and at N = 50 on the same
run should agree in RANK even though they will differ in value; if they do not,
the number is a property of the window rather than of the molecule.

NOT A RANK METRIC YET, AND HIGHER IS BETTER. `rank_shortlist.LOWER_IS_BETTER` is
the registry that catches direction mix-ups, and it only knows lower-is-better
columns; `agreement` must not be passed to `rank_by` as if it were one. And per
#13 nothing here labels anything downstream until it has been scored against the
covalent validation set.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from shared import pose_vector as pv

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pre-registered constants.
# ---------------------------------------------------------------------------

# Two poses "agree" when their reactive regions land within this of each other.
# NOT INVENTED: D0062 measured reactive-region placement against crystal poses
# at a 1.0 A radius (59.8% of poses inside it), so this is the tolerance at
# which "the reactive region is in the same place" has already been given a
# meaning in this project. Reused rather than re-chosen, so agreement here and
# placement accuracy there are talking about the same distance.
#
# It doubles as the single-linkage clustering threshold, so there is one knob
# rather than two that can drift apart.
DEFAULT_TOLERANCE_A = 1.0

# Below three poses there is nothing to be said about spread. Two poses give
# exactly one pairwise distance: no IQR, no jackknife, and no way to tell a
# genuine agreement from a coincidence of two. A molecule with fewer than this
# has NOT BEEN MEASURED, which is a different fact from "its poses disagree",
# and the two must not collide on the same number.
MIN_POSES_FOR_CONSENSUS = 3

# Tolerance for the <= comparison at the boundary. A pose pair built to sit
# exactly on the tolerance is a chemical statement, and which side of the last
# ulp the arithmetic lands on must not decide it. Same reasoning, same value, as
# `nac_criterion._TOL`.
_TOL = 1e-9


class ConsensusError(ValueError):
    """Consensus could not be MEASURED — never silently reported as 0.0.

    Raised for too few poses, an unusable pose set, or an incoherent request.
    "Could not be measured" and "the poses disagree" are different facts about a
    molecule and a caller must be able to tell them apart: the first says look
    at the docking, the second says look at the molecule.
    """


@dataclass(frozen=True)
class ReactivePose:
    """One docked pose, reduced to the two things consensus needs.

    `energy` is the docking score in kcal/mol, LOWER IS BETTER. It is used only
    to choose which poses enter the comparison.

    `reactive_xyz` is (n_reactive, 3): the coordinates of the atoms matched by
    the warhead class's `reactive_atom_smarts`, in SMARTS order — the same
    array `nac_criterion.measure` is handed, so both components of the ranking
    read the same atoms.

    `atom_ids` are those atoms' indices in the molecule. They travel WITH the
    coordinates for the reason `PoseVector.resi` does: two pose sets built from
    different SMARTS matches are not comparable, and a bare coordinate array
    gives a caller no way to notice. A symmetric warhead matched at two
    different reactive centres produces two different `atom_ids` and will refuse
    to compare rather than averaging over the ambiguity — which is the same
    resolution D0065 reached by docking symmetric warheads once per centre.
    """

    energy: float
    reactive_xyz: np.ndarray
    atom_ids: tuple[int, ...]

    def __post_init__(self) -> None:
        xyz = np.asarray(self.reactive_xyz, dtype=float)
        if xyz.ndim != 2 or xyz.shape[1] != 3:
            raise ConsensusError(
                f"reactive_xyz must be (n_reactive, 3), got {xyz.shape}")
        if xyz.shape[0] == 0:
            raise ConsensusError("pose has no reactive atoms")
        if len(self.atom_ids) != xyz.shape[0]:
            raise ConsensusError(
                f"{len(self.atom_ids)} atom ids but {xyz.shape[0]} coordinates")
        if not np.isfinite(xyz).all():
            raise ConsensusError("reactive_xyz contains non-finite coordinates")
        if not np.isfinite(self.energy):
            raise ConsensusError(
                f"energy is {self.energy!r}; a pose with no energy cannot be "
                "placed in a top-N by energy")


@dataclass(frozen=True)
class ConsensusResult:
    """A molecule's pose agreement, and everything needed to audit it.

    `agreement` is the headline: the fraction of pose PAIRS among the top
    `n_poses` whose reactive regions land within `tolerance_a` of each other.
    1.0 is perfect agreement, 0.0 is none. It is never quoted alone —
    `agreement_jackknife` is (min, max) of the same quantity recomputed with
    each pose left out in turn, so a value that rests on one pose is visible.

    `median_rmsd` with `iqr_rmsd` and `max_rmsd` describe the pairwise
    reactive-region RMSD distribution the fraction was cut from, so the
    tolerance can be re-drawn after the fact without re-docking. This is the
    same reason `NACResult` carries its raw distance and angle beside its
    verdict.

    `n_modes` and `dominant_mode_fraction` come from single-linkage clustering
    at the same tolerance. They are DESCRIPTIVE. `dominant_mode_fraction` is a
    population count and is only meaningful when the poses came from
    independent runs — see the module docstring.

    `representative_index` indexes the ORIGINAL `poses` list, not the top-N
    slice, and is a real pose. When `n_modes > 1` it describes ONE mode and
    says nothing about the others; `agreement` is the field that says so.

    `top_n_requested` and `n_poses` are both recorded because they can differ:
    a molecule that produced fewer poses than requested is measured on what
    exists, and the number it gets is not comparable with one measured on more.
    `require_same_n` is the guard for that.
    """

    top_n_requested: int
    n_poses: int
    tolerance_a: float
    agreement: float
    agreement_jackknife: tuple[float, float]
    median_rmsd: float
    iqr_rmsd: tuple[float, float]
    max_rmsd: float
    n_modes: int
    dominant_mode_fraction: float
    representative_index: int

    @property
    def is_multimodal(self) -> bool:
        """More than one binding mode at `tolerance_a`.

        When True, no single pose and no single central value describes this
        molecule's answer — it has more than one.
        """
        return self.n_modes > 1

    def __str__(self) -> str:  # pragma: no cover - display only
        lo, hi = self.agreement_jackknife
        q1, q3 = self.iqr_rmsd
        short = "" if self.n_poses == self.top_n_requested else \
            f" (asked for {self.top_n_requested})"
        return (f"agreement={self.agreement:.3f} [{lo:.3f}, {hi:.3f}] "
                f"at N={self.n_poses}{short}, tol={self.tolerance_a:.2f} A; "
                f"pairwise RMSD median={self.median_rmsd:.2f} "
                f"IQR=[{q1:.2f}, {q3:.2f}] max={self.max_rmsd:.2f} A; "
                f"{self.n_modes} mode(s), "
                f"dominant={self.dominant_mode_fraction:.2f}")


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------

def _as_pose_vectors(poses: list[ReactivePose]) -> list[pv.PoseVector]:
    """Flatten reactive-atom coordinates so `pose_vector`'s machinery applies.

    A PoseVector here holds 3*n_reactive coordinates rather than a per-residue
    contact profile, and its basis is each atom id repeated three times. The
    Euclidean distance `pose_vector.distance_matrix` computes over that is
    exactly RMSD * sqrt(n_reactive), so clustering and medoid selection are the
    same operations on the same ordering — and `require_same_basis` becomes the
    check that two pose sets describe the same reactive atoms.

    PRIVATE ON PURPOSE. These are not contact profiles and must never reach
    `pose_vector.fit_score` or `reference_profile`, which would compare them
    against a residue-distance reference and return a plausible number from a
    meaningless comparison.
    """
    return [pv.PoseVector(tuple(a for a in p.atom_ids for _ in range(3)),
                          np.asarray(p.reactive_xyz, dtype=float).ravel())
            for p in poses]


def reactive_rmsd(a: ReactivePose, b: ReactivePose) -> float:
    """RMSD between two poses' reactive regions, in the receptor frame.

    No superposition and no symmetry perception: see the module docstring for
    both. Raises if the two poses describe different atoms.
    """
    va, vb = _as_pose_vectors([a, b])
    if va.resi != vb.resi:
        raise ConsensusError(
            "these poses' reactive regions are different atoms "
            f"({a.atom_ids} vs {b.atom_ids}); they cannot be compared")
    n = len(a.atom_ids)
    return float(np.linalg.norm(va.values - vb.values) / np.sqrt(n))


def rmsd_matrix(poses: list[ReactivePose]) -> np.ndarray:
    """Pairwise reactive-region RMSD, (n, n)."""
    if not poses:
        raise ConsensusError("no poses")
    vectors = _as_pose_vectors(poses)
    try:
        pv.require_same_basis(vectors)
    except ValueError as exc:
        raise ConsensusError(
            "these poses' reactive regions are not the same atoms; a consensus "
            "over them would compare different parts of different molecules"
        ) from exc
    n_atoms = len(poses[0].atom_ids)
    return pv.distance_matrix(vectors) / np.sqrt(n_atoms)


# ---------------------------------------------------------------------------

def _agreement(d: np.ndarray, tolerance_a: float) -> float:
    """Fraction of distinct pose pairs within tolerance."""
    iu = np.triu_indices(d.shape[0], k=1)
    return float(np.mean(d[iu] <= tolerance_a + _TOL))


def consensus(poses: list[ReactivePose], *, top_n: int,
              tolerance_a: float = DEFAULT_TOLERANCE_A) -> ConsensusResult:
    """Agreement among a molecule's `top_n` lowest-energy poses.

    `top_n` is REQUIRED and has no default. A consensus is not a property of a
    molecule alone — it is a property of a molecule at a stated N — and D0068 is
    what happens when a number's defining parameter is left implicit.

    Poses are ordered by energy ascending, ties broken by their position in
    `poses` so the selection is deterministic and re-running cannot silently
    swap a tied pair.

    Raises `ConsensusError` when fewer than `MIN_POSES_FOR_CONSENSUS` poses are
    available. A molecule that could not be measured must not be handed back a
    0.0 that reads as "its poses disagree".

    When fewer than `top_n` poses exist but at least the minimum do, the
    consensus is measured on what exists and `n_poses` records the difference.
    Such a result is not comparable with one measured at the full N — see
    `require_same_n`.
    """
    if top_n < MIN_POSES_FOR_CONSENSUS:
        raise ConsensusError(
            f"top_n={top_n} is below {MIN_POSES_FOR_CONSENSUS}; with fewer "
            "poses there is one pairwise distance at most and no spread to "
            "report, so the result would be a central value with nothing "
            "beside it")
    if tolerance_a <= 0:
        raise ConsensusError(
            f"tolerance_a={tolerance_a} must be positive; a zero tolerance "
            "makes every pose its own mode by construction")
    if len(poses) < MIN_POSES_FOR_CONSENSUS:
        raise ConsensusError(
            f"{len(poses)} pose(s) is too few to measure consensus (need "
            f"{MIN_POSES_FOR_CONSENSUS}). This molecule has NOT been measured; "
            "it has not been shown to have disagreeing poses")

    order = sorted(range(len(poses)), key=lambda i: (poses[i].energy, i))
    if len(poses) < top_n:
        log.warning(
            "asked for the top %d poses but only %d exist; measuring at N=%d. "
            "This value is not comparable with one measured at N=%d.",
            top_n, len(poses), len(poses), top_n)
    chosen = order[:top_n]
    selected = [poses[i] for i in chosen]
    n = len(selected)

    d = rmsd_matrix(selected)
    iu = np.triu_indices(n, k=1)
    pair = d[iu]

    # Jackknife on the SELECTED set: how far does the headline move if any one
    # of these poses had not been found? Deliberately not a re-selection of the
    # top-N from a smaller pool, which would answer a different question.
    jack = []
    for i in range(n):
        keep = [j for j in range(n) if j != i]
        jack.append(_agreement(d[np.ix_(keep, keep)], tolerance_a))

    n_atoms = len(selected[0].atom_ids)
    labels = pv.cluster(_as_pose_vectors(selected),
                        threshold=tolerance_a * np.sqrt(n_atoms) + _TOL)
    sizes = np.bincount(labels)

    q1, q3 = (float(x) for x in np.percentile(pair, [25.0, 75.0]))
    return ConsensusResult(
        top_n_requested=top_n,
        n_poses=n,
        tolerance_a=float(tolerance_a),
        agreement=_agreement(d, tolerance_a),
        agreement_jackknife=(float(min(jack)), float(max(jack))),
        median_rmsd=float(np.median(pair)),
        iqr_rmsd=(q1, q3),
        max_rmsd=float(pair.max()),
        n_modes=int(len(sizes)),
        dominant_mode_fraction=float(sizes.max() / n),
        # Map the medoid back onto the caller's own indexing, so what comes out
        # addresses a pose they hold rather than a position in a slice this
        # function made.
        representative_index=int(chosen[pv.representative(
            _as_pose_vectors(selected))]),
    )


def require_same_n(results: list[ConsensusResult]) -> tuple[int, float]:
    """Refuse to compare consensus values measured at different N or tolerance.

    D0068 consequence 2 — "`nrun` becomes part of the metric's definition, not
    a tuning knob" — applied to this metric's own defining parameters. Two
    molecules whose consensus was measured over 10 poses and over 40 are not
    ranked against each other by this module's numbers, and a caller ranking a
    column without checking has no way to notice: both values are populated,
    plausible, and in [0, 1].

    Returns the (N, tolerance) they share.
    """
    if not results:
        raise ConsensusError("no results")
    ns = {r.n_poses for r in results}
    tols = {r.tolerance_a for r in results}
    if len(ns) != 1:
        raise ConsensusError(
            f"these consensus values were measured at different N ({sorted(ns)}) "
            "and are not comparable; re-measure them at one N")
    if len(tols) != 1:
        raise ConsensusError(
            f"these consensus values used different tolerances ({sorted(tols)}) "
            "and are not comparable")
    return ns.pop(), tols.pop()
