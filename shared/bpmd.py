"""
Purpose: rank candidates by whether the warhead STAYS in position, not by how often docking put it there.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: a solvated complex (from gromacs_explicit.solvate) + the warhead atom index
Output: a per-pose stability score, and the per-replica spread behind it

WHY THIS REPLACES THE FREQUENCY SCORE.

`nac_criterion` scores a molecule by the FRACTION of docking runs that reached a
near-attack conformation. D0068 measured that this does not converge: at 10x the
search effort the same molecules fall from 2.91x to 0.96x, and the
crystallographic positives -- never selected on score -- fall identically. More
search finds lower-energy poses that are less reaction-competent, so **the
frequency is partly a property of how hard you looked.**

Stability is not. You hand this module a pose and it measures whether that pose
survives thermal motion. Whatever search produced it is irrelevant by then. That
is the whole reason @tt8804 specified binding-pose metadynamics in #14 before we
substituted a cheaper proxy, and the proxy is the part that broke.

    THE SCREEN STILL DOES ITS JOB. The 200-run screen's top 300 sit at 0.96x
    against 0.99x for known crystallographic binders and 0.76x for random
    inactives -- it CONCENTRATES active-like molecules and cannot ORDER them.
    So it stays, as the filter that makes this affordable, and this module does
    the ordering.

THE COLLECTIVE VARIABLE IS THE WARHEAD DISTANCE, NOT LIGAND RMSD.

Standard BPMD biases along heavy-atom RMSD to the starting pose, which asks "does
the whole molecule stay put". That is the wrong question here and would actively
mislead: a molecule whose scaffold drifts while the warhead stays locked on
Cys113 is a GOOD answer, and RMSD would punish it. D0062 measured the same thing
from the other side -- whole-molecule RMSD correlates only rho=0.43 with reactive-
region placement, and 52.7% of poses get the reactive region right while failing
the RMSD test.

So the CV is **d(warhead reactive atom, Cys113 SG)**, and the score is how hard
the metadynamics has to push to get the warhead out of the near-attack window.

WELL-TEMPERED, AND WHY. Plain metadynamics keeps depositing bias forever and the
free-energy estimate never settles; well-tempered metadynamics decays the hill
height so the estimate converges. Given D0068, a method whose answer depends on
how long you ran it is exactly what we cannot afford twice.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

# --------------------------------------------------------------------------
# Protocol constants. The replica count and length are DELIBERATELY not fixed
# here -- `scripts/bpmd_convergence.py` measures what is sufficient on a known
# active before any production run commits to them. Fixing them by assumption is
# how D0068 happened.
# --------------------------------------------------------------------------

# The near-attack window from nac_criterion: a van der Waals contact, not a bond.
NAC_MIN_NM = 0.28
NAC_MAX_NM = 0.42

# Hill deposition. HEIGHT is in kJ/mol; SIGMA is the CV's resolution in nm and
# should be roughly the thermal fluctuation of the distance in an unbiased run --
# too narrow wastes compute, too wide blurs the barrier being measured.
HILL_HEIGHT_KJ = 0.3
HILL_SIGMA_NM = 0.02
HILL_PACE = 500                    # steps between depositions
BIASFACTOR = 10.0                  # well-tempered: bias decays, estimate converges

# Beyond this the warhead has unambiguously left; no point simulating further.
UNBOUND_NM = 1.0

# --------------------------------------------------------------------------
# THE CV MUST BE BOUNDED BY A WALL, NOT BY THE COMMITTOR.
#
# The first convergence run died in every one of its replicas:
#
#     An error occurred while PLUMED was calculating the forces
#     (tools/Grid.cpp:111) PLMD::GridBase::index_t ... getIndex
#     An error happened while calculating metad
#
# That is METAD indexing its bias grid outside the grid's own bounds. The
# protocol relied on `COMMITTOR` to end the run once the warhead was gone, so
# the grid only had to reach a little past UNBOUND_NM and GRID_MAX was 2.0 nm.
#
# COMMITTOR FIRED AND THE RUN DID NOT STOP. PLUMED's own output records
# "SET COMMITTED TO BASIN 1" from t = 681 ps onward and keeps printing it for
# the next two nanoseconds: PLUMED raised its stop flag every step and GROMACS
# ignored it. The CV then wandered to 1.97 nm and the next step took it off the
# grid. Every replica crashed at ~3 ns of a 10 ns run, so the convergence run
# that was supposed to fix the protocol's parameters produced nothing at all.
#
# THE FAILURE LOOKED LIKE A TOPOLOGY PROBLEM. It arrived immediately after the
# "GPU update kernel rejected this atom ordering; retrying with -update cpu"
# warning, which is unrelated (GROMACS always refuses the GPU update path when
# PLUMED is driving forces) and reads like the cause.
#
# So the CV is now confined by an UPPER_WALLS restraint, which is a force and
# cannot be ignored by the host MD engine. The wall sits BEYOND UNBOUND_NM: by
# the time it is felt the warhead has already left by any definition this
# module uses, so it cannot alter the escape barrier being measured — it only
# stops an escaped ligand from diffusing off the grid. GRID_MAX keeps clear
# headroom above the wall so even a hard overshoot stays indexable.
WALL_NM = 1.5                      # 0.5 nm past "unambiguously gone"
WALL_KAPPA_KJ_PER_NM2 = 2000.0
GRID_MIN_NM = 0.1
GRID_MAX_NM = 2.5                  # 1.0 nm of headroom above the wall
GRID_SPACING_NM = 0.005            # ~4 bins per SIGMA


class BPMDError(RuntimeError):
    """The run could not be set up or its output could not be read."""


@dataclass(frozen=True)
class ReplicaResult:
    """One metadynamics replica."""

    replica: int
    max_cv_nm: float               # furthest the warhead got
    frac_in_window: float          # time fraction inside the near-attack window
    bias_at_exit_kj: float         # bias needed before it left, NaN if it never did
    escaped: bool


@dataclass(frozen=True)
class PoseStability:
    """A pose's stability, and the spread across replicas that produced it.

    `spread` is reported and never hidden. Replica-to-replica variation IS the
    uncertainty on this number, and a mean without it would repeat the mistake
    of quoting a single unqualified value.
    """

    n_replicas: int
    mean_frac_in_window: float
    spread_frac: float
    median_bias_to_escape_kj: float
    n_escaped: int

    @property
    def score(self) -> float:
        """Higher is more stable. NaN-free by construction.

        Two components deliberately kept multiplicative: a pose must BOTH stay in
        the window most of the time AND be expensive to push out. A pose that sits
        in the window only because nothing pushed it is not stable, and one that
        resists a large bias but spends most of its time outside was never in a
        near-attack conformation to begin with.
        """
        cost = (self.median_bias_to_escape_kj
                if np.isfinite(self.median_bias_to_escape_kj) else 0.0)
        return float(self.mean_frac_in_window * (1.0 + cost))


def plumed_input(warhead_idx: int, sg_idx: int, *, stride: int = 500,
                 hill_height: float = HILL_HEIGHT_KJ,
                 sigma: float = HILL_SIGMA_NM,
                 pace: int = HILL_PACE,
                 biasfactor: float = BIASFACTOR,
                 temp_k: float = 300.0) -> str:
    """A PLUMED input biasing the warhead-to-sulfur distance.

    Atom indices are **1-based**, which is PLUMED's convention and not Python's.
    Passing a 0-based index here biases the wrong atom and the run still
    completes and still produces a plausible number -- the exact failure shape
    this project keeps finding -- so the caller is required to have converted,
    and `bpmd_atom_indices` is the only supported way to get them.
    """
    if warhead_idx < 1 or sg_idx < 1:
        raise BPMDError(f"PLUMED indices are 1-based; got warhead={warhead_idx}, "
                        f"sg={sg_idx}. Use bpmd_atom_indices().")
    nbin = int(round((GRID_MAX_NM - GRID_MIN_NM) / GRID_SPACING_NM))
    return f"""# Binding-pose metadynamics, biased along the WARHEAD-SULFUR distance.
# Not ligand RMSD: a scaffold that drifts while the warhead stays locked is a
# good answer, and RMSD would punish it (D0062).
d: DISTANCE ATOMS={warhead_idx},{sg_idx}

metad: METAD ARG=d ...
   PACE={pace}
   HEIGHT={hill_height}
   SIGMA={sigma}
   BIASFACTOR={biasfactor}
   TEMP={temp_k}
   GRID_MIN={GRID_MIN_NM}
   GRID_MAX={GRID_MAX_NM}
   GRID_BIN={nbin}
   FILE=HILLS
...

# KEEPS THE CV ON THE GRID. This is a force, so the MD engine cannot decline it
# the way it declined COMMITTOR's stop flag. It sits {WALL_NM - UNBOUND_NM:.1f} nm beyond
# UNBOUND_NM, so it acts only on ligands that have already left.
uwall: UPPER_WALLS ARG=d AT={WALL_NM} KAPPA={WALL_KAPPA_KJ_PER_NM2}

# INFORMATIONAL ONLY — GROMACS DOES NOT HONOUR THIS. PLUMED raises its stop
# flag and the run continues; the convergence run's replicas all logged
# "SET COMMITTED TO BASIN 1" for 2 ns before dying off the end of the grid.
# Escape is therefore read from the COLVAR by `read_colvar`, never inferred
# from a run having stopped early, and cost is bounded by the run LENGTH.
COMMITTOR ARG=d BASIN_LL1={UNBOUND_NM} BASIN_UL1=99.0

PRINT ARG=d,metad.bias,uwall.bias FILE=COLVAR STRIDE={stride}
"""


def bpmd_atom_indices(gro_path: Path, warhead_serial: int,
                      sg_serial: int) -> tuple[int, int]:
    """Convert 0-based atom serials to PLUMED's 1-based indices, checking they exist.

    Exists so the conversion happens in exactly one place, with a bounds check.
    An off-by-one here biases a neighbouring atom and the simulation runs
    perfectly.
    """
    n = 0
    for i, line in enumerate(gro_path.read_text().splitlines()):
        if i == 1:
            n = int(line.strip())
            break
    for name, s in (("warhead", warhead_serial), ("SG", sg_serial)):
        if not 0 <= s < n:
            raise BPMDError(f"{name} serial {s} outside the system's {n} atoms")
    return warhead_serial + 1, sg_serial + 1


def read_colvar(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """(distance in nm, metadynamics bias in kJ/mol) from a PLUMED COLVAR file."""
    if not path.is_file():
        raise BPMDError(f"no COLVAR at {path}")
    d, b = [], []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            d.append(float(parts[1]))
            b.append(float(parts[2]))
        except ValueError:
            continue
    if not d:
        raise BPMDError(f"{path} contained no readable frames")
    return np.array(d), np.array(b)


def analyse_replica(colvar: Path, replica: int) -> ReplicaResult:
    """Reduce one replica's trajectory to the four numbers that matter."""
    d, bias = read_colvar(colvar)
    in_win = ((d >= NAC_MIN_NM) & (d <= NAC_MAX_NM)).mean()
    escaped = bool((d >= UNBOUND_NM).any())
    if escaped:
        # The bias standing when it first left. That is the depth of the well it
        # had to climb out of, which is the quantity being ranked.
        bias_at_exit = float(bias[np.argmax(d >= UNBOUND_NM)])
    else:
        bias_at_exit = float("nan")
    return ReplicaResult(replica=replica, max_cv_nm=float(d.max()),
                         frac_in_window=float(in_win),
                         bias_at_exit_kj=bias_at_exit, escaped=escaped)


def combine(results: list[ReplicaResult]) -> PoseStability:
    """Pool replicas into one stability, keeping the spread.

    A replica that never escaped contributes its **final** bias as a LOWER BOUND
    on the escape cost rather than being dropped. Dropping non-escapers would
    bias the median toward the weak poses -- the ones that left -- which is
    precisely backwards.
    """
    if not results:
        raise BPMDError("no replicas to combine")
    fracs = np.array([r.frac_in_window for r in results])
    costs = np.array([r.bias_at_exit_kj for r in results], dtype=float)
    finite = costs[np.isfinite(costs)]
    return PoseStability(
        n_replicas=len(results),
        mean_frac_in_window=float(fracs.mean()),
        spread_frac=float(fracs.std(ddof=1)) if len(fracs) > 1 else 0.0,
        median_bias_to_escape_kj=float(np.median(finite)) if finite.size else float("nan"),
        n_escaped=int(sum(r.escaped for r in results)),
    )
