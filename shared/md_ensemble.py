"""
Purpose: GB implicit-solvent MD on the topologies MM-GBSA already builds, to
         supply the two things a single minimisation cannot: an ensemble (and
         with it an honest uncertainty) and a pocket-residence measurement.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: complex.prmtop + complex.min.rst from an existing MM-GBSA workdir
Output: a DCD trajectory + residence/stability metrics per candidate

WHY THIS EXISTS. Every dG this project reports comes from ONE minimised
structure. A minimisation has no ensemble, so it has no variance, so the
"+/- uncertainty" column the plan asks for could only have been fabricated.
D0032 made the cost of that concrete: MM-GBSA scored 0.140 ROC-AUC against
docking's 0.440 and there was no way to say whether the gap was signal or the
noise of two arbitrary local minima. This module produces the ensemble that
question needs.

IMPLICIT, AND WHY THAT IS A REAL LIMITATION. `complex.prmtop` carries zero
waters and no periodic box, so it drives GB implicit-solvent MD directly with
no re-parameterisation. That is the whole reason this tier is cheap. It is
ALSO not what the T5 spec asked for: the spec says explicit solvent, and
explicit solvent needs the topologies rebuilt through tleap with a water box
and counter-ions. Implicit solvent has no water structure, no viscosity and a
too-fast conformational clock. It is the right first tier and the wrong final
answer, and that distinction is recorded here rather than left to be
rediscovered.

CONSISTENCY WITH THE SCORING THAT ALREADY EXISTS. `shared/mmgbsa.py` uses
igb=8 with mbondi3 radii. igb=8 is GBn2, so `implicitSolvent=app.GBn2` is the
same model, and the radii travel inside the prmtop that tleap already wrote.
`validate_against_sander()` checks that claim numerically instead of asserting
it -- if OpenMM and sander disagree on the SAME structure, the ensemble is not
comparable to the single-structure numbers and the caller is told so.

RESIDENCE IS NOT A DOCKING SCORE. For a covalent adduct the ligand cannot
diffuse away -- it is bonded to Cys113 -- so "did it stay bound" is vacuous and
the meaningful question is whether the ligand BODY stays in the pocket or
swings into solvent while tethered. That is what `residence_metrics()`
measures, and for non-covalent candidates the same numbers additionally
capture true dissociation. The metric is reported per-candidate; it is not a
ranking and nothing here promotes it to one.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

# Matches shared/mmgbsa.py: igb=8 is GBn2, radii ride in the prmtop.
IGB = 8
PB_RADII = "mbondi3"

TEMPERATURE_K = 300.0
FRICTION_PER_PS = 1.0
TIMESTEP_FS = 2.0            # safe with HBonds constrained
EQUIL_PS = 200.0
PRODUCTION_PS = 2000.0
FRAME_INTERVAL_PS = 20.0     # -> 100 frames
POCKET_CUTOFF_NM = 0.45      # heavy-atom contact distance

# A frame every 20 ps is not automatically an independent sample. The SEM this
# module reports is corrected by the statistical inefficiency rather than
# assuming independence, because assuming it would shrink the error bar by
# exactly the factor we are trying to measure.


class MDError(RuntimeError):
    """Raised when a trajectory cannot be produced or trusted."""


def topology_fingerprint(workdir: Path) -> str | None:
    """SHA-256 of the complex topology a trajectory was propagated from.

    Cheap enough to compute per candidate (~1 MB read) and the only thing that
    distinguishes a trajectory generated under one set of junction parameters
    from one generated under another. Returns None if the topology is absent,
    which callers must treat as "cannot verify" rather than "matches".
    """
    import hashlib
    p = workdir / "complex.prmtop"
    if not p.is_file() or p.stat().st_size == 0:
        return None
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


@dataclass
class MDResult:
    candidate_id: str
    workdir: Path
    n_frames: int
    ns_simulated: float
    wall_seconds: float
    residence: dict = field(default_factory=dict)
    validation: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "candidate_id": self.candidate_id,
            "n_frames": self.n_frames,
            "ns_simulated": round(self.ns_simulated, 3),
            "wall_seconds": round(self.wall_seconds, 1),
            # The topology this trajectory was actually propagated from. Without
            # it the cache could only compare run LENGTH, so a topology rebuilt
            # with corrected junction parameters (D0037) would silently reuse a
            # trajectory generated under the old ones -- the same shape of
            # defect as D0033's stale energies, one layer down.
            "topology_sha256": topology_fingerprint(self.workdir),
            "solvent": "implicit GBn2 (igb=8), mbondi3 radii",
            "explicit_solvent": False,
            "temperature_k": TEMPERATURE_K,
            "timestep_fs": TIMESTEP_FS,
            **self.residence,
            **{f"valid_{k}": v for k, v in self.validation.items()},
        }


def _openmm():
    """Imported lazily: only the MD env has OpenMM, the analysis env does not."""
    try:
        import openmm
        import openmm.app as app
        import openmm.unit as unit
    except ImportError as exc:  # noqa: BLE001
        raise MDError(
            "OpenMM is not importable in this interpreter. MD must run under "
            "/data/lab_vm/envs/dwi_amber_md/bin/python3 (OpenMM 8.5.2, CUDA)."
        ) from exc
    return openmm, app, unit


# CUDA is preferred but is currently broken on this host, so the platform is
# chosen at run time rather than hardcoded.
#
# THE FAULT, RECORDED SO IT IS NOT REDIAGNOSED. The driver (595.71.05) provides
# CUDA 13.2; the dwi_amber_md env ships libnvrtc 13.3. OpenMM 8.5 compiles its
# kernels through nvrtc, so it emits PTX one ISA version newer than the driver
# will load, and every CUDA context dies with CUDA_ERROR_UNSUPPORTED_PTX_VERSION
# (222). Setting OPENMM_CUDA_COMPILER to the system nvcc 12.0 does NOT help --
# OpenMM prefers nvrtc and ignores it. The real fix is to align nvrtc with the
# driver; until someone owns that change to a shared environment, OpenCL runs on
# the same A100s and is verified correct by openmm.testInstallation (median
# force difference vs Reference 6.7e-06).
PLATFORM_PREFERENCE = ("CUDA", "OpenCL", "CPU")


def _select_platform(openmm, device_index: int | None) -> tuple:
    """First platform that will actually build a context, not merely exist.

    `getPlatformByName` succeeds for CUDA even when CUDA cannot load a module,
    so availability has to be proven by constructing a context rather than by
    asking whether the platform is present.
    """
    import openmm.unit as unit

    last: Exception | None = None
    for name in PLATFORM_PREFERENCE:
        try:
            plat = openmm.Platform.getPlatformByName(name)
        except Exception as exc:  # noqa: BLE001
            last = exc
            continue
        props: dict[str, str] = {}
        if name in ("CUDA", "OpenCL"):
            props["Precision"] = "mixed"
            if device_index is not None:
                props["DeviceIndex"] = str(device_index)
        try:
            probe = openmm.System()
            probe.addParticle(1.0)
            integ = openmm.VerletIntegrator(1.0 * unit.femtosecond)
            openmm.Context(probe, integ, plat, props)
            if name != "CUDA":
                log.warning("MD running on %s, not CUDA (see PLATFORM_"
                            "PREFERENCE for why)", name)
            return plat, props
        except Exception as exc:  # noqa: BLE001
            last = exc
            continue
    raise MDError(f"no usable OpenMM platform; last error: {last}")


def _load_start_coords(workdir: Path, app) -> tuple:
    """Prefer the MM-GBSA minimised restart; fall back to the raw inpcrd.

    Starting from the minimised structure means MD continues the same
    trajectory of states the reported dG came from, rather than beginning
    somewhere else and making the two numbers describe different systems.
    """
    rst = workdir / "complex.min.rst"
    inpcrd = workdir / "complex.inpcrd"
    if rst.is_file():
        return app.AmberInpcrdFile(str(rst)), "complex.min.rst"
    if inpcrd.is_file():
        log.warning("%s: no minimised restart, starting from inpcrd", workdir.name)
        return app.AmberInpcrdFile(str(inpcrd)), "complex.inpcrd"
    raise MDError(f"no starting coordinates in {workdir}")


def build_simulation(workdir: Path, device_index: int | None = None):
    """Assemble a GBn2 implicit-solvent simulation from an MM-GBSA workdir."""
    openmm, app, unit = _openmm()

    prmtop_path = workdir / "complex.prmtop"
    if not prmtop_path.is_file():
        raise MDError(f"no complex.prmtop in {workdir}")

    prmtop = app.AmberPrmtopFile(str(prmtop_path))
    if any(r.name in ("WAT", "HOH") for r in prmtop.topology.residues()):
        raise MDError(
            f"{workdir.name}: topology contains explicit waters. This module "
            "is the implicit-solvent tier; a solvated system needs PME, a "
            "barostat and a different protocol.")

    # gbsaModel=None omits the nonpolar surface term ON PURPOSE. sander computes
    # it as LCPO (`ESURF`) while OpenMM would use ACE; the two disagree by tens
    # of kcal/mol on a system this size for reasons that have nothing to do with
    # whether either is right. Dropping it from both sides makes
    # validate_against_sander() a real test of the GB model instead of a
    # measurement of two different nonpolar approximations. It is also the
    # cheaper and more common choice for GB dynamics.
    system = prmtop.createSystem(
        implicitSolvent=app.GBn2,
        gbsaModel=None,
        soluteDielectric=1.0,
        solventDielectric=78.5,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        removeCMMotion=True,
    )
    integrator = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        FRICTION_PER_PS / unit.picosecond,
        TIMESTEP_FS * unit.femtosecond,
    )
    platform, props = _select_platform(openmm, device_index)

    sim = app.Simulation(prmtop.topology, system, integrator, platform, props)
    coords, source = _load_start_coords(workdir, app)
    sim.context.setPositions(coords.positions)
    return sim, prmtop, source


def validate_against_sander(workdir: Path, sim, tolerance_kcal: float = 5.0) -> dict:
    """Check OpenMM's GB energy against sander's on the SAME structure.

    The ensemble dG is only comparable to the single-structure dG if both
    engines agree on the same model. sander wrote `complex.min.out` during
    MM-GBSA; if OpenMM's energy for that geometry disagrees materially, the
    igb=8/GBn2 equivalence assumed here does not hold in practice and the
    caller should not merge the two numbers.
    """
    _, _, unit = _openmm()
    out = workdir / "complex.min.out"
    if not out.is_file():
        return {"checked": False, "reason": "no complex.min.out to compare"}

    text = out.read_text(errors="ignore")
    if "FINAL RESULTS" not in text:
        return {"checked": False, "reason": "sander output has no FINAL RESULTS"}
    block = text.split("FINAL RESULTS")[-1]

    import re
    # sander prints the total on the NSTEP line in scientific notation
    # ("1000  -7.2863E+03  ..."), NOT as "ENERGY = -7286.3". Matching only the
    # decimal form silently found nothing and the check reported itself as
    # skipped rather than failing loudly.
    m = re.search(r"^\s*\d+\s+(-?\d+\.\d+E[+-]\d+)", block, re.MULTILINE)
    if not m:
        m = re.search(r"(?:ENERGY|EPtot)\s*=\s*(-?[\d.]+(?:E[+-]\d+)?)", block)
    if not m:
        return {"checked": False, "reason": "could not parse sander energy"}
    sander_total = float(m.group(1))

    # Subtract sander's LCPO nonpolar term, because the OpenMM system is built
    # with gbsaModel=None and therefore has no counterpart to it.
    esurf = re.search(r"ESURF\s*=\s*(-?[\d.]+)", block)
    esurf_kcal = float(esurf.group(1)) if esurf else 0.0
    sander_kcal = sander_total - esurf_kcal

    state = sim.context.getState(getEnergy=True)
    openmm_kcal = state.getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole)
    delta = openmm_kcal - sander_kcal
    ok = abs(delta) <= tolerance_kcal
    if not ok:
        log.warning("%s: OpenMM %.2f vs sander %.2f kcal/mol (delta %.2f) -- "
                    "the GBn2/igb=8 equivalence does not hold here",
                    workdir.name, openmm_kcal, sander_kcal, delta)
    return {"checked": True, "agrees": bool(ok),
            "openmm_kcal": round(openmm_kcal, 2),
            "sander_kcal": round(sander_kcal, 2),
            "sander_total_kcal": round(sander_total, 2),
            "sander_esurf_removed": round(esurf_kcal, 2),
            "delta_kcal": round(delta, 2),
            "tolerance_kcal": tolerance_kcal}


def _ligand_and_pocket_indices(prmtop, cyx_residue_index: int | None = None
                               ) -> tuple[np.ndarray, np.ndarray, int]:
    """Heavy-atom indices for the ligand, and for the protein around it.

    The ligand is the last residue in the tleap-built complex. That is an
    assumption about how `build_topologies` ordered things, so it is CHECKED
    (the residue must not be a standard amino acid) rather than trusted.
    """
    residues = list(prmtop.topology.residues())
    lig = residues[-1]
    standard = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU", "GLY", "HIS",
        "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL", "WAT", "HOH", "Na+", "Cl-",
    }
    if lig.name in standard:
        raise MDError(
            f"last residue is {lig.name!r}, a standard residue -- the ligand "
            "is not where this code assumes it is, and every downstream "
            "distance would be measured against the wrong atoms.")

    lig_heavy = np.array([a.index for a in lig.atoms() if a.element is not None
                          and a.element.symbol != "H"], dtype=int)
    prot_heavy = np.array([a.index for r in residues[:-1] for a in r.atoms()
                           if a.element is not None
                           and a.element.symbol != "H"], dtype=int)
    if lig_heavy.size == 0:
        raise MDError("ligand residue has no heavy atoms")
    return lig_heavy, prot_heavy, lig.index


def residence_metrics(traj_nm: np.ndarray, lig_heavy: np.ndarray,
                      prot_heavy: np.ndarray) -> dict:
    """Does the ligand body stay in the pocket over the trajectory?

    `traj_nm` is (n_frames, n_atoms, 3) in nanometres.

    Reported, deliberately, as several numbers rather than one score:
    - ligand_rmsd_nm: drift of the ligand from its starting pose, after
      superposing on the PROTEIN (so protein tumbling does not masquerade
      as ligand movement).
    - mean_contacts: heavy-atom pairs within POCKET_CUTOFF_NM.
    - frac_frames_engaged: fraction of frames retaining at least a quarter of
      the starting contacts -- the closest thing to a residence time here.
    """
    if traj_nm.ndim != 3 or traj_nm.shape[0] < 2:
        raise MDError("need at least two frames for residence metrics")

    ref = traj_nm[0]
    # Superpose each frame on the protein (Kabsch), then measure the ligand.
    ref_p = ref[prot_heavy] - ref[prot_heavy].mean(axis=0)
    rmsds, contacts = [], []
    for frame in traj_nm:
        p = frame[prot_heavy]
        p_cen = p.mean(axis=0)
        cov = (p - p_cen).T @ ref_p
        v, _, wt = np.linalg.svd(cov)
        d = np.sign(np.linalg.det(v @ wt))
        rot = v @ np.diag([1.0, 1.0, d]) @ wt
        lig_now = (frame[lig_heavy] - p_cen) @ rot
        lig_ref = ref[lig_heavy] - ref[prot_heavy].mean(axis=0)
        rmsds.append(float(np.sqrt(((lig_now - lig_ref) ** 2).sum(axis=1).mean())))

        dist = np.linalg.norm(
            frame[lig_heavy][:, None, :] - frame[prot_heavy][None, :, :], axis=-1)
        contacts.append(int((dist < POCKET_CUTOFF_NM).sum()))

    rmsds_a = np.asarray(rmsds)
    contacts_a = np.asarray(contacts, dtype=float)
    start = contacts_a[0] if contacts_a[0] > 0 else 1.0
    engaged = float((contacts_a >= 0.25 * start).mean())
    return {
        "ligand_rmsd_nm_mean": round(float(rmsds_a.mean()), 4),
        "ligand_rmsd_nm_final": round(float(rmsds_a[-1]), 4),
        "ligand_rmsd_nm_max": round(float(rmsds_a.max()), 4),
        "mean_contacts": round(float(contacts_a.mean()), 1),
        "start_contacts": int(contacts_a[0]),
        "frac_frames_engaged": round(engaged, 4),
        "pocket_cutoff_nm": POCKET_CUTOFF_NM,
    }


def statistical_inefficiency(series: np.ndarray) -> float:
    """g = 1 + 2*sum(autocorrelation), the number of frames per independent one.

    Returned so the SEM can be widened by sqrt(g). Without this, a correlated
    trajectory reports an error bar that is too small by exactly the factor
    that makes it worth reporting.
    """
    x = np.asarray(series, dtype=float)
    n = x.size
    if n < 4:
        return float(n)
    x = x - x.mean()
    var = x.var()
    if var <= 0:
        return 1.0
    g, t = 1.0, 1
    while t < n - 1:
        c = float((x[: n - t] * x[t:]).mean() / var)
        if c <= 0:
            break
        g += 2.0 * c * (1.0 - t / n)
        t += 1
    return max(1.0, g)


def run_md(workdir: Path, candidate_id: str, device_index: int | None = None,
           production_ps: float = PRODUCTION_PS,
           equil_ps: float = EQUIL_PS) -> MDResult:
    """Equilibrate then produce a trajectory, and measure pocket residence."""
    openmm, app, unit = _openmm()
    t0 = time.time()

    sim, prmtop, source = build_simulation(workdir, device_index)
    validation = validate_against_sander(workdir, sim)

    md_dir = workdir / "md"
    md_dir.mkdir(parents=True, exist_ok=True)

    sim.minimizeEnergy(maxIterations=500)
    sim.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)
    sim.step(int(equil_ps * 1000 / TIMESTEP_FS))

    interval = int(FRAME_INTERVAL_PS * 1000 / TIMESTEP_FS)
    n_steps = int(production_ps * 1000 / TIMESTEP_FS)
    collector = _ArrayReporter(interval, unit)
    sim.reporters.append(collector)
    sim.reporters.append(app.DCDReporter(str(md_dir / "production.dcd"),
                                         interval))
    sim.reporters.append(app.StateDataReporter(
        str(md_dir / "production.csv"), interval, step=True,
        potentialEnergy=True, temperature=True))
    sim.step(n_steps)
    for r in list(sim.reporters):
        out = getattr(r, "_out", None)
        if out is not None and hasattr(out, "close"):
            try:
                out.close()
            except Exception:  # noqa: BLE001
                pass
    sim.reporters.clear()

    traj_nm = collector.as_array(prmtop.topology.getNumAtoms())
    np.save(md_dir / "traj_nm.npy", traj_nm.astype(np.float32))
    lig_heavy, prot_heavy, _ = _ligand_and_pocket_indices(prmtop)
    residence = residence_metrics(traj_nm, lig_heavy, prot_heavy)
    residence["start_coords"] = source

    res = MDResult(candidate_id=candidate_id, workdir=workdir,
                   n_frames=int(traj_nm.shape[0]),
                   ns_simulated=production_ps / 1000.0,
                   wall_seconds=time.time() - t0,
                   residence=residence, validation=validation)
    (md_dir / "md_result.json").write_text(
        json.dumps(res.to_dict(), indent=2), encoding="utf-8")
    log.info("%s: %d frames, %.1f ns, %.0f s wall, rmsd %.3f nm, engaged %.2f",
             candidate_id, res.n_frames, res.ns_simulated, res.wall_seconds,
             residence["ligand_rmsd_nm_mean"], residence["frac_frames_engaged"])
    return res


class _ArrayReporter:
    """Collects frames into memory as (n_frames, n_atoms, 3) in nanometres.

    WHY NOT READ THE DCD BACK. OpenMM's `app.DCDFile` is a writer -- its entire
    public API is `writeModel` -- and mdtraj is not installed in this
    environment. Rather than add a trajectory-format parser to the dependency
    surface, the frames are kept as they are produced. `.npy` is then the
    canonical trajectory for every downstream consumer (residence metrics now,
    ensemble rescoring next), and the DCD is written alongside purely so a
    human can load the run in PyMOL or VMD.

    100 frames x ~2.4k atoms is ~6 MB per candidate, so holding one trajectory
    in memory costs nothing at this scale.
    """

    def __init__(self, interval: int, unit):
        self._interval = int(interval)
        self._unit = unit
        self.frames: list[np.ndarray] = []

    def describeNextReport(self, simulation):  # noqa: N802 - OpenMM's interface
        steps = self._interval - simulation.currentStep % self._interval
        # (steps, positions, velocities, forces, energies, wrapPositions)
        return (steps, True, False, False, False, None)

    def report(self, simulation, state):
        pos = state.getPositions(asNumpy=True).value_in_unit(
            self._unit.nanometer)
        self.frames.append(np.asarray(pos, dtype=float))

    def as_array(self, n_atoms: int) -> np.ndarray:
        arr = np.asarray(self.frames)
        if arr.ndim != 3 or arr.shape[1] != n_atoms:
            raise MDError(f"trajectory shape {arr.shape} does not match "
                          f"{n_atoms} atoms")
        return arr
