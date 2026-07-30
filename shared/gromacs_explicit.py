"""
Purpose: Explicit-solvent MD in GROMACS for the non-covalent approaches.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-29
Input: an existing MM-GBSA workdir (ligand.mol2, ligand.frcmod, receptor_cys.pdb)
Output: a solvated, equilibrated production trajectory + per-stage logs

WHY THIS EXISTS. Every dynamics run in this project so far has been GB
IMPLICIT solvent, which has no water molecules at all -- it approximates their
averaged effect with a formula. That is the one fidelity axis D0036 could not
exclude, and it is most suspect exactly here: T_1 and T_2 produce flexible,
lipophilic, non-covalent binders, and the force that drives a greasy molecule
into a pocket is water being displaced from it. Implicit solvent cannot
represent a water bridging ligand and protein, or an ordered water that has to
be paid for on binding.

IT IS NOT A RANKING. Nothing here feeds the enrichment gate or reorders a
shortlist. D0036 established that better sampling did not rescue the ranking,
and this tier is a different question: does the complex behave the same way when
the water is real? Reported per candidate, not as a score.

SOLVATE ON THE AMBER SIDE, DELIBERATELY. `gmx solvate` inserts SPC water
coordinates, while an Amber topology expects TIP3P; the two differ in geometry
and charges. tleap adds TIP3P water and neutralising counter-ions with the
parameters that belong to them, and parmed then converts the whole solvated
system. This also reuses the existing GAFF2 ligand and ff19SB protein
parameterisation rather than re-deriving it, so the explicit-solvent run
describes the same molecule the rest of the project scored.

THE GPU NEEDED A DIFFERENT BUILD. The conda-forge GROMACS in `dwi_amber_md` is
compiled with OpenCL, and modern GROMACS refuses NVIDIA devices under OpenCL
(CUDA superseded it) -- `mdrun -gpu_id 0` fails with "incompatible devices",
which reads like a hardware problem and is not one. A CUDA build
(`dwi_gromacs_cuda`, GROMACS 2026.3) sees all 8 A100s and runs the same 30k-atom
system at 740 ns/day against 85 ns/day on 16 CPU cores. Recorded because the
OpenCL binary is still first on PATH and will mislead the next person.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)

AMBER_ENV = Path("/data/lab_vm/envs/dwi_amber_md")
GMX_ENV = Path("/data/lab_vm/envs/dwi_gromacs_cuda")   # CUDA, not the OpenCL one

WATER_MODEL = "TIP3PBOX"
BOX_PADDING_A = 12.0
TEMPERATURE_K = 300.0
PRESSURE_BAR = 1.0

# Stage lengths. Production defaults to 10 ns, which at ~740 ns/day is ~20 min
# per candidate on one A100 -- short for a binding study, ample for asking
# whether a docked pose survives contact with real water.
NVT_PS = 100.0
NPT_PS = 200.0
PRODUCTION_PS = 10000.0
FRAME_INTERVAL_PS = 10.0


class GromacsError(RuntimeError):
    """A stage failed, or produced output that cannot be trusted."""


@dataclass
class StageResult:
    name: str
    seconds: float
    ns_per_day: float | None = None


def _bin(env: Path, tool: str) -> str:
    p = env / "bin" / tool
    if not p.is_file():
        raise GromacsError(f"{tool} not found at {p}")
    return str(p)


def _run(cmd: list[str], cwd: Path, log_name: str, timeout: int = 86400) -> str:
    """Run one command, tee its output, and fail loudly with the tail."""
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                          timeout=timeout)
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    (cwd / log_name).write_text(out, encoding="utf-8")
    if proc.returncode != 0:
        raise GromacsError(f"{Path(cmd[0]).name} failed ({proc.returncode}); "
                           f"see {cwd/log_name}\n{out[-800:]}")
    return out


SOLVATE_LEAP = """source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams ligand.frcmod
LIG = loadmol2 ligand.mol2
rec = loadpdb receptor_cys.pdb
cpx = combine {{rec LIG}}
solvatebox cpx {water} {pad}
addions cpx Na+ 0
addions cpx Cl- 0
saveamberparm cpx solv.prmtop solv.inpcrd
quit
"""


def solvate(src: Path, wd: Path) -> dict:
    """Build an explicitly solvated complex with tleap, and verify it."""
    wd.mkdir(parents=True, exist_ok=True)
    for f in ("ligand.mol2", "ligand.frcmod", "receptor_cys.pdb"):
        s = src / f
        if not s.is_file():
            raise GromacsError(f"{src.name}: missing {f}")
        shutil.copy2(s, wd / f)

    (wd / "solvate.leap").write_text(
        SOLVATE_LEAP.format(water=WATER_MODEL, pad=BOX_PADDING_A),
        encoding="utf-8")
    _run([_bin(AMBER_ENV, "tleap"), "-f", "solvate.leap"], wd, "tleap.log",
         timeout=3600)

    prm, crd = wd / "solv.prmtop", wd / "solv.inpcrd"
    if not prm.is_file() or prm.stat().st_size == 0:
        raise GromacsError(f"{wd.name}: tleap produced no solvated topology")

    # A solvated system that somehow has no water, or a non-neutral charge, is
    # a silent failure mode: tleap reports Errors = 0 and writes a file either
    # way. Both are checked rather than assumed.
    import parmed
    p = parmed.load_file(str(prm), xyz=str(crd))
    waters = sum(1 for r in p.residues if r.name in ("WAT", "HOH"))
    ions = sum(1 for r in p.residues if r.name in ("Na+", "Cl-"))
    charge = sum(a.charge for a in p.atoms)
    if waters == 0:
        raise GromacsError(f"{wd.name}: solvated system contains no water")
    if p.box is None:
        raise GromacsError(f"{wd.name}: solvated system has no periodic box")
    if abs(charge) > 0.05:
        raise GromacsError(
            f"{wd.name}: net charge {charge:+.3f} after addions; the system is "
            "not neutral and PME will silently absorb the excess")

    p.save(str(wd / "sys.top"), overwrite=True)
    p.save(str(wd / "sys.gro"), overwrite=True)
    return {"atoms": len(p.atoms), "waters": waters, "ions": ions,
            "net_charge": round(float(charge), 4),
            "box_a": [round(float(x), 2) for x in p.box[:3]]}


_COMMON = """cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
constraints     = h-bonds
nstlist         = 40
"""

# REPLICATE SEEDS ARE EXPLICIT. GROMACS defaults gen-seed = -1, a
# pseudo-random value taken from the process, so two runs of one system diverge
# and NEITHER can be reproduced. That is exactly what happened to the implicit
# tier: the same molecule gave a mean ligand RMSD of 9.0 nm once and 1.75 nm the
# next time, and there was no way to go back and see which was typical
# (D0038). A replicate must be independent of its siblings and reproducible on
# its own, which requires the seed to be chosen rather than drawn.
BASE_SEED = 20260730


def replicate_seed(candidate_id: str, replicate: int) -> int:
    """A stable, distinct velocity seed for (candidate, replicate).

    Derived from the candidate id so two candidates never share a seed, and
    from the replicate index so siblings differ. Rerunning replicate 3 of a
    given candidate reproduces replicate 3, not a new trajectory.
    """
    import zlib
    h = zlib.crc32(candidate_id.encode()) & 0x7FFFFFFF
    return (BASE_SEED + h + replicate * 7919) % 2147483647


MDP = {
    "min": "integrator = steep\nnsteps = 5000\nemtol = 500.0\n" + _COMMON,
}


def nvt_mdp(seed: int) -> str:
    return (f"integrator = md\ndt = 0.002\nnsteps = {int(NVT_PS*500)}\n"
            f"tcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\n"
            f"ref-t = {TEMPERATURE_K}\ngen-vel = yes\n"
            f"gen-temp = {TEMPERATURE_K}\ngen-seed = {seed}\n"
            f"ld-seed = {seed}\n" + _COMMON)


def npt_mdp(seed: int) -> str:
    return (f"integrator = md\ndt = 0.002\nnsteps = {int(NPT_PS*500)}\n"
            f"tcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\n"
            f"ref-t = {TEMPERATURE_K}\n"
            f"pcoupl = C-rescale\npcoupltype = isotropic\ntau-p = 2.0\n"
            f"ref-p = {PRESSURE_BAR}\ncompressibility = 4.5e-5\n"
            f"continuation = yes\nld-seed = {seed}\n" + _COMMON)


def production_mdp(ps: float, seed: int = BASE_SEED) -> str:
    return (f"integrator = md\ndt = 0.002\nnsteps = {int(ps*500)}\n"
            f"ld-seed = {seed}\n"
            f"tcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\n"
            f"ref-t = {TEMPERATURE_K}\n"
            f"pcoupl = C-rescale\npcoupltype = isotropic\ntau-p = 2.0\n"
            f"ref-p = {PRESSURE_BAR}\ncompressibility = 4.5e-5\n"
            f"continuation = yes\n"
            f"nstxout-compressed = {int(FRAME_INTERVAL_PS*500)}\n"
            f"nstenergy = {int(FRAME_INTERVAL_PS*500)}\n" + _COMMON)


def _stage(wd: Path, name: str, mdp: str, start_gro: str, start_cpt: str | None,
           gpu_id: int | None, threads: int) -> StageResult:
    import time
    (wd / f"{name}.mdp").write_text(mdp, encoding="utf-8")
    grompp = [_bin(GMX_ENV, "gmx"), "grompp", "-f", f"{name}.mdp",
              "-c", start_gro, "-p", "sys.top", "-o", f"{name}.tpr",
              "-maxwarn", "5"]
    if start_cpt:
        grompp += ["-t", start_cpt]
    _run(grompp, wd, f"grompp_{name}.log", timeout=3600)

    t0 = time.time()
    mdrun = [_bin(GMX_ENV, "gmx"), "mdrun", "-deffnm", name,
             "-ntmpi", "1", "-ntomp", str(threads)]
    if gpu_id is not None:
        mdrun += ["-gpu_id", str(gpu_id)]
    try:
        _run(mdrun, wd, f"mdrun_{name}.log")
    except GromacsError:
        # GROMACS' GPU update kernel requires constrained atoms to be
        # contiguous, and an Amber-ordered topology need not satisfy that:
        # "Update groups can not be used for this system because atoms that are
        # (in)directly constrained together are interdispersed with other
        # atoms", followed by a segfault (signal 11) rather than a clean exit.
        # It is a kernel limitation, not a bad system -- the same topology runs
        # correctly with the update on the CPU. Retried once that way rather
        # than dropping the candidate, because a silent gap in the set would
        # look like a candidate that failed chemically.
        text = (wd / f"mdrun_{name}.log").read_text(errors="ignore")
        if "Update groups can not be used" not in text:
            raise
        log.warning("%s/%s: GPU update kernel rejected this atom ordering; "
                    "retrying with -update cpu", wd.name, name)
        _run(mdrun + ["-update", "cpu"], wd, f"mdrun_{name}_updatecpu.log")
    secs = time.time() - t0

    nsday = None
    logf = wd / f"{name}.log"
    if logf.is_file():
        for line in logf.read_text(errors="ignore").splitlines():
            if line.strip().startswith("Performance:"):
                try:
                    nsday = float(line.split()[1])
                except (IndexError, ValueError):
                    pass
    return StageResult(name=name, seconds=round(secs, 1), ns_per_day=nsday)


def run_pipeline(src: Path, wd: Path, gpu_id: int | None = None,
                 threads: int = 8, production_ps: float = PRODUCTION_PS,
                 replicate: int = 1, candidate_id: str | None = None) -> dict:
    """Solvate once, then run ONE replicate from that shared system.

    Layout: solvation and minimisation live in `wd`; each replicate gets
    `wd/repN` for its own NVT, NPT and production. Water placement and the
    minimised structure are therefore identical across replicates and only the
    VELOCITIES differ, which is what makes them replicates of one system rather
    than five loosely related simulations.

    Solvation is skipped when it already exists, so replicates 2..N cost only
    their own dynamics.
    """
    cid = candidate_id or wd.name
    seed = replicate_seed(cid, replicate)

    if (wd / "sys.top").is_file() and (wd / "sys.gro").is_file():
        import parmed
        pm = parmed.load_file(str(wd / "solv.prmtop"), xyz=str(wd / "solv.inpcrd"))
        info = {"atoms": len(pm.atoms),
                "waters": sum(1 for r in pm.residues
                              if r.name in ("WAT", "HOH")),
                "ions": sum(1 for r in pm.residues
                            if r.name in ("Na+", "Cl-")),
                "box_a": [round(float(x), 2) for x in pm.box[:3]],
                "solvation_reused": True}
    else:
        info = solvate(src, wd)
        info["solvation_reused"] = False
        log.info("%s: solvated %d atoms (%d waters, %d ions)", cid,
                 info["atoms"], info["waters"], info["ions"])

    stages = []
    if not (wd / "min.gro").is_file():
        stages.append(_stage(wd, "min", MDP["min"], "sys.gro", None, gpu_id,
                             threads))

    rep = wd / f"rep{replicate}"
    rep.mkdir(parents=True, exist_ok=True)
    # grompp resolves -p and -c relative to its own cwd, so the shared files are
    # linked into the replicate directory rather than copied: five replicates of
    # a 30k-atom system would otherwise duplicate ~30 MB of topology each.
    for f in ("sys.top", "min.gro"):
        link = rep / f
        if not link.exists():
            link.symlink_to(wd / f)

    stages.append(_stage(rep, "nvt", nvt_mdp(seed), "min.gro", None, gpu_id,
                         threads))
    stages.append(_stage(rep, "npt", npt_mdp(seed), "nvt.gro", "nvt.cpt",
                         gpu_id, threads))
    stages.append(_stage(rep, "prod", production_mdp(production_ps, seed),
                         "npt.gro", "npt.cpt", gpu_id, threads))

    traj = rep / "prod.xtc"
    if not traj.is_file() or traj.stat().st_size == 0:
        raise GromacsError(f"{cid} rep{replicate}: production produced no "
                           "trajectory")

    return {
        **info,
        "replicate": replicate,
        "velocity_seed": seed,
        "production_ps": production_ps,
        "frame_interval_ps": FRAME_INTERVAL_PS,
        "water_model": "TIP3P",
        "explicit_solvent": True,
        "gromacs": "2026.3 CUDA (dwi_gromacs_cuda)",
        "stages": {s.name: {"seconds": s.seconds, "ns_per_day": s.ns_per_day}
                   for s in stages},
        "trajectory": str(traj),
        "not_a_ranking": "explicit-solvent behaviour of one complex; feeds no "
                         "gate and reorders no shortlist (D0036)",
    }
