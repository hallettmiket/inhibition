"""
Purpose: the cheap attack-geometry gate — free pose check, then a 10 ns sweep.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: candidates + their elevated poses
Output: 00_outputs/blacksmith/attack_sweep/attack_sweep_<N>.csv, ranked

Implements #32 and `docs/prereg_attack_sweep.md`. @tt8804: *"100% residence can
still result in a warhead that does not stay near the attack distance or angle...
we should do 10 ns sweeps to first filter for that ability to start at cys113 and
stay near before comitting to 100 ns MD."*

WHY THE FUNNEL NEEDED THIS. Residence and attack geometry are nearly independent
across the six molecules measured so far -- rho = +0.319, p = 0.538. The two with
PERFECT residence are attack-ready 0.4% and 1.7% of frames; the one that
dissociated at 81 ns is attack-ready 55.2%. Ranking on residence selects against
the thing we want.

TWO STAGES, AND THE FIRST IS FREE.

  STAGE 0  the elevated pose's OWN geometry, before any simulation. Only 1 of 6
           poses elevated in the bornite cohort started inside the attack window
           at all, and it is the only one that got anywhere. If that separation
           holds it is a filter that costs nothing, and no sweep should be run
           on a pose that starts 9 A away pointing the wrong direction.

  STAGE 1  10 ns of unbiased MD, then attack geometry per frame. ~0.4 GPU-h
           against ~4 for a full run, so it pays for itself on one rejection.

BOTH OBSERVABLES ARE REPORTED, BECAUSE THEY MAY DIVERGE.
  frac_attack_ready  how LONG the molecule is competent
  n_visits           how MANY independent excursions it makes into attack
                     geometry -- and a covalent reaction needs ONE good approach,
                     not sustained occupancy, so this is the more mechanistically
                     honest quantity. On the six existing trajectories the two
                     rank identically (rho = +1.000); that may not survive a
                     cohort with more range, which is why both are kept.

THIS RANKS, IT DOES NOT THRESHOLD. There is no evidence for a cut-off value --
one molecule above 5% and five below is not a threshold, it is one molecule. The
output is an ordering; how deep to elevate is a separate decision with its own
cost argument.
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import md_movie as mov                  # noqa: E402
from shared import nac_criterion as nac             # noqa: E402
from shared import outputs as sout                  # noqa: E402
from shared import target_config as tc              # noqa: E402
from shared import md_movie as mdm
from shared import run_paths as rp        # noqa: E402

log = logging.getLogger("attack-sweep")
MD = rp.residence_work()
#: The sweep writes to its OWN root. The workdir is <root>/<candidate>/md/rep1
#: regardless of tag, so a 10 ns sweep and a later 100 ns run of the same
#: molecule would otherwise collide -- and the 100 ns run would find a finished
#: 10 ns prod.xtc sitting there and skip itself.
SWEEP_ROOT = rp.sweep_work()
POSES = rp.BLACKSMITH
OUT = sout.Topic("blacksmith", rp.sweep_topic())
PY = Path.home() / ".micromamba/envs/dwi_reactive/bin/python"

#: Triage-sweep length, from config. D0085 set it to 8 ns -- the plateau, with
#: 7 ns the hard floor and nothing bought past 10 -- but this stayed a 10 ns
#: literal, so the constant and the spec disagreed by 2 ns on every sweep: 25%
#: more GPU time per mode than the experiment concluded was needed.
SWEEP_PS = tc.md_sweep_ps()
FRAMES = 500               # 20 ps resolution over the sweep

#: A visit has to LAST to count as a visit (#34, adversary audit).
#:
#: Measured on the six 100 ns trajectories, the median attack-ready episode is
#: one or two frames for five of the six -- 100% of two molecules' episodes are
#: <=200 ps. Those are boundary recrossings, not resolved approaches. With a mean
#: episode of ~1 frame the visit count is algebraically almost
#: `frac_attack_ready x n_frames`, which is why the two ranked identically
#: (rho = +1.000) and why calling them separate observables was not defensible.
#:
#: IN PICOSECONDS, NOT FRAMES, because a frame rule makes the count a function of
#: the save interval: one molecule gave 57/26/14/7/3 visits at
#: 100 ps/200 ps/500 ps/1 ns/2 ns. The sweep saves every 20 ps and the validation
#: ran at 100 ps, so a frame-based count would not be comparable between them.
MIN_DWELL_PS = 100.0


def _mp():
    spec = importlib.util.spec_from_file_location(
        "mdprio_report", REPO / "scripts" / "mdprio_report.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["mdprio_report"] = m
    spec.loader.exec_module(m)
    return m


def competent(angle: np.ndarray, kind: str) -> np.ndarray:
    """Angular competence by the mechanism's own criterion, never one constant."""
    if "anti" in (kind or ""):
        return angle >= nac.SN2_ANGLE_MIN
    return angle <= nac.PERPENDICULAR_MAX_OFF_NORMAL


def _episodes(ready: np.ndarray) -> list[tuple[int, int]]:
    """(start, length) of every contiguous run of True."""
    out, i, n = [], 0, len(ready)
    while i < n:
        if ready[i]:
            j = i
            while j < n and ready[j]:
                j += 1
            out.append((i, j - i))
            i = j
        else:
            i += 1
    return out


def attack_ready_max_a() -> float:
    """Upper distance bound for "attack ready", from config.

    NOT `nac.NAC_DIST_MAX`. That is the near-attack WINDOW (4.2 A) and using it
    here meant a mode whose warhead sat at a trajectory median of 3.6 A scored
    93% "attack ready" -- true by that definition, and not what the number is
    read as. The worklist selects modes at < 3.0 A; the readout now agrees with
    the selection (@twu383, 2026-09-02).
    """
    try:
        from shared import target_config as _tc
        v = _tc.get("md.attack_ready_max_a", default=None)
        if v:
            return float(v)
    except Exception:                                      # noqa: BLE001
        pass
    return float(nac.NAC_DIST_MAX)


def _cfg(key: str, default):
    """One config read, failing to the value that reproduces the old behaviour."""
    try:
        from shared import target_config as _tc
        v = _tc.get(f"md.{key}", default=None)
        if v is not None:
            return v
    except Exception:                                      # noqa: BLE001
        pass
    return default


def attack_ready_use_angle() -> bool:
    """Whether the angular criterion is part of "attack ready".

    FALSE since 2026-09-02 (@twu383). Not because the angle is meaningless --
    it is measured and reported either way -- but because it is not what limits
    the answer: at a 3.0 A cutoff, adding it moved the count of discriminating
    modes from 1 to 0. And on a set already selected on distance it is the term
    D0110 showed is not class-neutral, so gating on it re-weights the campaign
    towards BDHI for a steric reason rather than a chemical one.

    Defaults TRUE so a caller with no config reproduces the previous numbers.
    """
    return bool(_cfg("attack_ready_use_angle", True))


def elevation_thresholds() -> dict:
    """The 100 ns gate, from config. One place, so the GUI cannot drift from it."""
    return {
        "occupancy_min": float(_cfg("elevate_occupancy_min", 0.60)),
        "rmsd_max_a": float(_cfg("elevate_rmsd_max_a", 3.5)),
        "rmsd_mean_a": float(_cfg("elevate_rmsd_mean_a", 3.0)),
    }


def rmsd_stats(rep: Path) -> dict:
    """Ligand RMSD over the production run, in ANGSTROM, from `rmsd.xvg`.

    ALREADY ON DISK for every completed sweep -- `attack_sweep` runs `gmx rms`
    beside the trajectory -- so this is arithmetic over existing files, not a
    reason to re-simulate anything.

    UNITS ARE THE TRAP. GROMACS writes nm and every threshold the user states is
    in Angstrom; a missing factor of 10 makes every run pass a 3.5 bar by a mile
    and reads as a spectacular result. The conversion happens HERE, once, and
    the columns are named `_a` so a nm value cannot be mistaken for one.

    Returns {} when the file is absent or unparseable -- an ABSENT reading, not
    a passing one. `elevation_verdict` refuses to elevate on a missing RMSD
    rather than treating it as zero.
    """
    f = Path(rep) / "rmsd.xvg"
    if not f.is_file():
        return {}
    ys = []
    try:
        for ln in f.read_text(errors="replace").splitlines():
            ln = ln.strip()
            if not ln or ln[0] in "#@":
                continue
            parts = ln.split()
            if len(parts) >= 2:
                ys.append(float(parts[1]))
    except (OSError, ValueError):
        return {}
    if not ys:
        return {}
    a = np.asarray(ys, dtype=float) * 10.0                 # nm -> Angstrom
    return {
        "rmsd_max_a": float(a.max()),
        "rmsd_mean_a": float(a.mean()),
        "rmsd_median_a": float(np.median(a)),
        "rmsd_final_a": float(a[-1]),
        "rmsd_frames": int(a.size),
    }


def elevation_verdict(rec: dict) -> dict:
    """Does this mode earn a 100 ns run? (@twu383, 2026-09-02.)

    Two independent conditions, BOTH required:

      1. the pose held      -- rmsd_max_a < 3.5 OR rmsd_mean_a < 3.0
      2. the warhead stayed -- frac_attack_ready >= 0.60 within 3.5 A

    The OR in (1) is the spike allowance the request asked for: a run that
    touches 4 A briefly and sits low otherwise passes on the mean; one that
    never exceeds 3.5 passes without it.

    They are separate because they can disagree. Occupancy is about the WARHEAD
    reaching Cys113; RMSD is about the whole LIGAND staying where it was docked.
    A molecule can pivot its warhead into place while its scaffold walks off,
    and a molecule can sit perfectly still facing the wrong way. Either alone
    would elevate one of those.

    A MISSING READING NEVER ELEVATES. If `rmsd_max_a` is absent the verdict is
    "unknown", not "pass" -- the failure mode this project keeps producing is a
    guard that passes because the thing it inspects is not there.
    """
    th = elevation_thresholds()
    occ = rec.get("frac_attack_ready")
    mx, mn = rec.get("rmsd_max_a"), rec.get("rmsd_mean_a")

    out = {"elevate_occupancy_min": th["occupancy_min"],
           "elevate_rmsd_max_a": th["rmsd_max_a"],
           "elevate_rmsd_mean_a": th["rmsd_mean_a"]}
    if occ is None or mx is None or mn is None:
        out.update({"elevate": False, "elevate_why": "no reading",
                    "pose_held": None, "warhead_engaged": None})
        return out

    held = (float(mx) < th["rmsd_max_a"]) or (float(mn) < th["rmsd_mean_a"])
    engaged = float(occ) >= th["occupancy_min"]
    if held and engaged:
        why = "pass"
    elif not held and not engaged:
        why = "pose left and warhead not engaged"
    elif not held:
        why = f"pose left (max {float(mx):.2f} A, mean {float(mn):.2f} A)"
    else:
        why = f"warhead engaged {float(occ)*100:.0f}% < {th['occupancy_min']*100:.0f}%"
    out.update({"elevate": bool(held and engaged), "elevate_why": why,
                "pose_held": bool(held), "warhead_engaged": bool(engaged)})
    return out


def geometry_stats(dist: np.ndarray, angle: np.ndarray, kind: str,
                   frame_ps: float) -> dict:
    """Geometry readings for one trajectory.

    `frame_ps` is REQUIRED rather than defaulted: `n_visits` is meaningless
    without it (see MIN_DWELL_PS), and a default would silently produce a number
    that is not comparable to the one it is being validated against.
    """
    if not len(dist):
        # nac_criterion raises rather than returning a false verdict when it
        # cannot measure; this did the opposite and raised IndexError on
        # `ready[0]` from deep inside a worker, where it reads as a crash rather
        # than as "there is no trajectory here".
        raise ValueError("empty trajectory: no frames to score")
    # THE FLOOR IS THE PHYSICAL ONE, THE CEILING IS THE CAMPAIGN'S.
    # Below NAC_DIST_MIN the two atoms overlap (PoseBusters' C...S clash
    # threshold is 2.625 A), so a closer frame is a clash and not a better
    # approach. The ceiling comes from config, defaulting to the old window so
    # a run that does not set it reproduces the previous numbers exactly.
    hi = attack_ready_max_a()
    inw = (dist >= nac.NAC_DIST_MIN) & (dist < hi)
    # DISTANCE ONLY by default since 2026-09-02 (see `attack_ready_use_angle`).
    # The angular term is still computed and reported as
    # `frac_attack_ready_angle`, so dropping it from the gate does not drop it
    # from the record and the two can be compared on any finished run.
    comp = competent(angle, kind)
    use_ang = attack_ready_use_angle()
    ready = (inw & comp) if use_ang else inw
    # An excursion is a rising edge: not-ready -> ready. The first frame counts
    # as a visit if it is already ready, otherwise a pose that starts competent
    # and never leaves would be recorded as zero visits.
    raw_visits = int(np.sum(np.diff(ready.astype(int)) == 1) + (1 if ready[0] else 0))
    min_frames = max(1, int(round(MIN_DWELL_PS / frame_ps)))
    eps = _episodes(ready)
    visits = sum(1 for _, ln in eps if ln >= min_frames)
    return {
        "frac_in_window": float(inw.mean()),
        "frac_attack_ready": float(ready.mean()),
        "n_visits": visits,
        # Both are kept so the debounce is inspectable rather than assumed. A
        # large gap between them means the molecule is skimming the boundary,
        # not approaching.
        "n_visits_raw": raw_visits,
        "min_dwell_ps": MIN_DWELL_PS,
        "frame_ps": float(frame_ps),
        "median_episode_ps": float(np.median([ln for _, ln in eps]) * frame_ps)
                             if eps else 0.0,
        # STAMPED ON EVERY ROW, so a `frac_attack_ready` computed under one
        # definition can never be silently compared with one computed under
        # another. This is the column that says what the number means.
        "attack_ready_max_a": float(hi),
        "attack_ready_min_a": float(nac.NAC_DIST_MIN),
        "attack_ready_angle_deg": float(nac.PERPENDICULAR_MAX_OFF_NORMAL),
        # WHICH DEFINITION PRODUCED THE NUMBER, on the row itself. Without this
        # a distance-only fraction and a distance+angle one are two plausible
        # floats in the same column -- the exact shape this project keeps being
        # bitten by.
        "attack_ready_uses_angle": bool(use_ang),
        # The other definition, always measured. Free, and it is what makes the
        # choice reviewable later without re-reading 4,000 trajectories.
        "frac_attack_ready_angle": float((inw & comp).mean()),
        "start_dist_a": float(dist[0]),
        "start_angle_deg": float(angle[0]),
        "start_attack_ready": bool(ready[0]),
        "median_dist_a": float(np.median(dist)),
        "median_angle_deg": float(np.median(angle)),
        "min_dist_a": float(dist.min()),
    }


def _static_sg() -> np.ndarray:
    """Cys113 SG in the prepared receptor -- the SAME frame the MD uses.

    Verified 2026-09-02: `receptor_cys.pdb` in a sweep workdir and the docking
    receptor `3IKD_noligand.pdb` both put Cys113 SG at (13.385, 3.989, -2.040),
    residue 113 renumbering to 63 in the MD system. There is no frame
    transformation between docking and MD, which is why the drift seen at the
    start of production is PHYSICS and not a coordinate mismatch.
    """
    rec = rp.receptor_prep()
    for ln in Path(rec).read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        # BY RESIDUE NUMBER. 3IKD also has a Cys57, and taking the first CYS SG
        # would measure the approach to the wrong sulfur entirely.
        if (ln[17:20].strip() == "CYS" and ln[12:16].strip() == "SG"
                and ln[22:26].strip() == "113"):
            return np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    raise SweepError(f"no Cys113 SG in {rec}")


def _cand_meta(ident: str) -> dict:
    """The candidate's mechanism and reactive SMARTS, from the warhead library."""
    import pandas as _pd
    wh = _pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    frames = sorted((rp.DATA / "04_t4_combinatorial").glob("D4_*.parquet"),
                    key=lambda f: int(f.stem.split("_")[1]))
    df = _pd.read_parquet(frames[-1])
    row = df[df.candidate_id == ident]
    if row.empty:
        raise SweepError(f"{ident} is in no candidate frame")
    cls = str(row.iloc[0].warhead_class)
    w = wh[wh.class_id == cls]
    if w.empty:
        raise SweepError(f"warhead class {cls!r} not in the library")
    return {"mechanism": str(w.iloc[0].mechanism),
            "smarts": str(w.iloc[0].reactive_atom_smarts)}


def _pose_at_rank(pose_sdf: Path, pose_rank: int):
    """The pose carrying `pose_rank`, selected by its property (never position)."""
    from rdkit import Chem as _Chem
    for m in _Chem.SDMolSupplier(str(pose_sdf), removeHs=False):
        if m is None or not m.HasProp("pose_rank"):
            continue
        if int(m.GetProp("pose_rank")) == pose_rank:
            return m
    raise SweepError(f"no pose_rank {pose_rank} in {pose_sdf.name}")


def pose_mode(pose_sdf: Path, pose_rank: int) -> int | None:
    """The `mode` property of the pose with this `pose_rank`, or None.

    Read by identity: the pose whose own `pose_rank` matches, and then that
    pose's own `mode`. Nothing here counts positions in the file.
    """
    if not pose_sdf.is_file():
        return None
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False, sanitize=False):
            if m is None or not m.HasProp("pose_rank"):
                continue
            if int(m.GetProp("pose_rank")) == pose_rank:
                return int(m.GetProp("mode")) if m.HasProp("mode") else None
    except Exception as exc:                                # noqa: BLE001
        log.warning("%s: could not read mode (%s)", pose_sdf.name, exc)
    return None


class SweepError(RuntimeError):
    """A sweep could not be run or reused. Named so it cannot pass as a result."""


def _finished(rep: Path) -> bool:
    """True only if mdrun reported it finished this trajectory.

    Existence of `prod.xtc` is not completion -- see run_sweep. GROMACS writes
    "Finished mdrun" to the log on a clean exit, and that is the only marker
    here that distinguishes a finished run from one in progress.
    """
    log_f, xtc = rep / "prod.log", rep / "prod.xtc"
    if not (log_f.is_file() and xtc.is_file()):
        return False
    try:
        return "Finished mdrun" in log_f.read_text(errors="replace")
    except OSError:
        return False


class SweepAborted(SweepError):
    """The pose left the site during equilibration; production was not run.

    A distinct type because this is a RESULT, not a failure -- the molecule was
    measured and found to have gone. Recording it as `failed:` would make a
    real observation indistinguishable from a crashed GPU job.
    """

    def __init__(self, dist_a: float) -> None:
        super().__init__(f"warhead {dist_a:.2f} A from SG after equilibration")
        self.dist_a = dist_a


def _elev():
    """`elevation_run.read_gro` -- one .gro parser, not two."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "elev_for_sweep", REPO / "scripts" / "elevation_run.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def _gro_resids(gro: Path, n_atoms: int) -> list[int] | None:
    """Residue number per atom, in the same order `read_gro` returns atoms.

    A SEPARATE READ, so the two must be checked to agree. `read_gro` gives
    names and coordinates but not residue numbers; taking them from a second
    pass is only safe if both passes see the same atoms in the same order, so
    the count is asserted rather than assumed. A silent off-by-one here would
    shift every residue number by one and pick the wrong cysteine.
    """
    lines = gro.read_text(errors="replace").splitlines()
    if len(lines) < 3:
        return None
    body = lines[2:2 + n_atoms]
    if len(body) != n_atoms:
        log.warning("%s: %d atom lines for %d atoms; not measuring",
                    gro.name, len(body), n_atoms)
        return None
    try:
        return [int(l[0:5]) for l in body]
    except ValueError:
        log.warning("%s: unparseable residue numbers; not measuring", gro.name)
        return None


def _equil_distance(cand: str, rep: Path, pose: Path,
                    pose_rank: int) -> float | None:
    """Warhead-SG distance in the equilibrated frame, or None if not certain.

    RETURNS None RATHER THAN A GUESS, and the caller then runs production. An
    abort must be taken only on a distance we are sure of: skipping work on a
    bad number silently discards a molecule, where doing the work costs six
    minutes. Fail safe means DO the sweep.

    PERIODIC IMAGES ARE THE WHOLE REASON THIS FUNCTION IS CAREFUL. The first
    version read the .gro naively and reported 51.18 A for a ligand in a ~7 nm
    box -- the box length minus the real distance, which is large, finite and
    completely wrong. A pose sitting perfectly in the pocket but wrapped across
    a boundary would have been aborted. Minimum image is applied here, the box
    is required to be orthorhombic, and a result beyond half the box diagonal is
    refused outright as still-suspect.
    """
    gro = rep / "npt.gro"
    if not gro.is_file():
        log.warning("%s: no npt.gro to measure; running production", cand)
        return None
    try:
        E = _elev()
        names, xyz, box = E.read_gro(gro)
        if box.size >= 6 and float(np.abs(box[3:]).max()) > 1e-6:
            log.warning("%s: triclinic box; not measuring", cand)
            return None

        c = _cand_meta(cand)
        m0 = _pose_at_rank(pose, pose_rank)
        mt = m0.GetSubstructMatches(Chem.MolFromSmarts(c["smarts"]))[0]

        # The ligand's atoms, by RESIDUE NAME. The atom count must match the
        # pose or the index `mt[0]` addresses a different atom -- and a distance
        # from the wrong atom is a plausible wrong number, not an error.
        lig = [k for k, (resn, _) in enumerate(names)
               if resn.strip() in {"LIG", "UNL", "MOL"}]
        if len(lig) != m0.GetNumAtoms():
            log.warning("%s: gro has %d ligand atoms, pose has %d; not "
                        "measuring", cand, len(lig), m0.GetNumAtoms())
            return None
        # CYS113 BY RESIDUE NUMBER, THROUGH THE SHARED CONSTANT.
        #
        # 3IKD has Cys57 as well as Cys113 and BOTH are reduced CYS in the MD
        # system, so matching on resname alone finds two and cannot choose --
        # the first version correctly refused rather than guess, which is why
        # it never measured anything. The MD system renumbers from 1
        # (`md_movie.PIN1_OFFSET`), putting Cys113 at residue 63.
        #
        # The offset is not trusted: the residue found there must BE a cysteine
        # or this returns None. An offset that slipped by one would name a
        # different residue with a distance that still looks plausible, which is
        # the failure `elevation_report`'s own numbering guard exists to stop.
        resids = _gro_resids(gro, len(names))
        if resids is None:
            return None
        target = mdm.CYS113_RESI
        sg = [k for k, (resn, at) in enumerate(names)
              if resids[k] == target and at.strip() == "SG"
              and resn.strip() == "CYS"]
        if len(sg) != 1:
            log.warning("%s: %d CYS SG at residue %d (offset %d); not "
                        "measuring", cand, len(sg), target, mdm.PIN1_OFFSET)
            return None

        v = xyz[lig[mt[0]]] - xyz[sg[0]]
        mic = float(np.linalg.norm(v - box[:3] * np.round(v / box[:3]))) * 10.0
        half_diag = float(np.linalg.norm(box[:3])) * 10.0 / 2.0
        if mic > half_diag:
            log.warning("%s: %.1f A exceeds half the box diagonal (%.1f); "
                        "measurement suspect, running production",
                        cand, mic, half_diag)
            return None
        return mic
    except Exception as exc:                              # noqa: BLE001
        log.warning("%s: could not measure equilibrated frame (%s); "
                    "running production", cand, exc)
        return None


def run_sweep(cand: str, pose: Path, pose_rank: int, gpu: int,
              ps: float, net_charge: int | None,
              abort_above: float | None = None) -> Path | None:
    """10 ns of MD through the production script, so the physics is identical."""
    # ONE WORKDIR PER (MOLECULE, MODE). md_residence names the directory after
    # the candidate, so sweeping two modes of one molecule would put them in the
    # same place and the second would find the first's finished trajectory and
    # skip itself -- reporting mode 0's result as mode 4's.
    # THE LENGTH IS PART OF THE PATH, not just the tag.
    #
    # The resume guard checked only that `prod.xtc` existed, so a finished 200 ps
    # trajectory silently satisfied a 10 ns request -- caught tonight when a
    # smoke test's 200 ps run was reused as if it were the real sweep. That is
    # the same defect class the per-mode workdir fixed for pose_rank, with the
    # length dimension left open. A 10 ns answer read off 200 ps of dynamics
    # would be indistinguishable from a real one in every artefact downstream.
    #
    # AND THE GUARD MUST CHECK COMPLETION, NOT EXISTENCE. `prod.xtc` exists the
    # moment mdrun starts writing it, so an IN-PROGRESS run satisfied this test
    # and a second process read the partial trajectory as a finished one. Caught
    # 2026-08-11 the first time two workers were run: `t4_e0b03662d460_m1` was
    # 3.4 ns into its 10 ns when a second invocation analysed the same directory
    # and wrote `status: ok, frac_attack_ready 0.0` in four seconds. A partial
    # answer stamped ok is indistinguishable from a real one downstream -- the
    # same defect this docstring already records for the 200 ps case, with the
    # completeness dimension left open instead of the length one.
    root = SWEEP_ROOT / f"rank{pose_rank}_{int(ps)}ps"
    rep = root / cand / "md" / "rep1"
    if _finished(rep):
        log.info("%s: %d ps trajectory already complete, not re-running",
                 cand, int(ps))
        return rep
    if (rep / "prod.xtc").is_file():
        # Present but unfinished: either another process is running it right now,
        # or one died partway. Neither is a result. Refuse rather than analyse it.
        raise SweepError(
            f"{cand}: an unfinished {int(ps)} ps trajectory is already in "
            f"{rep} — another worker may be running it. Not analysing a partial "
            f"run as a complete one.")
    base = [str(PY), str(REPO / "scripts/md_residence_3ikd.py"),
            "--candidate", cand, "--pose", str(pose),
            "--pose-rank", str(pose_rank), "--production-ps", str(int(ps)),
            "--gpu", str(gpu), "--keep", "--tag", f"sweep_r{pose_rank}",
            "--work-root", str(root)]
    if net_charge is not None:
        base += ["--net-charge", str(net_charge)]

    # ---- EARLY GIVE-UP, after equilibration and before production ----------
    #
    # WHY HERE AND NOT EARLIER. The DOCKED geometry does not predict the sweep:
    # over the first nac_v8 modes, five poses all sitting at 2.78-2.97 A gave
    # outcomes from 0.000 to 0.926 attack-ready. What separates them is the
    # distance AFTER 300 ps of unrestrained equilibration -- 3.58 A gave 0.926,
    # 6.46 A gave 0.000 -- and that cannot be known without running it.
    #
    # Equilibration is 300 ps of the 1,500 a sweep runs, so a pose that has
    # already left the site is dropped for a fifth of the cost.
    #
    # THE THRESHOLD IS DELIBERATELY LOOSE AND THE REASON IS n. It is set from
    # FIVE observations, which is exactly the basis this project has a record
    # about not trusting (D0094, D0085). 6.0 A is well past the 4.2 A window and
    # past every pose that has scored above zero so far, so it catches only the
    # unambiguously departed. `equil_dist_a` is recorded on EVERY row whether it
    # aborts or not, so the cut can be re-derived from hundreds of runs later
    # rather than from these five.
    if abort_above:
        r0 = subprocess.run(base + ["--stop-after", "npt"],
                            capture_output=True, text=True, timeout=7200)
        if r0.returncode != 0:
            raise SweepError(f"{cand}: equilibration rc={r0.returncode} "
                             f"{(r0.stderr or '').strip()[-200:]}")
        d_eq = _equil_distance(cand, rep, pose, pose_rank)
        if d_eq is not None and d_eq > abort_above:
            log.info("%s: warhead %.2f A from SG after equilibration "
                     "(> %.1f) — skipping production", cand, d_eq, abort_above)
            raise SweepAborted(d_eq)
        log.info("%s: %.2f A after equilibration — running production",
                 cand, d_eq if d_eq is not None else float("nan"))
        cmd = base + ["--reuse-equilibration"]
    else:
        cmd = base
    log.info("%s: %.0f ps sweep on GPU %d", cand, ps, gpu)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
    if r.returncode != 0:
        log.warning("%s: sweep failed rc=%d %s", cand, r.returncode, r.stderr[-300:])
        raise SweepError(f"{cand}: rc={r.returncode} "
                         f"{(r.stderr or r.stdout or '').strip()[-200:]}")
    if not (rep / "prod.xtc").is_file():
        # EXIT 0 AND NO TRAJECTORY IS STILL A FAILURE, AND IT HAS A REASON. This
        # returned None, which the caller recorded as a bare "sweep failed" -- so
        # a molecule that could not be parameterised looked exactly like one whose
        # GPU run crashed, and neither said why. The reason is in the child's own
        # output; carry it instead of discarding it.
        tail = (r.stdout or "").strip().splitlines()
        raise SweepError(f"{cand}: exited 0 but wrote no prod.xtc"
                         + (f" — {tail[-1][:160]}" if tail else ""))
    return rep


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", nargs="+", required=True)
    # FROM `run.topic`, NOT A LITERAL. This defaulted to `nac_v3_poses`, so a
    # sweep launched after the 3.0.0 screen would have simulated 2.2.0's poses
    # while reporting them against 3.0.0's ranking -- 10 ns of GPU per mode spent
    # on the wrong structure, with nothing in the output to say which run it came
    # from. Same defect as D0080, one stage further downstream.
    ap.add_argument("--pose-dir", default=None,
                    help="representative poses; defaults to <run.topic>_poses")
    ap.add_argument("--pose-rank", type=int, default=1)
    ap.add_argument("--gpu", type=int, default=1)
    ap.add_argument("--sweep-ps", type=float, default=SWEEP_PS)
    ap.add_argument("--stage0-only", action="store_true",
                    help="report starting geometry only — costs no GPU at all")
    # A STAGE-0 ROW IS NOT A RESULT, AND MUST BE ABLE TO LAND SOMEWHERE ELSE.
    #
    # `--stage0-only` writes a row with `status = "stage0 only"` and no
    # measurements, into the SAME table as finished sweeps. Anything that
    # resumes by asking "is this ident present?" then treats a free geometry
    # probe as a completed simulation -- which is exactly what happened on
    # 2026-09-02: 12 modes were marked done by a probe and dropped out of the
    # worklist. `sweep_supervisor.done_tasks` now requires `status == "ok"`,
    # but the deeper fix is that a probe should not have to share the results
    # directory at all.
    # DEFAULT OFF until the equilibrated-frame measurement has been verified on
    # a real trajectory. It was briefly 6.0 with a naive .gro read that ignored
    # periodic images, and it wrongly aborted t4_b49ffa60a11a_m113 at a reported
    # 56.6 A -- the box length minus the real distance. An abort DISCARDS work,
    # so it stays opt-in: the caller asks for it, having read this.
    ap.add_argument("--abort-above-a", type=float, default=0.0, metavar="A",
                    help="after equilibration, if the warhead is further than "
                         "this from Cys113 SG, skip production and record the "
                         "distance. 0 disables. Default 6.0 -- LOOSE on "
                         "purpose: it is set from five observations and only "
                         "catches poses that have unambiguously left")
    ap.add_argument("--out-topic", default=None, metavar="NAME",
                    help="write rows to attack_sweep_<NAME> instead of the "
                         "current run's table; use for stage-0 probes so they "
                         "cannot be mistaken for finished sweeps")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.out_topic:
        global OUT
        OUT = sout.Topic("blacksmith", f"attack_sweep_{args.out_topic}")
        log.info("rows -> attack_sweep_%s (NOT the run's results table)",
                 args.out_topic)
    if args.pose_dir is None:
        from shared import target_config as tc
        # THROUGH THE RESOLVER, so a corrected pose set is picked up rather
        # than the superseded one this literal named. `rp.poses_dir` takes the
        # highest integer version present; the base directory has no suffix and
        # still resolves when no correction exists.
        args.pose_dir = str(rp.poses_dir(tc.get('run.topic')))
    log.info("poses from %s", Path(args.pose_dir).name)
    mp = _mp()

    rows = []
    for cand in args.candidates:
        pose = Path(args.pose_dir) / f"{cand}.sdf"
        # THE MODE IS READ FROM THE POSE, NOT DERIVED FROM ITS RANK.
        #
        # This was `mode = pose_rank - 1`, which assumes the exported poses are
        # in mode order. `export_nac_poses` sorts NAC-viable first by energy and
        # then the rest by energy, so that ordering is not guaranteed by
        # construction -- it happens to hold for all 1,751 poses on disk today,
        # which is exactly the kind of invariant that holds until it does not.
        # The SDF carries an explicit `mode` property; read it (#53).
        mode = pose_mode(pose, args.pose_rank)
        # `ident` carries the MODE, so two modes of one molecule are two rows
        # rather than one overwriting the other downstream.
        #
        # ALWAYS `_m<mode>`, INCLUDING MODE 0. It used to write the bare ident
        # for pose_rank 1, so mode 0 was `t4_x` in this table and `t4_x_m0` in
        # the rank table. Every obvious join on `ident` therefore dropped exactly
        # the modes that were simulated, which is how #53 stayed invisible and
        # how a wrong correlation reached #36. Join on (parent_ident, mode).
        rec = {"ident": f"{cand}_m{mode}" if mode is not None else cand,
               "parent_ident": cand, "pose_rank": args.pose_rank,
               "mode": mode, "sweep_ps": args.sweep_ps}
        if not pose.is_file():
            rec["status"] = f"no pose at {pose}"
            rows.append(rec); log.warning("%s: %s", cand, rec["status"]); continue

        if args.stage0_only:
            # STAGE 0 NOW ACTUALLY MEASURES SOMETHING (2026-09-02).
            #
            # It used to write this row and `continue` -- no geometry at all --
            # while the module docstring advertised "the elevated pose's OWN
            # geometry, before any simulation ... a filter that costs nothing".
            # The flag existed, produced a row, computed nothing, and those rows
            # were counted as finished sweeps by anything resuming on ident.
            #
            # AND IT DOES NOT PREDICT THE SWEEP. Measured over the first five
            # nac_v8 sweeps, all five had a DOCKED warhead-SG distance of
            # 2.78-2.97 A and outcomes from 0.000 to 0.926 attack-ready. What
            # separated them was the distance AFTER the 300 ps unrestrained
            # equilibration (3.58 A -> 0.926, 6.46 A -> 0.000), which is the
            # tier-1 quantity and is not obtainable without running it. So this
            # is honest bookkeeping, not a triage filter -- do not use it as one.
            rec["status"] = "stage0 only"
            try:
                c = _cand_meta(cand)
                m0 = _pose_at_rank(pose, args.pose_rank)
                mt = m0.GetSubstructMatches(Chem.MolFromSmarts(c["smarts"]))[0]
                pos = m0.GetConformer().GetPositions()[list(mt)]
                r0 = nac.measure(c["mechanism"], pos, _static_sg())
                rec.update({"start_dist_a": r0.distance,
                            "start_angle_deg": r0.angle,
                            "start_attack_ready": bool(r0.viable),
                            "stage0_frame": "docked pose vs static Cys113 SG"})
            except Exception as exc:                      # noqa: BLE001
                rec["status"] = f"stage0 failed: {type(exc).__name__}: {exc}"[:180]
            rows.append(rec); continue

        try:
            rep = run_sweep(cand, pose, args.pose_rank, args.gpu, args.sweep_ps,
                            None, abort_above=(args.abort_above_a or None))
        except SweepAborted as exc:
            # A RESULT, NOT A FAILURE. The pose was measured and had gone, so
            # the row carries the distance and a status of its own. `frac_*`
            # stay NaN because nothing was simulated -- writing 0.0 would be a
            # measurement nobody made.
            rec["status"] = "aborted: left during equilibration"
            rec["equil_dist_a"] = exc.dist_a
            rows.append(rec)
            log.info("%s: %s", cand, exc)
            continue
        except SweepError as exc:
            # Recorded, never dropped, and never as a number: the row says why
            # there is no reading rather than carrying a value from a partial run.
            rec["status"] = f"skipped: {exc}"
            rows.append(rec); log.warning("%s: %s", cand, exc); continue
        if rep is None:
            rec["status"] = "sweep failed"
            rows.append(rec); continue

        dense = rep / "sweep_dense.pdb"
        if not dense.is_file():
            mov.build_movie_pdb(rep, dense, n_frames=FRAMES)
        s = mp.nac_series(cand, rep, dense)
        if s is None:
            rec["status"] = "no attack-geometry series"
            rows.append(rec); continue

        rec["status"] = "ok"
        rec["mechanism"] = s["mechanism"]
        # The frame spacing is DERIVED from the run, not assumed: the movie is
        # built with n_frames=FRAMES over `ps` picoseconds, and `n_visits` is
        # meaningless without it (see MIN_DWELL_PS). Taking the actual series
        # length rather than FRAMES, because build_movie_pdb can return fewer
        # frames than asked for and a stale divisor would silently rescale every
        # dwell time.
        n_fr = max(1, len(s["dist"]))
        rec.update(geometry_stats(s["dist"], s["angle"], s["kind"],
                                  frame_ps=float(args.sweep_ps) / n_fr))
        # POSE STABILITY, from the `rmsd.xvg` this sweep already wrote. Half of
        # the 100 ns gate, and it costs a file read -- computing it here rather
        # than at elevation time means the verdict travels with the row instead
        # of being re-derived by whoever reads it next.
        rec.update(rmsd_stats(rep))
        rec.update(elevation_verdict(rec))
        rows.append(rec)
        log.info("%s: engaged %.1f%% within %.1f A, rmsd max %s mean %s -> %s",
                 cand, rec["frac_attack_ready"] * 100, attack_ready_max_a(),
                 f"{rec['rmsd_max_a']:.2f}" if "rmsd_max_a" in rec else "?",
                 f"{rec['rmsd_mean_a']:.2f}" if "rmsd_mean_a" in rec else "?",
                 "ELEVATE" if rec.get("elevate") else rec.get("elevate_why"))

    t = pd.DataFrame(rows)
    ok = t[t.status == "ok"] if "status" in t.columns else t
    if not ok.empty:
        t = t.sort_values("frac_attack_ready", ascending=False, na_position="last")
    dest = OUT.write("attack_sweep", ".csv")
    t.to_csv(dest, index=False)

    print(f"\nAttack-geometry sweep — {args.sweep_ps/1000:.0f} ns, ranked\n")
    if ok.empty:
        print("  no molecules completed"); print(f"\n  -> {dest}"); return
    cols = ["ident", "start_attack_ready", "frac_in_window",
            "frac_attack_ready", "n_visits", "min_dist_a"]
    print(t[[c for c in cols if c in t.columns]].round(4).to_string(index=False))
    # REJECT-ONLY, NOT A RANKING (#34).
    #
    # Agreement between the 10 ns sweep and the 100 ns run is rho = +0.60 on
    # `frac_attack_ready` -- the quantity printed here -- against +0.83 for mere
    # proximity, which is the number older docs quote. The pre-registered reading
    # table licenses +0.60 for discarding the bottom and explicitly NOT for
    # ordering the middle, and this printout used to announce the opposite.
    #
    # The order below is still by attack-readiness, because something has to be
    # printed in some order and the elevation queue has to take a best-first
    # view when capacity binds. What has changed is the claim attached to it.
    print("\n  ORDERED, NOT RANKED. rho(10 ns, 100 ns) = +0.60 on this reading —")
    print("  enough to reject the bottom, NOT enough to order the middle. Treat")
    print("  the ordering as a queue, not as a statement that #2 beats #5.")
    print("  A survivor is a molecule with a SUSTAINED episode (n_visits > 0);")
    print("  a single 20 ps touch is not evidence of anything.")
    print(f"\n  readings fixed in advance: docs/prereg_attack_sweep.md")
    print(f"  -> {dest}")


if __name__ == "__main__":
    main()
