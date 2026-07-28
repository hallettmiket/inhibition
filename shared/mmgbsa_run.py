"""
Purpose: The non-covalent MM-GBSA run, shared by T_1 and T_2.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: an approach's ranked frame (with `shortlist`) + its docked poses
Output: the frame with dG columns; per-candidate Amber working directories

WHY SHARED. T_1 and T_2 are the two reversible approaches and the integration
phase pools them, so their dG must mean the same thing. Same three legs, same
implicit solvent, same minimiser, same protonation model. Two scripts that begin
identical do not stay that way — the argument that moved docking and ranking
into shared modules applies here unchanged.

PROTONATE BEFORE PARAMETERISING. Generators emit neutral SMILES; at pH 7.4 a
carboxylic acid is a carboxylate, and Pin1's pocket is cationic (Arg69 6.5 A,
Lys63 7.3 A from Cys113). antechamber assigns AM1-BCC charges under the total it
is GIVEN, so a carboxylate passed as neutral redistributes a whole electron
across the molecule rather than merely losing one.

THE POSES WERE DOCKED NEUTRAL. Protonation moves hydrogens only, so the docked
heavy-atom geometry stays usable — but Vina scored the neutral form, and that is
recorded in the manifest rather than quietly corrected.

WORKERS SHARE ONE PROCESS. Candidates are scored in a pool and ONE process
writes the frame at the end, so parallel work cannot race on frame versions the
way a multi-process fan-out does.
"""

from __future__ import annotations

import json
import logging
import os
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from . import io as dio
from . import mmgbsa as mg
from . import mmgbsa_noncovalent as mgn
from . import protonation as prot

log = logging.getLogger(__name__)

DATA_ROOT = Path("/data/lab_vm/append_only/inhibition")
NICE = 19

# ONE THREAD PER WORKER. antechamber's AM1-BCC step shells out to `sqm`, which is
# OpenMP/MKL-threaded and, left alone, takes every core it can see — a single sqm
# was measured at 7,543% CPU (≈75 cores), and six "workers" consumed 204 of the
# machine's 224 cores while other researchers were logged in. Concurrency here is
# the worker count, NOT the thread count, so each worker is pinned to one thread
# and `workers` means what it says.
_SINGLE_THREAD_ENV = {
    "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


def _pin_to_one_thread() -> None:
    """Called at the top of every worker, before any Amber tool is invoked."""
    os.environ.update(_SINGLE_THREAD_ENV)


def score_one(job: dict) -> dict:
    """Build, minimise and score one candidate's three legs."""
    os.nice(NICE)
    _pin_to_one_thread()
    cid = job["candidate_id"]
    wd = Path(job["work_root"]) / cid
    wd.mkdir(parents=True, exist_ok=True)

    done = wd / "result.json"
    if done.is_file():
        return json.loads(done.read_text())

    try:
        pose = Path(job["pose_root"]) / f"{cid}_out.pdbqt"
        if not pose.is_file():
            pose = Path(job["pose_root"]) / f"{cid}.pdbqt"
        if not pose.is_file():
            raise mg.MMGBSAError(f"no docked pose for {cid} under {job['pose_root']}")

        sdf = mgn.pose_to_sdf(pose, wd, job["smiles"])
        # prepare_receptor emits the apo CYS form beside the CYX one, so the
        # covalent and non-covalent paths share one preparation and cannot drift
        # on histidine states.
        _cyx, cys, _idx, _n = mg.prepare_receptor(wd)
        mol2, frcmod = mgn.parameterize_ligand(sdf, wd,
                                               net_charge=int(job["formal_charge"]))
        legs = mgn.build_topologies(wd, mol2, frcmod, cys)
        verified = mgn.verify_complex(legs)

        energies = {leg: mg.minimize_and_score(wd, leg)
                    for leg in ("complex", "receptor", "ligand")}
        result = {"candidate_id": cid, **mgn.delta_g(energies),
                  **{f"topology_{k}": v for k, v in verified.items()}}
        done.write_text(json.dumps(result, indent=2), encoding="utf-8")
        return result
    except Exception as exc:  # noqa: BLE001 - one failure must not end the run
        (wd / "traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return {"candidate_id": cid, "mmgbsa_error": str(exc)[:300]}


def run(*, experiment: str, approach: str, work_dirname: str = "mmgbsa_2",
        workers: int = 6, limit: int | None = None):
    """Score one approach's shortlist and merge the result onto its frame."""
    data = DATA_ROOT / experiment
    work_root = data / work_dirname
    pose_root = data / "docking" / "poses"
    work_root.mkdir(parents=True, exist_ok=True)

    df, frame_path = dio.latest_frame(experiment, approach)
    if "shortlist" not in df.columns:
        raise SystemExit("frame has no shortlist — run the ranking stage first")

    todo = df[df["shortlist"].fillna(False)].copy()
    if limit:
        todo = todo.head(limit)
    log.info("[%s] MM-GBSA on %d shortlisted candidates from %s (%d workers)",
             approach, len(todo), frame_path.name, workers)

    todo = prot.protonate_frame(todo)
    changed = int(todo["charge_changed"].sum())
    unsure = int((~todo["protonation_confident"]).sum())
    log.info("[%s] protonation at pH %s: %d/%d change charge; %d carry a "
             "borderline group and are flagged", approach, prot.PHYSIOLOGICAL_PH,
             changed, len(todo), unsure)

    jobs = [{"candidate_id": r["candidate_id"],
             "smiles": r["protonated_smiles"],
             "formal_charge": int(r["protonated_charge"]),
             "work_root": str(work_root), "pose_root": str(pose_root)}
            for _, r in todo.iterrows()]

    results, failures = [], []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(score_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            (failures if "mmgbsa_error" in r else results).append(r)
            log.info("[%s] [%d/%d] %s -> %s", approach, i, len(jobs),
                     r["candidate_id"],
                     f"dG {r['dG_kcal']:.2f}" if "dG_kcal" in r
                     else f"FAILED {r['mmgbsa_error'][:80]}")

    if not results and not failures:
        raise SystemExit("nothing to do")

    cols = pd.DataFrame(results + failures)
    keep = [c for c in ("candidate_id", "dG_kcal", "G_complex", "G_receptor",
                        "G_ligand", "mmgbsa_error") if c in cols.columns]
    stale = [c for c in keep if c != "candidate_id" and c in df.columns]
    if stale:
        df = df.drop(columns=stale)
    merged = df.merge(cols[keep].drop_duplicates("candidate_id"),
                      on="candidate_id", how="left")
    if len(merged) != len(df):
        raise RuntimeError(f"merge changed row count {len(df)} -> {len(merged)}")

    if limit:
        log.warning("--limit %d: NOT writing a frame; a partial run must not "
                    "become the latest frame for the next stage to read.", limit)
        return merged, None, results, failures, changed

    out = dio.write_full_frame(
        merged, approach=approach, experiment=experiment,
        stage=f"{approach}_mmgbsa",
        params={"tier": "T5 physics rescoring",
                "scheme": "non-covalent 3-leg",
                "igb": mg.IGB, "pb_radii": mg.PB_RADII,
                "ensemble_averaged": False,
                "uncertainty_reported": False,
                "protonation_ph": prot.PHYSIOLOGICAL_PH,
                "n_charge_changed": changed,
                "docking_used_neutral_forms": True,
                "not_implemented": ["explicit-solvent MD", "AI cofold pose",
                                    "anti-target / selectivity docking"],
                "comparable": "across the whole approach (no link-atom constant)",
                "n_scored": len(results), "n_failed": len(failures)},
        inputs={"frame": frame_path})
    return merged, out, results, failures, changed


# --------------------------------------------------------------------------
# covalent
# --------------------------------------------------------------------------

def score_one_covalent(job: dict) -> dict:
    """Build, minimise and score one covalent adduct's three legs."""
    os.nice(NICE)
    _pin_to_one_thread()
    from . import covalent_adduct as cad
    from . import warhead_library as wl

    cid, did = job["candidate_id"], job["dock_id"]
    wd = Path(job["work_root"]) / cid
    wd.mkdir(parents=True, exist_ok=True)

    done = wd / "result.json"
    if done.is_file():
        cached = json.loads(done.read_text())
        cached.setdefault("dock_id", did)
        return cached

    try:
        pose = Path(job["pose_root"]) / f"{did}_docked.sdf"
        if not pose.is_file():
            raise mg.MMGBSAError(f"no docked pose at {pose}")

        lib = wl.load()
        smarts = cad.adduct_attachment_smarts(job["warhead_class"], library=lib)
        cyx, cys, cyx_idx, n_res = mg.prepare_receptor(wd)
        mol2, frcmod, att, cap, q = mg.parameterize_ligand(
            pose, wd, smarts, net_charge=0)
        legs = mg.build_topologies(wd, mol2, frcmod, cyx, cys, cyx_idx,
                                   n_res + 1, att, cap, q)
        verified = mg.verify_complex(legs["complex"][0], cyx_idx, att)

        energies = {leg: mg.minimize_and_score(wd, leg)
                    for leg in ("complex", "receptor", "ligand")}
        result = {"candidate_id": cid, "dock_id": did,
                  "warhead_class": job["warhead_class"],
                  **mg.delta_g(energies),
                  **{f"topology_{k}": v for k, v in verified.items()}}
        done.write_text(json.dumps(result, indent=2), encoding="utf-8")
        return result
    except Exception as exc:  # noqa: BLE001
        (wd / "traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return {"candidate_id": cid, "dock_id": did,
                "warhead_class": job["warhead_class"],
                "mmgbsa_error": str(exc)[:300]}


def run_covalent(*, experiment: str, approach: str, work_dirname: str = "mmgbsa",
                 workers: int = 6, limit: int | None = None,
                 classes: set[str] | None = None):
    """Score one covalent approach's shortlist, one system per MOLECULE.

    ONE SYSTEM PER MOLECULE, NOT PER ROUTE (D0029). Several warhead classes can
    reach an identical adduct, so a shortlist may carry the same molecule more
    than once. Scoring each row minimised the SAME ~2,400-atom system
    repeatedly; results are computed per `dock_id` and mapped back to every
    route afterwards.

    NO PROTONATION STEP HERE, deliberately. The covalent ligand is already
    bonded to Cys113 and its parameterisation runs through the link-atom scheme,
    which sets the attachment charge explicitly. Running an ionizer over an
    adduct would perturb exactly the atom the junction parameters describe.
    """
    data = DATA_ROOT / experiment
    work_root = data / work_dirname
    pose_root = data / "docking"
    work_root.mkdir(parents=True, exist_ok=True)

    df, frame_path = dio.latest_frame(experiment, approach)
    for col in ("shortlist", "dock_id"):
        if col not in df.columns:
            raise SystemExit(f"frame has no {col!r} — run the earlier stages first")

    todo = df[df["shortlist"].fillna(False)].copy()
    if classes:
        todo = todo[todo["warhead_class"].isin(classes)]
    n_rows = len(todo)
    todo = todo.drop_duplicates("dock_id")
    if len(todo) != n_rows:
        log.info("[%s] %d shortlisted rows -> %d distinct molecules; %d duplicate "
                 "route(s) reuse their molecule's result (D0029)",
                 approach, n_rows, len(todo), n_rows - len(todo))
    if limit:
        todo = todo.head(limit)
    log.info("[%s] covalent MM-GBSA on %d molecules from %s (%d workers)",
             approach, len(todo), frame_path.name, workers)

    jobs = [{"candidate_id": r["candidate_id"], "dock_id": r["dock_id"],
             "warhead_class": r["warhead_class"],
             "work_root": str(work_root), "pose_root": str(pose_root)}
            for _, r in todo.iterrows()]

    results, failures = [], []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(score_one_covalent, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            (failures if "mmgbsa_error" in r else results).append(r)
            log.info("[%s] [%d/%d] %s (%s) -> %s", approach, i, len(jobs),
                     r["candidate_id"], r["warhead_class"],
                     f"dG {r['dG_kcal']:.2f}" if "dG_kcal" in r
                     else f"FAILED {r['mmgbsa_error'][:70]}")

    if not results and not failures:
        raise SystemExit("nothing to do")

    cols = pd.DataFrame(results + failures)
    keep = [c for c in ("dock_id", "dG_kcal", "G_complex", "G_receptor",
                        "G_ligand", "mmgbsa_error") if c in cols.columns]
    stale = [c for c in keep if c != "dock_id" and c in df.columns]
    if stale:
        df = df.drop(columns=stale)
    merged = df.merge(cols[keep].drop_duplicates("dock_id"),
                      on="dock_id", how="left")
    if len(merged) != len(df):
        raise RuntimeError(f"merge changed row count {len(df)} -> {len(merged)}")

    if limit or classes:
        log.warning("partial run (limit/classes): NOT writing a frame; it would "
                    "hold only part of the shortlist and become the latest.")
        return merged, None, results, failures

    out = dio.write_full_frame(
        merged, approach=approach, experiment=experiment,
        stage=f"{approach}_mmgbsa",
        params={"tier": "T5 physics rescoring",
                "scheme": "link-atom 3-leg, cut at Cys113 SG-C",
                "igb": mg.IGB, "pb_radii": mg.PB_RADII,
                "ensemble_averaged": False,
                "uncertainty_reported": False,
                "one_system_per_molecule": True,
                "comparable": "within warhead class only (D0020, D0023)",
                "n_scored": len(results), "n_failed": len(failures)},
        inputs={"frame": frame_path})
    return merged, out, results, failures
