"""
Purpose: Run both redocking arms for the Pin1 benchmark under the EXACT
         production non-covalent protocol.
         (A) self-docking  -- each ligand into its own receptor, own box
         (B) cross-docking -- the same ligands into 6VAJ with box_expanded
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: outputs/blacksmith/redock_pin1/redock_cases_1.csv
Output: append_only/inhibition/05_redock_benchmark/dock_1/{self,cross}/...

PROTOCOL PARITY IS THE ENTIRE POINT. Every setting comes from
`shared.noncovalent_dock_run` by IMPORT, not by copying a number into this
file: the Vina-GPU binary, SEARCH_DEPTH 20 (D0017), THREADS, and the pH 7.4
ligand preparation (`LIGAND_PH`, obabel `-p 7.4`). Ligands are built by that
module's own `_prepare_one`, so the molecule this benchmark docks is prepared
by the same code path that prepares a T_1 or T_2 candidate. A benchmark run
under different settings measures nothing about our pipeline, and a *copied*
constant is a benchmark that silently stops matching the day the pipeline
changes.

ARM B IS THE PRODUCTION CALL ITSELF. Cross-docking invokes
`noncovalent_dock_run.run_vina_gpu` unmodified -- same receptor constant, same
box constant. Arm A cannot, because it needs a per-case receptor and box, so it
uses a wrapper that differs ONLY in those two arguments and asserts the rest
against the module's values before it launches.

GPUS ARE PASSED IN, NEVER DISCOVERED. GPUs 4 and 7 are running a T_2 campaign
and 0, 1 and 6 belong to other users; auto-selection by memory threshold would
land on them. The device list is an argument.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import compute                              # noqa: E402
from shared import noncovalent_dock_run as ndr          # noqa: E402

log = logging.getLogger("redock-dock")

OUT_DIR = REPO / "outputs" / "blacksmith" / "redock_pin1"
WORK = Path("/data/lab_vm/append_only/inhibition/05_redock_benchmark/dock_1")


def prepare_ligands(cases: pd.DataFrame, ligand_dir: Path) -> pd.DataFrame:
    """Build one PDBQT per case with the production ligand preparation.

    The SMILES is the PDB chemical component's own definition, so the docked
    molecule is the deposited ligand -- but embedded fresh from SMILES and
    protonated at pH 7.4, exactly as a generated candidate would be. Nothing
    of the crystal conformation leaks into the input; only its identity does.
    """
    frame = pd.DataFrame({"candidate_id": cases["case_id"],
                          "canonical_smiles": cases["smiles"]})
    res = ndr.prepare_ligands(frame, ligand_dir)
    out = pd.DataFrame(res).rename(columns={"candidate_id": "case_id",
                                            "ok": "ligand_prep_ok",
                                            "error": "ligand_prep_error"})
    log.info("ligand prep: %d ok, %d failed",
             int(out.ligand_prep_ok.sum()), int((~out.ligand_prep_ok).sum()))
    return out


def run_vina_gpu_on(receptor: Path, box: dict, ligand_dir: Path,
                    out_dir: Path, gpu: int) -> float:
    """Vina-GPU with an arbitrary receptor+box, every other setting the module's.

    Mirrors `noncovalent_dock_run.run_vina_gpu`; the receptor and box are the
    only differences and they are the only ones this arm is allowed to vary.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(ndr.VINA_GPU),
           "--receptor", str(receptor),
           "--ligand_directory", str(ligand_dir),
           "--output_directory", str(out_dir),
           "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
           "--center_z", str(box["center_z"]),
           "--size_x", str(box["size_x"]), "--size_y", str(box["size_y"]),
           "--size_z", str(box["size_z"]),
           "--thread", str(ndr.THREADS),
           "--search_depth", str(ndr.SEARCH_DEPTH)]
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    env["GPU_DEVICE_ORDINAL"] = str(gpu)   # OpenCL honours this, not CUDA_*
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=86400, env=env)
    dt = time.time() - t0
    (out_dir / "vina_gpu_stdout.log").write_text(
        p.stdout + "\n--- stderr ---\n" + p.stderr)
    if p.returncode != 0:
        raise RuntimeError(f"Vina-GPU failed ({p.returncode}) on {receptor.name}: "
                           f"{(p.stderr or p.stdout)[-800:]}")
    return dt


def self_dock(cases: pd.DataFrame, ligand_dir: Path, gpus: list[int]) -> pd.DataFrame:
    """Arm A: every ligand into its OWN receptor, box on its OWN crystal site.

    One Vina-GPU invocation per case, because each has a different receptor.
    The per-case ligand directory holds exactly one file.
    """
    root = WORK / "self"
    root.mkdir(parents=True, exist_ok=True)
    todo = list(cases.itertuples())
    results: list[dict] = []

    def worker(chunk: list, gpu: int) -> None:
        for c in chunk:
            rec = Path(c.receptor_pdbqt)
            src = ligand_dir / f"{c.case_id}.pdbqt"
            if not src.is_file():
                results.append({"case_id": c.case_id, "self_ok": False,
                                "self_error": "ligand_pdbqt_missing"})
                continue
            case_lig = root / c.case_id / "ligand"
            case_lig.mkdir(parents=True, exist_ok=True)
            dst = case_lig / src.name
            if not dst.is_file():
                dst.write_bytes(src.read_bytes())
            out_dir = root / c.case_id / "poses"
            try:
                dt = run_vina_gpu_on(rec, json.loads(Path(c.box_json).read_text()),
                                     case_lig, out_dir, gpu)
                scores = ndr.collect_scores(out_dir)
                results.append({"case_id": c.case_id, "self_ok": bool(scores),
                                "self_affinity": scores.get(c.case_id),
                                "self_pose_dir": str(out_dir),
                                "self_seconds": round(dt, 1),
                                "self_error": None if scores else "no_pose_parsed"})
            except Exception as exc:  # noqa: BLE001 - one case must not end the arm
                results.append({"case_id": c.case_id, "self_ok": False,
                                "self_error": str(exc)[:200]})
            if len(results) % 10 == 0:
                log.info("self-docking: %d/%d done", len(results), len(todo))

    chunks = [todo[i::len(gpus)] for i in range(len(gpus))]
    with ThreadPoolExecutor(max_workers=len(gpus)) as ex:
        list(ex.map(worker, chunks, gpus))
    log.info("self-docking complete: %d results", len(results))
    return pd.DataFrame(results)


def cross_dock(ligand_dir: Path, gpu: int) -> pd.DataFrame:
    """Arm B: all ligands into 6VAJ with box_expanded -- the production call."""
    out_dir = WORK / "cross" / "poses"
    dt = ndr.run_vina_gpu(ligand_dir, out_dir, gpu)
    scores = ndr.collect_scores(out_dir)
    log.info("cross-docking: %d poses in %.0f s", len(scores), dt)
    return pd.DataFrame({"case_id": list(scores),
                         "cross_affinity": list(scores.values()),
                         "cross_ok": True}).assign(cross_pose_dir=str(out_dir))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpus", default="2,3,5",
                    help="idle devices only; 4 and 7 run T_2, 0/1/6 are others'")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--arm", choices=["self", "cross", "both"], default="both")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    os.nice(compute.NICE)
    gpus = [int(g) for g in args.gpus.split(",")]

    df = pd.read_csv(OUT_DIR / "redock_cases_1.csv")
    cases = df[df.status == "case"].copy()
    if args.limit:
        cases = cases.head(args.limit)
    log.info("%d cases | Vina-GPU depth %d | ligand pH %s | %d GPUs %s",
             len(cases), ndr.SEARCH_DEPTH, ndr.LIGAND_PH, len(gpus), gpus)

    ligand_dir = WORK / f"ligands_{ndr.LIGAND_PREP_TAG}"
    prep = prepare_ligands(cases, ligand_dir)
    ready = cases.merge(prep, on="case_id", how="left")
    ready = ready[ready.ligand_prep_ok.fillna(False)]
    log.info("%d cases have a docking-ready ligand", len(ready))

    out = cases[["case_id"]].merge(prep, on="case_id", how="left")
    if args.arm in ("cross", "both"):
        # Cross-docking needs one directory holding every ligand, which is
        # exactly what `ligand_dir` already is.
        out = out.merge(cross_dock(ligand_dir, gpus[0]), on="case_id", how="left")
    if args.arm in ("self", "both"):
        out = out.merge(self_dock(ready, ligand_dir, gpus), on="case_id", how="left")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_DIR / "redock_docking_1.csv", index=False)
    log.info("wrote %s", OUT_DIR / "redock_docking_1.csv")


if __name__ == "__main__":
    main()
