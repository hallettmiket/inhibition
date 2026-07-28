"""
Purpose: The non-covalent docking run, shared by T_1 and T_2.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: an approach's latest frame (survivors carry `canonical_smiles`)
Output: the frame with `vina_affinity` merged back; ligand + pose PDBQTs

THE NON-COVALENT COUNTERPART TO `covalent_dock_run`. T_1 and T_2 are the two
reversible approaches and the integration phase pools them on one plot, so they
must dock through one engine, one box and one search depth for the same reason
T_3 and T_4 must. The loop lives here; each approach supplies only its identity.

ENGINE: VINA-GPU AT search_depth 20 (D0017). Not CPU Vina, and not gnina. D0017
adopted Vina-GPU only after a search_depth sweep showed the depth-10
disagreement with CPU Vina was search CONVERGENCE rather than implementation —
every metric improved monotonically with depth, and at 20 the enrichment ROC-AUC
drift is 0.005 against a 0.10 threshold. Below 20 the adoption evidence does not
hold, so the depth is pinned here rather than exposed as a tuning knob.

THE GATE ONLY TRANSFERS IF THE ENGINE MATCHES. The non-covalent enrichment gate
measured ROC-AUC 0.535 (D0016) using this receptor and this box. That verdict
describes T_1 and T_2 only insofar as they dock the same way, which is why the
expanded box is read from the shared receptor directory rather than restated.

WHY THE EXPANDED BOX. `box_expanded.json` (26 A) is declared `used_by: [t1, t2]`
against the covalent box's 20 A. T_1 and T_2 place whole molecules with no
anchor at Cys113, so they need room the covalent approaches do not.

GPUS ARE ALLOCATED EXPLICITLY, NOT INFERRED. Vina-GPU runs in virtual-screening
mode: one process, one GPU, a whole ligand directory. The caller passes the
device id. This is deliberate — `covalent_dock_run.select_gpus` infers idleness
from a 1024 MiB memory threshold, and gnina occupies only ~500 MiB, so a job
that auto-selected would land on GPUs already running a covalent dock while
idle ones sat unused.
"""

from __future__ import annotations

import json
import logging
import os
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from . import io as dio

RDLogger.DisableLog("rdApp.*")

log = logging.getLogger(__name__)

DATA_ROOT = Path("/data/lab_vm/append_only/inhibition")
RECEPTOR_ROOT = Path("/data/lab_vm/immutable/inhibition/receptor")
RECEPTOR_PDBQT = RECEPTOR_ROOT / "6VAJ_prepared.pdbqt"
BOX_EXPANDED = RECEPTOR_ROOT / "box_expanded.json"
VINA_GPU = Path("/data/lab_vm/envs/dwi_vinagpu/bin/vina-gpu")
OBABEL = "/data/lab_vm/envs/dwi_cheminf/bin/obabel"

SEARCH_DEPTH = 20      # D0017 — the adoption evidence does not hold below this
THREADS = 8000
NICE = 19
# Ligand prep is embarrassingly parallel and would happily take all 224 cores.
# The project budget lives in shared/compute.py (raised to 50 by @mhallet on
# 2026-07-28).
from . import compute                             # noqa: E402
MAX_PREP_WORKERS = compute.MAX_CPU_WORKERS

_RESULT_RE = re.compile(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)")


def _prepare_one(job: dict) -> dict:
    """Embed one ligand and convert it to PDBQT. Runs in a worker process."""
    out = Path(job["ligand_dir"]) / f"{job['candidate_id']}.pdbqt"
    if out.is_file():
        return {"candidate_id": job["candidate_id"], "ok": True}
    sdf = out.with_suffix(".sdf")
    try:
        m = Chem.MolFromSmiles(job["smiles"])
        if m is None:
            return {"candidate_id": job["candidate_id"], "ok": False,
                    "error": "unparseable"}
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
            return {"candidate_id": job["candidate_id"], "ok": False,
                    "error": "embed_failed"}
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        except Exception:  # noqa: BLE001 - MMFF can fail; the pose is still usable
            pass
        w = Chem.SDWriter(str(sdf))
        w.write(m)
        w.close()
        conv = subprocess.run([OBABEL, str(sdf), "-O", str(out)],
                              capture_output=True, text=True, timeout=120)
        sdf.unlink(missing_ok=True)
        if conv.returncode != 0 or not out.is_file():
            return {"candidate_id": job["candidate_id"], "ok": False,
                    "error": "pdbqt_conversion_failed"}
        return {"candidate_id": job["candidate_id"], "ok": True}
    except Exception as exc:  # noqa: BLE001 - one bad ligand must not end prep
        return {"candidate_id": job["candidate_id"], "ok": False,
                "error": str(exc)[:200]}


def prepare_ligands(survivors: pd.DataFrame, ligand_dir: Path,
                    smiles_col: str = "canonical_smiles") -> list[dict]:
    """Embed + convert every survivor, capped at the lab's CPU limit."""
    ligand_dir.mkdir(parents=True, exist_ok=True)
    jobs = [{"candidate_id": r["candidate_id"], "smiles": r[smiles_col],
             "ligand_dir": str(ligand_dir)}
            for _, r in survivors.iterrows()]
    results = []
    with ProcessPoolExecutor(max_workers=MAX_PREP_WORKERS) as ex:
        for i, res in enumerate(ex.map(_prepare_one, jobs, chunksize=16), 1):
            results.append(res)
            if i % 500 == 0:
                log.info("prepared %d/%d ligands", i, len(jobs))
    bad = [r for r in results if not r["ok"]]
    if bad:
        log.warning("%d ligand(s) failed preparation and will not be docked; "
                    "first: %s", len(bad), bad[0])
    return results


def run_vina_gpu(ligand_dir: Path, out_dir: Path, gpu: int) -> float:
    """Vina-GPU in virtual-screening mode on ONE explicitly chosen GPU."""
    out_dir.mkdir(parents=True, exist_ok=True)
    b = json.loads(BOX_EXPANDED.read_text())
    cmd = [str(VINA_GPU),
           "--receptor", str(RECEPTOR_PDBQT),
           "--ligand_directory", str(ligand_dir),
           "--output_directory", str(out_dir),
           "--center_x", str(b["center_x"]), "--center_y", str(b["center_y"]),
           "--center_z", str(b["center_z"]),
           "--size_x", str(b["size_x"]), "--size_y", str(b["size_y"]),
           "--size_z", str(b["size_z"]),
           "--thread", str(THREADS), "--search_depth", str(SEARCH_DEPTH)]
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    env["GPU_DEVICE_ORDINAL"] = str(gpu)      # OpenCL honours this, not CUDA_*
    log.info("Vina-GPU on GPU %d, depth %d, %d ligands",
             gpu, SEARCH_DEPTH, len(list(ligand_dir.glob("*.pdbqt"))))
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=86400, env=env)
    dt = time.time() - t0
    (out_dir / "vina_gpu_stdout.log").write_text(
        p.stdout + "\n--- stderr ---\n" + p.stderr)
    if p.returncode != 0:
        raise RuntimeError(f"Vina-GPU failed ({p.returncode}); see "
                           f"{out_dir/'vina_gpu_stdout.log'}\n"
                           f"{(p.stderr or p.stdout)[-1500:]}")
    log.info("Vina-GPU finished in %.1f s", dt)
    return dt


def collect_scores(out_dir: Path) -> dict[str, float]:
    """Best score per ligand from Vina-GPU's output PDBQTs."""
    scores: dict[str, float] = {}
    for f in out_dir.glob("*.pdbqt"):
        m = _RESULT_RE.search(f.read_text(errors="replace"))
        if m:
            scores[f.stem.replace("_out", "")] = float(m.group(1))
    log.info("parsed %d Vina-GPU scores", len(scores))
    return scores


def run(*, experiment: str, approach: str, frame_prefix: str, gpu: int,
        limit: int | None = None):
    """Dock one approach's survivors and merge the result onto its frame."""
    os.nice(NICE)
    work = DATA_ROOT / experiment / "docking"
    ligand_dir, out_dir = work / "ligands", work / "poses"

    frame_path = dio.latest(DATA_ROOT / experiment, frame_prefix, ".parquet")
    if frame_path is None:
        raise SystemExit(f"no {frame_prefix} frame for {experiment}")
    df = dio.read_frame(frame_path)

    survivors = df[df["rejected_at"].isna()].copy()
    if limit:
        survivors = survivors.head(limit)
    log.info("[%s] %d survivors (of %d in the frame)",
             approach, len(survivors), len(df))

    prep = prepare_ligands(survivors, ligand_dir)
    n_ready = sum(1 for r in prep if r["ok"])
    log.info("[%s] %d/%d ligands ready to dock", approach, n_ready, len(prep))
    if not n_ready:
        raise SystemExit("no ligand survived preparation")

    elapsed = run_vina_gpu(ligand_dir, out_dir, gpu)
    scores = collect_scores(out_dir)

    # Drop a stale score column BEFORE merging, so a re-run does not produce
    # vina_affinity_x / _y and silently lose the column downstream.
    if "vina_affinity" in df.columns:
        log.info("dropping stale vina_affinity before merge")
        df = df.drop(columns=["vina_affinity"])
    scored = pd.DataFrame({"candidate_id": list(scores),
                           "vina_affinity": list(scores.values())})
    merged = df.merge(scored, on="candidate_id", how="left")
    if len(merged) != len(df):
        raise RuntimeError(
            f"merge changed row count {len(df)} -> {len(merged)}; duplicate "
            "candidate_id in the frame or the results")

    n_docked = int(merged["vina_affinity"].notna().sum())

    # See covalent_dock_run: a --limit run must not become the latest frame.
    if limit:
        log.warning("--limit %d: NOT writing a frame. A partial run must not "
                    "become the latest frame for the next stage to read.", limit)
        return merged, None, survivors, n_docked, elapsed

    out = dio.write_full_frame(
        merged, approach=approach, experiment=experiment,
        stage=f"{approach}_noncovalent_dock",
        params={"engine": "Vina-GPU 2.1 (D0017)",
                "search_depth": SEARCH_DEPTH,
                "threads": THREADS,
                "box": BOX_EXPANDED.name,
                "gpu": gpu,
                "rank_metric": "vina_affinity (kcal/mol, lower better)",
                "enrichment_caveat":
                    "D0016: non-covalent enrichment on this pocket is ROC-AUC "
                    "0.535 (CI 0.215-0.855, EF1% 0.0). These rankings are "
                    "weakly supported and must be shown with that verdict.",
                "elapsed_s": round(elapsed, 1),
                "n_docked": n_docked},
        inputs={"frame": frame_path, "receptor": RECEPTOR_PDBQT})
    return merged, out, survivors, n_docked, elapsed
