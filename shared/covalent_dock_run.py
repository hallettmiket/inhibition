"""
Purpose: The covalent docking run itself, shared by T_3 and T_4.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: an approach's latest frame (survivors carry `canonical_smiles` + `warhead_class`)
Output: the frame with dock columns merged back; poses + a resumable JSONL

WHY THIS IS SHARED, NOT COPIED (control S3). T_3 and T_4 are the two covalent
approaches and the integration phase offers a within-covalent comparison. That
comparison is only defensible if both docked through one protocol — and the
surest way to lose that is two scripts that start identical and drift. T_3's
config says as much about the adduct transform ("a second code path is exactly
how two approaches drift apart"); the same argument applies to the loop around
it, so the loop lives here and each approach supplies only its own identity.

WHAT EACH APPROACH STILL OWNS: the experiment directory, the frame prefix, and
which rows are survivors. Nothing about how a ligand becomes a pose.

THE ADDUCT IS WHAT GETS DOCKED (D0022, D0030). gnina replaces an implicit
hydrogen on the matched atom; it does not remove a leaving group and it does not
saturate an alkene. Every candidate is converted to its post-reaction form first
and docked through the SMARTS valid on that product.

ONE DOCK PER DISTINCT (ADDUCT, SMARTS). Candidates that collapse to the same
product are docked once and the result mapped back to every candidate sharing
it. For T_4 that merges the three SN2 acetamides (D0029); for T_3 it merges
duplicate decorations the generator proposed more than once.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from . import covalent_adduct as cad
from . import covalent_protocol as cp
from . import io as dio
from . import warhead_library as wl

RDLogger.DisableLog("rdApp.*")

log = logging.getLogger(__name__)

DATA_ROOT = Path("/data/lab_vm/append_only/inhibition")

MAX_GPUS_IDLE = 6
MAX_GPUS_SHARED = 4
GPU_BUSY_MEMORY_MIB = 1024
NICE = 19
CHUNK = 60          # ligands between GPU re-evaluations

DOCK_COLS = ("cnn_affinity", "cnn_score", "affinity_kcal",
             "cnn_uncalibrated_for_covalent", "protocol_fingerprint",
             "pose_path")

ADDUCT_COLS = ("candidate_id", "dock_id", "adduct_smiles",
               "adduct_attachment_smarts", "leaving_group_smiles",
               "adduct_atoms_removed", "adduct_approximation",
               "adduct_degenerate_attachment")


def select_gpus(explicit: list[int] | None = None) -> list[int]:
    """Idle GPUs, yielding to other researchers — re-checked between chunks.

    `explicit` overrides the search entirely. Prefer it whenever another docking
    job is already running: the auto-search calls a GPU busy above
    GPU_BUSY_MEMORY_MIB, and gnina occupies only ~500 MiB, so two auto-selecting
    covalent jobs both pick the SAME first-N devices and leave the rest idle.
    Vina-GPU (~18 GB) is visible to the threshold; gnina is not.
    """
    if explicit:
        log.info("GPUs: using %s (explicitly allocated)", explicit)
        return explicit
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=index,memory.used",
             "--format=csv,noheader,nounits"], text=True, timeout=30)
    except Exception as exc:  # noqa: BLE001
        log.warning("cannot query GPUs (%s); using GPU 0", exc)
        return [0]
    free, busy = [], []
    for line in out.strip().splitlines():
        idx, used = (x.strip() for x in line.split(","))
        (busy if int(used) > GPU_BUSY_MEMORY_MIB else free).append(int(idx))
    cap = MAX_GPUS_SHARED if busy else MAX_GPUS_IDLE
    chosen = free[:cap] or [0]
    log.info("GPUs: using %s (%d idle, %d busy elsewhere)",
             chosen, len(free), len(busy))
    return chosen


def _dock_one(job: dict) -> dict:
    """Embed and covalently dock one adduct on its assigned GPU."""
    os.environ["CUDA_VISIBLE_DEVICES"] = str(job["gpu"])
    wd = Path(job["workdir"])
    lig = wd / f"{job['dock_id']}.sdf"
    try:
        m = Chem.MolFromSmiles(job["adduct_smiles"])
        if m is None:
            return {**job, "error": "unparseable"}
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
            return {**job, "error": "embed_failed"}
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        except Exception:  # noqa: BLE001
            pass
        lig.parent.mkdir(parents=True, exist_ok=True)
        w = Chem.SDWriter(str(lig))
        w.write(m)
        w.close()
        r = cp.dock(lig, wd / f"{job['dock_id']}_docked.sdf",
                    job["warhead_class"], timeout_s=900)
        return {**job, **{k: r[k] for k in DOCK_COLS}}
    except Exception as exc:  # noqa: BLE001 - one bad ligand must not kill the run
        return {**job, "error": str(exc)[:300]}


def build_adducts(survivors: pd.DataFrame, lib) -> pd.DataFrame:
    """Post-reaction forms for every survivor, keyed by distinct product."""
    rows, failures = [], []
    for _, r in survivors.iterrows():
        try:
            a = cad.to_adduct_form(r["canonical_smiles"], r["warhead_class"],
                                   library=lib)
        except cad.AdductError as exc:
            failures.append((r["candidate_id"], str(exc)[:160]))
            continue
        rows.append({"candidate_id": r["candidate_id"],
                     "warhead_class": r["warhead_class"], **a.as_dict()})
    if failures:
        log.warning("%d candidate(s) have no well-defined adduct form and will "
                    "not be docked; first: %s", len(failures), failures[0])
    adducts = pd.DataFrame(rows)
    if adducts.empty:
        raise SystemExit("no candidate produced an adduct form")

    adducts["dock_id"] = adducts.apply(
        lambda r: "d_" + hashlib.sha256(
            f"{r['adduct_smiles']}|{r['adduct_attachment_smarts']}"
            .encode("utf-8")).hexdigest()[:12], axis=1)
    return adducts


def run(*, experiment: str, approach: str, frame_prefix: str,
        limit: int | None = None, results_name: str = "results_adduct.jsonl",
        gpus: list[int] | None = None):
    """Dock one approach's survivors and merge the result back onto its frame."""
    os.nice(NICE)
    out_root = DATA_ROOT / experiment / "docking"
    out_root.mkdir(parents=True, exist_ok=True)

    frame_path = dio.latest(DATA_ROOT / experiment, frame_prefix, ".parquet")
    if frame_path is None:
        raise SystemExit(f"no {frame_prefix} frame for {experiment}")
    df = dio.read_frame(frame_path)

    # Pin the protocol ONCE, up front. If it cannot be built, nothing should dock.
    proto = cp.load()
    log.info("covalent protocol fingerprint %s", proto.fingerprint()[:16])

    survivors = df[df["rejected_at"].isna()].copy()
    if limit:
        survivors = survivors.head(limit)
    log.info("[%s] %d survivors (of %d in the frame)",
             approach, len(survivors), len(df))

    lib = wl.load()
    adducts = build_adducts(survivors, lib)
    unique = adducts.drop_duplicates("dock_id")
    log.info("[%s] %d survivors -> %d distinct adduct docks (%.1fx saving)",
             approach, len(adducts), len(unique),
             len(adducts) / max(len(unique), 1))
    shared = (adducts.groupby("dock_id")["warhead_class"].nunique() > 1).sum()
    if shared:
        log.info("%d dock(s) are shared by more than one warhead class — their "
                 "adducts are the same molecule (D0029)", shared)

    results_path = out_root / results_name
    done: set[str] = set()
    if results_path.is_file():
        for line in results_path.read_text().splitlines():
            try:
                done.add(json.loads(line)["dock_id"])
            except Exception:  # noqa: BLE001 - older rows have no dock_id
                continue
        log.info("resuming: %d already docked", len(done))

    jobs = [{"dock_id": r["dock_id"],
             "adduct_smiles": r["adduct_smiles"],
             "warhead_class": r["warhead_class"],
             "workdir": str(out_root)}
            for _, r in unique.iterrows() if r["dock_id"] not in done]
    log.info("[%s] %d ligands to dock", approach, len(jobs))

    n = 0
    with open(results_path, "a", encoding="utf-8") as fh:
        for start in range(0, len(jobs), CHUNK):
            chunk = jobs[start:start + CHUNK]
            chosen = select_gpus(gpus)
            for i, j in enumerate(chunk):
                j["gpu"] = chosen[i % len(chosen)]
            with ProcessPoolExecutor(max_workers=len(chosen)) as ex:
                for fut in as_completed({ex.submit(_dock_one, j): j for j in chunk}):
                    fh.write(json.dumps(fut.result()) + "\n")
                    fh.flush()
                    n += 1
                    if n % 50 == 0:
                        log.info("[%s] %d/%d docked", approach, n, len(jobs))

    rows = [json.loads(l) for l in results_path.read_text().splitlines() if l.strip()]
    docked = pd.DataFrame(rows)
    if "error" in docked:
        errs = docked[docked["error"].notna()]
        if len(errs):
            log.warning("%d ligand(s) failed to dock; their rows keep empty "
                        "dock columns", len(errs))

    keep_dock = ["dock_id"] + [c for c in DOCK_COLS if c in docked.columns]
    merged_adducts = (adducts[[c for c in ADDUCT_COLS if c in adducts.columns]]
                      .merge(docked[keep_dock].drop_duplicates("dock_id"),
                             on="dock_id", how="left"))

    # Drop stale dock columns BEFORE merging. Re-running this stage on a frame
    # that had been merged once produced affinity_kcal_x / _y, so
    # `merged["affinity_kcal"]` did not exist and the success counter reported
    # 0 on a run that had actually worked.
    stale = [c for c in merged_adducts.columns
             if c != "candidate_id" and c in df.columns]
    if stale:
        log.info("dropping %d stale dock column(s) before merge: %s",
                 len(stale), stale)
        df = df.drop(columns=stale)
    merged = df.merge(merged_adducts.drop_duplicates("candidate_id"),
                      on="candidate_id", how="left")
    if len(merged) != len(df):
        raise RuntimeError(
            f"merge changed row count {len(df)} -> {len(merged)}; duplicate "
            "candidate_id in the frame or the results")

    n_docked = int(merged["affinity_kcal"].notna().sum()) \
        if "affinity_kcal" in merged else 0

    # A --limit run is a smoke test: it docks a handful and nulls every other
    # row's dock columns. Writing that to the frame series makes it the LATEST
    # frame, and the next stage reads the latest. A 6-ligand smoke test did
    # exactly this and broke the MM-GBSA launch three minutes later, which
    # failed loudly on `None_docked.sdf`. Partial runs report and stop.
    if limit:
        log.warning("--limit %d: NOT writing a frame. A partial run must not "
                    "become the latest frame for the next stage to read.", limit)
        return merged, None, proto, survivors, n_docked

    out = dio.write_full_frame(
        merged, approach=approach, experiment=experiment,
        stage=f"{approach}_covalent_dock",
        params={"protocol_fingerprint": proto.fingerprint(),
                "gnina_version": proto.version,
                "rank_metric": "affinity_kcal (D0015, re-measured D0028)",
                "warhead_library": str(wl.DEFAULT_LIBRARY.name),
                "n_docked": n_docked},
        inputs={"frame": frame_path})
    return merged, out, proto, survivors, n_docked
