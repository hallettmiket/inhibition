"""
Purpose: Dock actives + decoys through the real protocols and grade the result (M3).
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: frozen actives; decoys_{stratum}_*.csv
Output: scored frames, gate results, and the enrichment_gate.token

RUNS THE EXACT DOWNSTREAM PROTOCOLS, not a proxy. A gate that validates a
different protocol than the one used downstream validates nothing — so the
covalent stratum goes through shared.covalent_protocol (the same pinned gnina
setup T_3 and T_4 import) and the non-covalent stratum through the same Vina
call T_1 and T_2 will use.

RESOURCE DISCIPLINE. Shared machine: 8 A100s and ~58 logged-in users. This takes
up to 6 GPUs while they are idle, drops to 4 the moment anyone else appears,
caps CPU workers at 20 of 224, and nices everything to 19.

Resumable: each ligand's result is appended to a JSONL as it completes, so a
killed run picks up where it stopped rather than restarting hours of docking.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
import sys

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import covalent_adduct as cad           # noqa: E402
from shared import covalent_protocol as cp          # noqa: E402
from shared import enrichment_gate as eg            # noqa: E402
from shared import reference_set as rs              # noqa: E402
from shared.manifest import Manifest                # noqa: E402

log = logging.getLogger("gate")

OUT_ROOT = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate/m3_enrichment")
DECOY_DIR = Path("/data/lab_vm/immutable/inhibition/decoys")
RECEPTOR_PDBQT = Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdbqt")
BOX_COVALENT = Path("/data/lab_vm/immutable/inhibition/receptor/box.json")
BOX_EXPANDED = Path("/data/lab_vm/immutable/inhibition/receptor/box_expanded.json")

# Shared machine: 8 A100s, ~58 logged-in users. Take 6 GPUs only while they are
# genuinely idle and fall back to 4 the moment anyone else appears, so a
# colleague starting a job is never squeezed by this run. Re-evaluated between
# chunks rather than fixed at launch — an 8-hour docking run outlives whatever
# the machine looked like when it started.
MAX_GPUS_IDLE = 6
MAX_GPUS_SHARED = 4
GPU_BUSY_MEMORY_MIB = 1024      # above this, someone else is using the card
CPU_WORKERS = 20                # of 224; the lab rule caps CC at 20 cores
NICE = 19
CHUNK = 60                      # ligands between GPU re-evaluations


def select_gpus() -> tuple[list[int], bool]:
    """Pick GPUs to use right now, yielding to other researchers.

    Returns
    -------
    (gpu_ids, shared)
        ``shared`` is True when other jobs are present, in which case the
        allocation drops to MAX_GPUS_SHARED.
    """
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=index,memory.used", "--format=csv,noheader,nounits"],
            text=True, timeout=30)
    except Exception as exc:  # noqa: BLE001 - no nvidia-smi means no GPU run
        log.warning("cannot query GPUs (%s); falling back to GPU 0", exc)
        return [0], True

    free, busy = [], []
    for line in out.strip().splitlines():
        idx, used = (x.strip() for x in line.split(","))
        (busy if int(used) > GPU_BUSY_MEMORY_MIB else free).append(int(idx))

    # "Busy" counts only cards WE are not already using; our own workers show up
    # here too, so this is evaluated before each chunk rather than mid-flight.
    shared = bool(busy)
    cap = MAX_GPUS_SHARED if shared else MAX_GPUS_IDLE
    chosen = free[:cap]
    if not chosen:
        chosen = free[:1] or [0]
    log.info("GPU allocation: using %s (%d idle, %d busy elsewhere) — %s",
             chosen, len(free), len(busy),
             "sharing, capped at %d" % MAX_GPUS_SHARED if shared
             else "idle, up to %d" % MAX_GPUS_IDLE)
    return chosen, shared

ACTIVE_WARHEAD_CLASS = {
    "Sulfopin": "chloroacetamide",
    "BJP-06-005-3": "chloroacetamide",
    "KPT-6566": "naphthoquinone_c2",
    "Juglone": "naphthoquinone_c2",
    "Reddi-2023-4d": "sulfamate_acetamide",
    "Reddi-2023-4g": "sulfamate_acetamide",
    # SNAr aryl chloride, NOT a chloroacetamide — different mechanism entirely.
    "Tian-chloropyrimidine-covalent-6a": "snar_chloroazine",
}


def embed_3d(smiles: str, out_sdf: Path, seed: int = 42) -> bool:
    """Embed a single conformer and minimise. Deterministic seed so runs replay."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, randomSeed=seed) != 0:
        return False
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:  # noqa: BLE001 - MMFF can fail on odd valences; pose still usable
        pass
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(out_sdf))
    w.write(m)
    w.close()
    return True


def _dock_covalent(job: dict) -> dict:
    """One covalent dock. Runs in a worker process, pinned to one GPU.

    THE GATE MUST DOCK WHAT THE APPROACHES DOCK (D0022). Actives and decoys are
    converted to adduct form here through the same `shared.covalent_adduct`
    transform T_3 and T_4 use. A gate that validated the pre-reaction form while
    the approaches docked adducts would be grading a protocol nobody runs --
    which is the same failure the module docstring warns about, arriving through
    the ligand rather than the tool.

    This is not optional bookkeeping: the protocol now serves the ADDUCT
    attachment SMARTS, which cannot match a pre-reaction ligand, so the old path
    would fail outright rather than quietly mis-score.
    """
    os.environ["CUDA_VISIBLE_DEVICES"] = str(job["gpu"])
    lig = Path(job["workdir"]) / f"{job['id']}.sdf"
    try:
        adduct = cad.to_adduct_form(job["smiles"], job["warhead_class"])
    except cad.AdductError as exc:
        return {**job, "error": f"adduct: {str(exc)[:200]}"}
    if not embed_3d(adduct.adduct_smiles, lig):
        return {**job, "error": "embed_failed"}
    try:
        r = cp.dock(lig, Path(job["workdir"]) / f"{job['id']}_docked.sdf",
                    job["warhead_class"], timeout_s=900)
        return {**job, **{k: r[k] for k in
                          ("cnn_affinity", "cnn_score", "affinity_kcal",
                           "cnn_uncalibrated_for_covalent", "protocol_fingerprint")}}
    except Exception as exc:  # noqa: BLE001 - one bad ligand must not kill the run
        return {**job, "error": str(exc)[:300]}


def _dock_vina(job: dict) -> dict:
    """One non-covalent Vina dock."""
    lig = Path(job["workdir"]) / f"{job['id']}.sdf"
    if not embed_3d(job["smiles"], lig):
        return {**job, "error": "embed_failed"}
    pdbqt = lig.with_suffix(".pdbqt")
    env = dict(os.environ)
    env["PATH"] = f"/data/lab_vm/envs/dwi_cheminf/bin:{env.get('PATH','')}"
    conv = subprocess.run(["obabel", str(lig), "-O", str(pdbqt)],
                          capture_output=True, text=True, env=env)
    if conv.returncode != 0 or not pdbqt.is_file():
        return {**job, "error": "pdbqt_conversion_failed"}
    b = json.loads(BOX_EXPANDED.read_text())
    # The `vina` pip package ships PYTHON BINDINGS, not a CLI binary. An earlier
    # version shelled out to $ENV/bin/vina, which does not exist — the whole
    # non-covalent stratum died with FileNotFoundError AFTER the covalent one
    # had already succeeded, and the run still exited 0.
    try:
        from vina import Vina

        v = Vina(sf_name="vina", cpu=1, seed=42, verbosity=0)
        v.set_receptor(str(RECEPTOR_PDBQT))
        v.set_ligand_from_file(str(pdbqt))
        v.compute_vina_maps(
            center=[b["center_x"], b["center_y"], b["center_z"]],
            box_size=[b["size_x"], b["size_y"], b["size_z"]])
        v.dock(exhaustiveness=16, n_poses=9)
        v.write_poses(str(lig.with_name(lig.stem + "_docked.pdbqt")),
                      n_poses=1, overwrite=True)
        energies = v.energies(n_poses=1)
        best = float(energies[0][0]) if len(energies) else None
    except Exception as exc:  # noqa: BLE001 - one bad ligand must not kill the run
        return {**job, "error": f"vina: {str(exc)[:280]}"}
    return {**job, "vina_affinity": best}


def build_jobs(stratum: str, workdir: Path) -> list[dict]:
    """Assemble the actives + decoys job list for one stratum."""
    master = rs.load().master
    mech = "covalent_cys113" if stratum == "covalent" else "non_covalent"
    actives = master[(master.mechanism == mech) &
                     (master.canonical_smiles != "UNVERIFIED")][
        ["name", "canonical_smiles"]].reset_index(drop=True)

    # D0031: the covalent stratum uses the CLASS-MATCHED set. decoys_covalent_2
    # pooled 104 acrylamide decoys (a class with no active at all) against
    # chloroacetamide, sulfamate and naphthoquinone actives, so D0028 found the
    # gate partly measuring chemotype rather than binding.
    # The non-covalent decoys moved OUT of the immutable tree (read-only by
    # project rule) and were rebuilt to span the actives' full mass range.
    # decoys_non_covalent_1 topped out at MW 478, so Liu-2024-C3 (MW 547) had
    # zero property matches and filter_adequately_matched would have dropped it
    # -- silently costing the SIXTH independent chemotype, which is the gate's
    # floor. Written out rather than searched for, so a run records which file
    # it consumed.
    if stratum == "covalent":
        # decoys_covalent_9: rebuilt after the snar_chloroazine class test was
        # relaxed. The narrow adduct pattern described Tian 6a specifically, so
        # the class held 3 decoys and its active was untestable; it now holds
        # 1449 and Tian gets 23 property-matched decoys.
        dfile = (Path("/data/lab_vm/append_only/inhibition/00_shared_substrate")
                 / "decoys" / "decoys_covalent_9.csv")
    else:
        dfile = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate"
                     "/decoys/decoys_non_covalent_2.csv")
    if not dfile.is_file():
        raise SystemExit(f"no decoy file at {dfile}")
    decoys = pd.read_csv(dfile)
    if stratum == "covalent":
        decoys = decoys.rename(columns={"active": "matched_active",
                                        "canonical_smiles": "smiles"})
        # `warhead_class` here is a CHEMOTYPE; docking needs a library class for
        # its SMARTS, so map through the chemotype table rather than guessing.
        from shared import decoys_classmatched as _dcm
        rep = dict(zip(_dcm.load_chemotypes()["chemotype"],
                       _dcm.load_chemotypes()["representative_class"]))
        decoys["warhead_classes"] = decoys["warhead_class"].map(rep)
    actives, decoys, excluded = eg.filter_adequately_matched(actives, decoys)
    log.info("[%s] %d actives, %d decoys (excluded: %s)",
             stratum, len(actives), len(decoys), excluded)

    jobs: list[dict] = []
    for i, r in actives.iterrows():
        jobs.append({"id": f"act_{i:03d}", "name": r["name"],
                     "smiles": r["canonical_smiles"], "label": 1,
                     "warhead_class": ACTIVE_WARHEAD_CLASS.get(r["name"], "chloroacetamide"),
                     "workdir": str(workdir), "stratum": stratum})
    for i, r in decoys.iterrows():
        wc = (str(r.get("warhead_classes", "")).split("|")[0]
              if stratum == "covalent" else "")
        jobs.append({"id": f"dec_{i:04d}", "name": r["chembl_id"],
                     "smiles": r["smiles"], "label": 0,
                     "warhead_class": wc or "chloroacetamide",
                     "workdir": str(workdir), "stratum": stratum})
    return jobs


def run_stratum(stratum: str) -> pd.DataFrame:
    """Dock a whole stratum, resumably, and return the scored frame."""
    workdir = OUT_ROOT / stratum
    workdir.mkdir(parents=True, exist_ok=True)
    # A NEW results file per ligand form (D0022). `results.jsonl` holds the
    # pre-reaction run; reusing it would silently resume onto the old docks and
    # report the old gate under a new protocol fingerprint. Append-only means it
    # stays where it is.
    results_path = workdir / f"results_{cp.LIGAND_FORM}_classmatched.jsonl"

    done: set[str] = set()
    if results_path.is_file():
        for line in results_path.read_text().splitlines():
            try:
                done.add(json.loads(line)["id"])
            except Exception:  # noqa: BLE001
                continue
        log.info("[%s] resuming: %d already done", stratum, len(done))

    jobs = [j for j in build_jobs(stratum, workdir) if j["id"] not in done]
    log.info("[%s] %d ligands to dock", stratum, len(jobs))

    done_n = 0
    with open(results_path, "a", encoding="utf-8") as fh:
        # Chunked so GPU availability is re-checked as the run proceeds. A
        # docking run lasting hours outlives whatever the machine looked like
        # when it was launched, and a colleague should not have to wait for it.
        for start in range(0, len(jobs), CHUNK):
            chunk = jobs[start:start + CHUNK]
            if stratum == "covalent":
                gpus, _shared = select_gpus()
                for i, j in enumerate(chunk):
                    j["gpu"] = gpus[i % len(gpus)]
                fn, workers = _dock_covalent, len(gpus)   # one worker per GPU
            else:
                fn, workers = _dock_vina, CPU_WORKERS

            with ProcessPoolExecutor(max_workers=workers) as ex:
                futs = {ex.submit(fn, j): j for j in chunk}
                for fut in as_completed(futs):
                    res = fut.result()
                    fh.write(json.dumps(res) + "\n")
                    fh.flush()
                    done_n += 1
                    if done_n % 25 == 0:
                        log.info("[%s] %d/%d", stratum, done_n, len(jobs))

    rows = [json.loads(l) for l in results_path.read_text().splitlines() if l.strip()]
    df = pd.DataFrame(rows)
    df = df.rename(columns={"smiles": "canonical_smiles"})
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stratum", choices=["covalent", "non_covalent", "both"])
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    os.nice(NICE)

    strata = ["covalent", "non_covalent"] if args.stratum == "both" else [args.stratum]
    thresholds = eg.load_thresholds()
    results: list[eg.GateResult] = []

    for stratum in strata:
        df = run_stratum(stratum)
        errs = df[df.get("error").notna()] if "error" in df else df.iloc[0:0]
        if len(errs):
            log.warning("[%s] %d ligand(s) failed to dock", stratum, len(errs))
        df.to_csv(OUT_ROOT / stratum / "scored.csv", index=False)

        metrics = ([("cnn_affinity", True), ("affinity_kcal", False)]
                   if stratum == "covalent" else [("vina_affinity", False)])
        for metric, higher in metrics:
            if metric not in df.columns:
                log.warning("[%s] metric %s absent — skipped", stratum, metric)
                continue
            results.append(eg.evaluate(df, metric=metric, stratum=stratum,
                                       higher_is_better=higher,
                                       thresholds=thresholds))

    eg.write_token(results)
    for r in results:
        print(f"\n[{r.stratum}/{r.metric}] {r.verdict}")
        print(f"  ROC-AUC {r.roc_auc:.3f}  CI[{r.roc_auc_ci[0]:.3f},{r.roc_auc_ci[1]:.3f}]"
              f"  EF1% {r.ef_1pct:.1f}  BEDROC {r.bedroc:.3f}")
        print(f"  {r.n_actives} actives / {r.n_decoys} decoys / {r.n_chemotypes} chemotypes")
        for reason in r.reasons:
            print(f"  - {reason}")


if __name__ == "__main__":
    main()
