"""
Purpose: T_4 step 6 — covalent docking of survivors through the shared protocol.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the latest D4 frame (post reactivity triage)
Output: D4 with dock columns; poses under append_only/

PARITY IS THE POINT (control S3). This imports `shared.covalent_protocol`
UNCHANGED — the same pinned gnina binary, version, reactive-atom protocol and
docking parameters that T_3 uses. The protocol fingerprint is recorded per
candidate so the integration phase can verify T_3 and T_4 really did use the
same setup before offering a within-covalent comparison. "We both ran gnina" is
not that claim.

THE ADDUCT IS WHAT GETS DOCKED (D0022). gnina replaces an implicit hydrogen on
the matched atom; it does not remove a leaving group. So each candidate is
converted to its post-reaction form first — chlorine gone, sulfamate's twelve
atoms gone — and the attachment SMARTS used is the one valid on that product.
Docking the pre-reaction molecule scored a species that does not exist — the
reactive carbon bonded to Cys113 while still holding its leaving group — and,
where that group sat on a rigid ring, drove the halogen into Cys113's sulfur at
distances down to 0.89 A.

ONE DOCK PER DISTINCT ADDUCT. chloroacetamide, sulfamate_acetamide and
sulfonate_acetamide are SN2 at the same CH2 and differ only in what leaves, so
all three give an IDENTICAL adduct — verified for all 198 R-groups. Docking each
separately would burn three times the GPU time to produce three copies of one
number, and any spread between them would be pure search noise masquerading as a
warhead difference. They are docked once and the result is mapped back to every
class that shares it. The three classes' scores are then identical BY
CONSTRUCTION, which is the honest representation: what distinguishes those
warheads is kinetics, and that is the reactivity window's business, not docking's.

RANK METRIC (D0015). gnina's Vina-style `affinity` (kcal/mol, LOWER better) is
the rank metric, NOT `CNNaffinity`. The enrichment gate measured both on 6
actives and 294 warhead-bearing decoys: affinity gave ROC-AUC 0.815 with a CI
excluding 0.5 and EF1% 16.7, while CNNaffinity's CI included 0.5 and its EF1%
was 0.0 — not one known active in the top 1%. gnina itself warns that CNN
scoring is uncalibrated for covalent docking. CNNaffinity is carried as an
advisory annotation only.

ONLY SURVIVORS ARE DOCKED. Candidates stamped `rejected_at` by an earlier tier
skip this stage — the gates throttle compute, not membership. Their rows remain
in the frame with empty dock columns, so a later reweighting can still see them.

Resumable: results append to a JSONL as they complete, so a killed run picks up
rather than repeating hours of GPU time.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import covalent_adduct as cad         # noqa: E402
from shared import covalent_protocol as cp        # noqa: E402
from shared import io as dio                      # noqa: E402
from shared import warhead_library as wl          # noqa: E402

log = logging.getLogger("t4-dock")

EXPERIMENT = "04_t4_combinatorial"
OUT_ROOT = Path("/data/lab_vm/append_only/inhibition") / EXPERIMENT / "docking"

MAX_GPUS_IDLE = 6
MAX_GPUS_SHARED = 4
GPU_BUSY_MEMORY_MIB = 1024
NICE = 19
CHUNK = 60          # ligands between GPU re-evaluations


def select_gpus() -> list[int]:
    """Idle GPUs, yielding to other researchers — re-checked between chunks."""
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
    log.info("GPUs: using %s (%d idle, %d busy elsewhere)", chosen, len(free), len(busy))
    return chosen


def _dock_one(job: dict) -> dict:
    """Embed and covalently dock one candidate on its assigned GPU."""
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
        return {**job, **{k: r[k] for k in
                          ("cnn_affinity", "cnn_score", "affinity_kcal",
                           "cnn_uncalibrated_for_covalent",
                           "protocol_fingerprint", "pose_path")}}
    except Exception as exc:  # noqa: BLE001 - one bad ligand must not kill the run
        return {**job, "error": str(exc)[:300]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None,
                    help="dock only the first N survivors (smoke testing)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    os.nice(NICE)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    d = Path("/data/lab_vm/append_only/inhibition") / EXPERIMENT
    frame_path = dio.latest(d, "D4", ".parquet")
    if frame_path is None:
        raise SystemExit("no D4 frame — run 01_enumerate.py and 02_reactivity_triage.py")
    df = pd.read_parquet(frame_path)

    # Pin the protocol ONCE, up front. If it cannot be built, nothing should dock.
    proto = cp.load()
    log.info("covalent protocol fingerprint %s", proto.fingerprint()[:16])

    survivors = df[df["rejected_at"].isna()].copy()
    if args.limit:
        survivors = survivors.head(args.limit)
    log.info("%d survivors (of %d in the frame)", len(survivors), len(df))

    # --- adduct forms (D0022) -------------------------------------------------
    lib = wl.load()
    adduct_rows, adduct_failures = [], []
    for _, r in survivors.iterrows():
        try:
            a = cad.to_adduct_form(r["canonical_smiles"], r["warhead_class"],
                                   library=lib)
        except cad.AdductError as exc:
            adduct_failures.append((r["candidate_id"], str(exc)[:160]))
            continue
        adduct_rows.append({"candidate_id": r["candidate_id"],
                            "warhead_class": r["warhead_class"], **a.as_dict()})
    if adduct_failures:
        log.warning("%d candidate(s) have no well-defined adduct form and will "
                    "not be docked; first: %s", len(adduct_failures),
                    adduct_failures[0])
    adducts = pd.DataFrame(adduct_rows)
    if adducts.empty:
        raise SystemExit("no candidate produced an adduct form")

    # One dock per (adduct, class). Class is part of the key because it selects
    # the attachment SMARTS — but the three SN2 acetamides share BOTH the adduct
    # and the SMARTS, so they collapse to a single dock and the result is mapped
    # back to all three.
    adducts["dock_id"] = adducts.apply(
        lambda r: "d_" + hashlib.sha256(
            f"{r['adduct_smiles']}|{r['adduct_attachment_smarts']}"
            .encode("utf-8")).hexdigest()[:12], axis=1)
    unique = adducts.drop_duplicates("dock_id")
    log.info("%d survivors -> %d distinct adduct docks (%.1fx saving)",
             len(adducts), len(unique), len(adducts) / max(len(unique), 1))
    shared = (adducts.groupby("dock_id")["warhead_class"].nunique() > 1).sum()
    if shared:
        log.info("%d dock(s) are shared by more than one warhead class — their "
                 "adducts are the same molecule", shared)

    # A NEW results file, not a rewrite of the old one. `results.jsonl` holds the
    # pre-D0022 run keyed on candidate_id; append-only means it stays where it is
    # and stays readable, and nothing here has to delete under the data root.
    results_path = OUT_ROOT / "results_adduct.jsonl"
    done: set[str] = set()
    if results_path.is_file():
        for line in results_path.read_text().splitlines():
            try:
                done.add(json.loads(line)["dock_id"])
            except Exception:  # noqa: BLE001 - pre-D0022 rows have no dock_id
                continue
        log.info("resuming: %d already docked", len(done))

    jobs = [{"dock_id": r["dock_id"],
             "adduct_smiles": r["adduct_smiles"],
             "warhead_class": r["warhead_class"],
             "workdir": str(OUT_ROOT)}
            for _, r in unique.iterrows() if r["dock_id"] not in done]
    log.info("%d ligands to dock", len(jobs))

    n = 0
    with open(results_path, "a", encoding="utf-8") as fh:
        for start in range(0, len(jobs), CHUNK):
            chunk = jobs[start:start + CHUNK]
            gpus = select_gpus()
            for i, j in enumerate(chunk):
                j["gpu"] = gpus[i % len(gpus)]
            with ProcessPoolExecutor(max_workers=len(gpus)) as ex:
                for fut in as_completed({ex.submit(_dock_one, j): j for j in chunk}):
                    fh.write(json.dumps(fut.result()) + "\n")
                    fh.flush()
                    n += 1
                    if n % 50 == 0:
                        log.info("%d/%d docked", n, len(jobs))

    rows = [json.loads(l) for l in results_path.read_text().splitlines() if l.strip()]
    docked = pd.DataFrame(rows)
    errs = docked[docked.get("error").notna()] if "error" in docked else docked.iloc[0:0]
    if len(errs):
        log.warning("%d ligand(s) failed to dock; their rows keep empty dock columns",
                    len(errs))

    dock_cols = ("cnn_affinity", "cnn_score", "affinity_kcal",
                 "cnn_uncalibrated_for_covalent", "protocol_fingerprint",
                 "pose_path")
    # Map each dock back onto every candidate that shares its adduct.
    keep_dock = ["dock_id"] + [c for c in dock_cols if c in docked.columns]
    adduct_cols = ["candidate_id", "dock_id", "adduct_smiles",
                   "adduct_attachment_smarts", "leaving_group_smiles",
                   "adduct_atoms_removed", "adduct_approximation",
                   "adduct_degenerate_attachment"]
    docked = (adducts[[c for c in adduct_cols if c in adducts.columns]]
              .merge(docked[keep_dock].drop_duplicates("dock_id"),
                     on="dock_id", how="left"))
    keep = ["candidate_id"] + [c for c in docked.columns if c != "candidate_id"]
    # Drop any dock columns already on the frame BEFORE merging. Re-running this
    # stage on a frame that had been merged once produced affinity_kcal_x and
    # affinity_kcal_y, so `merged["affinity_kcal"]` did not exist and the
    # success counter reported 0/1683 on a run that had actually worked.
    stale = [c for c in list(dock_cols) + [c for c in keep if c != "candidate_id"]
             if c in df.columns]
    if stale:
        log.info("dropping %d stale dock column(s) before merge: %s",
                 len(stale), stale)
        df = df.drop(columns=stale)
    merged = df.merge(docked[keep].drop_duplicates("candidate_id"),
                      on="candidate_id", how="left")
    if len(merged) != len(df):
        raise RuntimeError(
            f"merge changed row count {len(df)} -> {len(merged)}; duplicate "
            "candidate_id in the frame or the results")

    out = dio.write_full_frame(
        merged, approach="t4", experiment=EXPERIMENT, stage="t4_covalent_dock",
        params={"protocol_fingerprint": proto.fingerprint(),
                "gnina_version": proto.version,
                "rank_metric": "affinity_kcal (D0015)",
                "n_docked": int(docked["affinity_kcal"].notna().sum())
                if "affinity_kcal" in docked else 0},
        inputs={"d4_frame": frame_path})

    ok = merged["affinity_kcal"].notna().sum() if "affinity_kcal" in merged else 0
    print(f"\nT_4 covalent docking -> {out}")
    print(f"  docked successfully {int(ok)} / {len(survivors)}")
    print(f"  protocol            {proto.version.strip()}")
    print(f"  fingerprint         {proto.fingerprint()[:16]}")
    if ok:
        print("\n  best affinity per warhead class (kcal/mol, lower better):")
        for cls, g in merged.dropna(subset=["affinity_kcal"]).groupby("warhead_class"):
            print(f"    {cls:22s} best {g['affinity_kcal'].min():6.2f}   "
                  f"median {g['affinity_kcal'].median():6.2f}   n={len(g)}")


if __name__ == "__main__":
    main()
