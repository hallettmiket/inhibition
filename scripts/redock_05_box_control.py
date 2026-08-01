"""
Purpose: Control for the redocking benchmark -- repeat SELF-docking with a
         tight, ligand-sized box to separate "the production box is too large"
         from "the engine cannot dock this target's chemistry".
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: redock_cases_<latest>.csv
Output: 00_outputs/blacksmith/redock_pin1/redock_boxcontrol_<N>.csv

WHY THIS CONTROL EXISTS. The benchmark proper runs at the PRODUCTION box --
`box_expanded.json`, 26 A on a side (D0002) -- because a benchmark under
different settings says nothing about T_1 and T_2. But a 26 A cube around a
small, shallow, solvent-exposed domain contains a great deal of protein surface
that is not the pocket, and a low recovery rate at that box size has two very
different explanations:

  (a) the search volume is too large and the scoring function prefers one of
      the decoy surface pockets it now contains -- a PROTOCOL choice we could
      change tomorrow; or
  (b) the engine cannot reproduce this target's binding mode at all -- a much
      deeper problem that no box change fixes.

Re-running the identical protocol with only the box edge changed distinguishes
them. Everything else -- engine, SEARCH_DEPTH 20, pH 7.4 ligand prep, receptor
-- is held fixed, so the difference is attributable to the one varied term.

THE TIGHT BOX IS THE CONVENTIONAL REDOCKING BOX: the ligand's own bounding box
plus PAD_A on every side, floored at MIN_EDGE_A. That is roughly what the Astex
and PoseBusters redocking sets use, so a number measured here is comparable to
the ~60-80% those benchmarks report for a well-behaved target.
"""

from __future__ import annotations

import json
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import compute                              # noqa: E402
from shared import outputs as sout           # noqa: E402
from shared import noncovalent_dock_run as ndr          # noqa: E402

sys.path.insert(0, str(REPO / "scripts"))
import importlib.util                                   # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "redock_04_rmsd", REPO / "scripts" / "redock_04_rmsd.py")
r4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r4)

import redock_03_dock as d3                             # noqa: E402

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger("redock-boxctl")

# Analysis outputs live under the GOVERNED root, not in the repo
# (rules/data-storage.md). See shared/outputs.py for why, and for the
# versioned-write / resolve-latest policy the append-only tree needs.
OUT = sout.Topic("blacksmith", "redock_pin1")
OUT_DIR = OUT.dir
WORK = Path("/data/lab_vm/append_only/inhibition/05_redock_benchmark/boxctl_1")

PAD_A = 5.0          # padding beyond the ligand's own extent, per side
MIN_EDGE_A = 12.0    # floor, so a tiny fragment still gets a searchable volume


def tight_box(ref_pdb: Path) -> dict:
    """Ligand bounding box + PAD_A per side, floored at MIN_EDGE_A."""
    m = Chem.MolFromPDBBlock(ref_pdb.read_text(), sanitize=False, removeHs=True)
    c = m.GetConformer()
    x = np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    centre = (x.max(axis=0) + x.min(axis=0)) / 2.0
    edge = np.maximum(x.max(axis=0) - x.min(axis=0) + 2 * PAD_A, MIN_EDGE_A)
    return {"center_x": float(centre[0]), "center_y": float(centre[1]),
            "center_z": float(centre[2]), "size_x": float(edge[0]),
            "size_y": float(edge[1]), "size_z": float(edge[2])}


def best_of_modes(pose: Path, ref, tmpl) -> tuple[float | None, float | None]:
    """(top-1 RMSD, best-over-all-modes RMSD) for one output PDBQT."""
    text = pose.read_text(errors="replace")
    blocks = ["MODEL" + b.split("ENDMDL")[0] + "ENDMDL\n"
              for b in text.split("MODEL")[1:]]
    top, best = None, None
    for i, blk in enumerate(blocks):
        try:
            with tempfile.TemporaryDirectory() as td:
                p = Path(td) / "m.pdbqt"
                p.write_text(blk)
                s = Path(td) / "m.sdf"
                subprocess.run([r4.OBABEL, str(p), "-O", str(s)],
                               capture_output=True, timeout=120)
                raw = Chem.MolFromMolFile(str(s), sanitize=False, removeHs=False)
            raw.UpdatePropertyCache(strict=False)
            v = r4.symmetric_rmsd(ref, Chem.RemoveHs(raw, sanitize=False), tmpl)
        except Exception:  # noqa: BLE001
            continue
        if i == 0:
            top = v
        best = v if best is None else min(best, v)
    return top, best


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpus", default="2,3,5")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    import os
    os.nice(compute.NICE)
    gpus = [int(g) for g in args.gpus.split(",")]

    cases = pd.read_csv(OUT.latest("redock_cases", ".csv"))
    cases = cases[cases.status == "case"].copy()
    ligand_dir = Path("/data/lab_vm/append_only/inhibition/05_redock_benchmark/"
                      f"dock_1/ligands_{ndr.LIGAND_PREP_TAG}")
    log.info("%d cases; tight box = ligand extent + %.0f A/side (min %.0f A); "
             "engine/depth/pH unchanged", len(cases), PAD_A, MIN_EDGE_A)

    todo = list(cases.itertuples())
    results: list[dict] = []

    def worker(chunk, gpu):
        for c in chunk:
            src = ligand_dir / f"{c.case_id}.pdbqt"
            if not src.is_file():
                continue
            lig = WORK / c.case_id / "ligand"
            lig.mkdir(parents=True, exist_ok=True)
            dst = lig / src.name
            if not dst.is_file():
                dst.write_bytes(src.read_bytes())
            out = WORK / c.case_id / "poses"
            box = tight_box(Path(c.ref_pdb))
            try:
                d3.run_vina_gpu_on(Path(c.receptor_pdbqt), box, lig, out, gpu)
                results.append({"case_id": c.case_id, "ok": True,
                                "box_edge_a": round(float(np.mean(
                                    [box["size_x"], box["size_y"], box["size_z"]])), 1),
                                "pose_dir": str(out)})
            except Exception as exc:  # noqa: BLE001
                results.append({"case_id": c.case_id, "ok": False,
                                "error": str(exc)[:160]})
            if len(results) % 20 == 0:
                log.info("tight-box self-docking: %d/%d", len(results), len(todo))

    from concurrent.futures import ThreadPoolExecutor
    chunks = [todo[i::len(gpus)] for i in range(len(gpus))]
    with ThreadPoolExecutor(max_workers=len(gpus)) as ex:
        list(ex.map(worker, chunks, gpus))

    res = pd.DataFrame(results)
    df = cases.merge(res, on="case_id", how="left")
    rows = []
    for c in df.itertuples():
        rec = {"case_id": c.case_id, "tier": c.tier, "heavy_atoms": c.heavy_atoms,
               "box_edge_a": getattr(c, "box_edge_a", None)}
        try:
            tmpl = Chem.MolFromSmiles(c.smiles)
            ref = r4.reference_mol(Path(c.ref_pdb), c.smiles)
            top, best = best_of_modes(
                Path(c.pose_dir) / f"{c.case_id}_out.pdbqt", ref, tmpl)
            rec.update({"tight_top1_rmsd_a": top, "tight_bestofn_rmsd_a": best})
        except Exception as exc:  # noqa: BLE001
            rec["tight_error"] = str(exc)[:160]
        rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(OUT.write("redock_boxcontrol", ".csv"), index=False)

    for tier in ("all", "drug_like", "fragment"):
        a = out if tier == "all" else out[out.tier == tier]
        for col, lbl in (("tight_top1_rmsd_a", "top-1"),
                         ("tight_bestofn_rmsd_a", "best-of-9")):
            v = a[col].dropna()
            if len(v):
                log.info("[tight box | %-9s | %-9s] n=%2d  <=2A %5.1f%%  median %.2f A",
                         tier, lbl, len(v), 100 * (v <= 2.0).mean(), v.median())


if __name__ == "__main__":
    main()
