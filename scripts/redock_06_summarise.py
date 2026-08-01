"""
Purpose: Consolidate the Pin1 redocking benchmark into one table + one JSON:
         both arms, both box sizes, top-1 and best-of-9, with drop accounting.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: redock_cases_1.csv, redock_rmsd_1.csv, redock_boxcontrol_1.csv, poses
Output: 00_outputs/blacksmith/redock_pin1/redock_benchmark_final_<N>.csv
        00_outputs/blacksmith/redock_pin1/redock_benchmark_final_<N>.json
        00_outputs/blacksmith/redock_pin1/redock_per_case_<N>.csv

TOP-1 IS THE PROTOCOL'S ANSWER; BEST-OF-9 IS WHAT IT COULD HAVE SAID. Vina
returns nine ranked modes and the pipeline consumes mode 1. Reporting both
separates the two failure modes that a single success rate confounds:

  * low top-1 AND low best-of-9  -> the correct pose is never generated
    (a SAMPLING failure -- more search depth or a smaller box might help)
  * low top-1 BUT high best-of-9 -> the correct pose is generated and then
    ranked below a wrong one (a SCORING failure -- no amount of search fixes it)

The second is the one that bears on D0041, because it is the same scoring
function being asked a strictly easier question: not "rank this binder above
other molecules" but "rank this binder's true pose above its own wrong poses".
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
import importlib.util                                   # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "redock_04_rmsd", REPO / "scripts" / "redock_04_rmsd.py")
r4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r4)

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger("redock-summary")
# Analysis outputs live under the GOVERNED root, not in the repo
# (rules/data-storage.md). See shared/outputs.py for why, and for the
# versioned-write / resolve-latest policy the append-only tree needs.
OUT = sout.Topic("blacksmith", "redock_pin1")
OUT_DIR = OUT.dir
SUCCESS_A = 2.0


def all_mode_rmsds(pose: Path, ref, tmpl) -> list[float]:
    text = pose.read_text(errors="replace")
    vals = []
    for blk in text.split("MODEL")[1:]:
        b = "MODEL" + blk.split("ENDMDL")[0] + "ENDMDL\n"
        try:
            with tempfile.TemporaryDirectory() as td:
                p = Path(td) / "m.pdbqt"
                p.write_text(b)
                s = Path(td) / "m.sdf"
                subprocess.run([r4.OBABEL, str(p), "-O", str(s)],
                               capture_output=True, timeout=120)
                raw = Chem.MolFromMolFile(str(s), sanitize=False, removeHs=False)
            raw.UpdatePropertyCache(strict=False)
            vals.append(r4.symmetric_rmsd(ref, Chem.RemoveHs(raw, sanitize=False),
                                          tmpl))
        except Exception:  # noqa: BLE001
            vals.append(float("nan"))
    return vals


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    z, p = 1.959964, k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (round(max(0.0, c - h), 3), round(min(1.0, c + h), 3))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    cases = pd.read_csv(OUT.latest("redock_cases", ".csv"))
    live = cases[cases.status == "case"].copy()
    dock = pd.read_csv(OUT.latest("redock_docking", ".csv"))
    rmsd = pd.read_csv(OUT.latest("redock_rmsd", ".csv"))
    boxc = pd.read_csv(OUT.latest("redock_boxcontrol", ".csv"))

    df = (live.merge(dock[["case_id", "self_pose_dir", "cross_pose_dir"]],
                     on="case_id", how="left")
              .merge(rmsd[["case_id", "self_rmsd_a", "cross_rmsd_a",
                           "self_affinity", "cross_affinity", "n_rot_bonds"]],
                     on="case_id", how="left")
              .merge(boxc[["case_id", "tight_top1_rmsd_a", "tight_bestofn_rmsd_a",
                           "box_edge_a"]], on="case_id", how="left"))

    # best-of-9 for the two production-box arms
    rows = []
    for c in df.itertuples():
        rec = {"case_id": c.case_id}
        tmpl = Chem.MolFromSmiles(c.smiles)
        for arm, pdir, refcol in (("self", c.self_pose_dir, c.ref_pdb),
                                  ("cross", c.cross_pose_dir, c.ref_6vaj_pdb)):
            try:
                ref = r4.reference_mol(Path(refcol), c.smiles)
                v = all_mode_rmsds(Path(pdir) / f"{c.case_id}_out.pdbqt", ref, tmpl)
                v = [x for x in v if x == x]
                rec[f"{arm}_bestofn_rmsd_a"] = min(v) if v else None
                rec[f"{arm}_n_modes"] = len(v)
            except Exception:  # noqa: BLE001
                rec[f"{arm}_bestofn_rmsd_a"] = None
        rows.append(rec)
    df = df.merge(pd.DataFrame(rows), on="case_id", how="left")
    df.to_csv(OUT.write("redock_per_case", ".csv"), index=False)

    arms = [
        ("A_self_dock_production_box", "self_rmsd_a", "self_bestofn_rmsd_a",
         "own receptor, own site, 26 A production box"),
        ("A_self_dock_tight_box", "tight_top1_rmsd_a", "tight_bestofn_rmsd_a",
         "own receptor, own site, ligand-sized box (control)"),
        ("B_cross_dock_6VAJ", "cross_rmsd_a", "cross_bestofn_rmsd_a",
         "6VAJ + box_expanded -- the exact T_1/T_2 production protocol"),
    ]
    out = []
    for name, top_col, best_col, desc in arms:
        for tier in ("all", "drug_like", "fragment"):
            sub = df if tier == "all" else df[df.tier == tier]
            t, b = sub[top_col].dropna(), sub[best_col].dropna()
            if t.empty:
                continue
            k = int((t <= SUCCESS_A).sum())
            lo, hi = wilson(k, len(t))
            out.append({
                "arm": name, "description": desc, "tier": tier, "n": len(t),
                "top1_success_2A": round(float((t <= SUCCESS_A).mean()), 4),
                "top1_n_success": k, "top1_ci95_low": lo, "top1_ci95_high": hi,
                "top1_median_rmsd_a": round(float(t.median()), 2),
                "top1_q25": round(float(t.quantile(.25)), 2),
                "top1_q75": round(float(t.quantile(.75)), 2),
                "bestof9_success_2A": round(float((b <= SUCCESS_A).mean()), 4),
                "bestof9_median_rmsd_a": round(float(b.median()), 2),
            })
    summary = pd.DataFrame(out)
    summary.to_csv(OUT.write("redock_benchmark_final", ".csv"), index=False)

    drops = cases[cases.status != "case"]
    payload = {
        "target": "Pin1 (UniProt Q13526)",
        "protocol": {
            "engine": "Vina-GPU 2.1", "search_depth": 20,
            "ligand_prep": "RDKit embed + MMFF, obabel -p 7.4 (LIGAND_PH)",
            "receptor_prep": "strip -> reduce -BUILD -> obabel -xr",
            "production_box": "box_expanded.json, 26 A",
            "rmsd": "symmetry-corrected, in place (no superposition); "
                    "validated against rdMolAlign.CalcRMS, max |diff| 0.0000 A "
                    "on 79/82 cases",
            "success_criterion_a": SUCCESS_A,
        },
        "accounting": {
            "pdb_entries_for_uniprot": 190,
            "nonpolymer_components_seen": 382,
            "components_kept_as_ligands": int((cases.classification == "ligand").sum()),
            "cases_benchmarked": int(len(live)),
            "excluded_covalent": int((cases.status == "excluded_covalent").sum()),
            "dropped": int((cases.status == "drop").sum()),
            "drop_reasons": drops[drops.status == "drop"].drop_reason
                             .value_counts().to_dict(),
            "rmsd_uncomputable": int(df.self_rmsd_a.isna().sum()),
        },
        "results": summary.to_dict(orient="records"),
    }
    (OUT.write("redock_benchmark_final", ".json")).write_text(json.dumps(payload, indent=2))

    pd.set_option("display.width", 250)
    log.info("\n%s", summary[["arm", "tier", "n", "top1_success_2A",
                              "top1_ci95_low", "top1_ci95_high",
                              "top1_median_rmsd_a", "bestof9_success_2A",
                              "bestof9_median_rmsd_a"]].to_string(index=False))


if __name__ == "__main__":
    main()
