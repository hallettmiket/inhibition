"""
Purpose: Validate Vina-GPU against the CPU-Vina scores M3 already produced.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: M3's non-covalent ligand PDBQTs and results.jsonl
Output: per-ligand comparison, agreement statistics, and a re-run of the gate

WHY VALIDATE AT ALL. Vina-GPU is a different implementation, not a faster build
of the same binary — different search, different parallelism, different RNG
consumption. Swapping it in changes every non-covalent number downstream, so
"it is faster" is not sufficient grounds to adopt it.

WHAT WOULD COUNT AS PASSING. Two things, and the second matters more:

1. Score agreement — Pearson/Spearman against CPU Vina and the mean absolute
   difference in kcal/mol. Docking is stochastic, so exact equality is not
   expected; systematic offset or poor rank correlation is.

2. **The enrichment verdict must reproduce.** M3 found non-covalent docking
   barely enriches (AUC 0.535, D0016). If Vina-GPU reports materially better
   enrichment on the same ligands, that is not good news — it means the two
   engines disagree about which molecules bind, and the disagreement, not the
   speedup, becomes the finding. A speed win bought by changing the answer is
   not a speed win.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import enrichment_gate as eg          # noqa: E402
from shared.manifest import Manifest              # noqa: E402

log = logging.getLogger("vgpu-validate")

M3_DIR = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate/"
              "m3_enrichment/non_covalent")
OUT_ROOT = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate/"
                "vina_gpu_validation")
OUT_DIR = OUT_ROOT  # rebound per search_depth in main()
VINA_GPU = Path("/data/lab_vm/envs/dwi_vinagpu/bin/vina-gpu")
RECEPTOR = Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdbqt")
BOX = Path("/data/lab_vm/immutable/inhibition/receptor/box_expanded.json")

# Agreement thresholds. Chosen so a real regression fails and ordinary docking
# stochasticity does not.
MIN_SPEARMAN = 0.80          # rank agreement is what a ranking pipeline needs
MAX_MEAN_ABS_DIFF = 1.0      # kcal/mol
MAX_AUC_DRIFT = 0.10         # the verdict must not move materially


def stage_ligands(dest: Path) -> dict[str, dict]:
    """Copy M3's input PDBQTs into a flat directory for screening mode."""
    dest.mkdir(parents=True, exist_ok=True)
    rows = [json.loads(l) for l in (M3_DIR / "results.jsonl").read_text().splitlines()
            if l.strip()]
    meta: dict[str, dict] = {}
    for r in rows:
        if r.get("vina_affinity") is None:
            continue
        src = M3_DIR / f"{r['id']}.pdbqt"
        if not src.is_file():
            continue
        (dest / src.name).write_bytes(src.read_bytes())
        meta[r["id"]] = {"label": r["label"], "name": r["name"],
                         "canonical_smiles": r["smiles"],
                         "cpu_vina": r["vina_affinity"]}
    log.info("staged %d ligands already scored by CPU Vina", len(meta))
    return meta


def run_vina_gpu(ligand_dir: Path, out_dir: Path, *, threads: int = 8000,
                 search_depth: int = 10) -> float:
    """Run Vina-GPU in virtual-screening mode. Returns wall-clock seconds."""
    out_dir.mkdir(parents=True, exist_ok=True)
    b = json.loads(BOX.read_text())
    cmd = [str(VINA_GPU),
           "--receptor", str(RECEPTOR),
           "--ligand_directory", str(ligand_dir),
           "--output_directory", str(out_dir),
           "--center_x", str(b["center_x"]), "--center_y", str(b["center_y"]),
           "--center_z", str(b["center_z"]),
           "--size_x", str(b["size_x"]), "--size_y", str(b["size_y"]),
           "--size_z", str(b["size_z"]),
           "--thread", str(threads), "--search_depth", str(search_depth)]
    log.info("running Vina-GPU: %s", " ".join(cmd[-6:]))
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
    dt = time.time() - t0
    (OUT_DIR / "vina_gpu_stdout.log").write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
    if p.returncode != 0:
        raise RuntimeError(f"Vina-GPU failed ({p.returncode}); see vina_gpu_stdout.log\n"
                           f"{(p.stderr or p.stdout)[-1500:]}")
    log.info("Vina-GPU finished in %.1f s", dt)
    return dt


_RESULT_RE = re.compile(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)")


def collect_scores(out_dir: Path) -> dict[str, float]:
    """Best score per ligand from Vina-GPU's output PDBQTs."""
    scores: dict[str, float] = {}
    for f in out_dir.glob("*.pdbqt"):
        lig_id = f.stem.replace("_out", "")
        m = _RESULT_RE.search(f.read_text(errors="replace"))
        if m:
            scores[lig_id] = float(m.group(1))
    log.info("parsed %d Vina-GPU scores", len(scores))
    return scores


def compare(meta: dict[str, dict], gpu: dict[str, float], elapsed: float) -> dict:
    """Score agreement plus a re-run of the enrichment gate on GPU scores."""
    rows = [{**v, "id": k, "gpu_vina": gpu[k]} for k, v in meta.items() if k in gpu]
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("no ligands scored by both engines")

    pearson = float(df.cpu_vina.corr(df.gpu_vina))
    spearman = float(df.cpu_vina.corr(df.gpu_vina, method="spearman"))
    mad = float((df.cpu_vina - df.gpu_vina).abs().mean())
    bias = float((df.gpu_vina - df.cpu_vina).mean())

    thresholds = eg.load_thresholds()
    cpu_res = eg.evaluate(df.rename(columns={"cpu_vina": "m"}), metric="m",
                          stratum="non_covalent_cpu", higher_is_better=False,
                          thresholds=thresholds)
    gpu_res = eg.evaluate(df.rename(columns={"gpu_vina": "m"}), metric="m",
                          stratum="non_covalent_gpu", higher_is_better=False,
                          thresholds=thresholds)
    auc_drift = abs(gpu_res.roc_auc - cpu_res.roc_auc)

    checks = {
        "spearman_ok": spearman >= MIN_SPEARMAN,
        "mean_abs_diff_ok": mad <= MAX_MEAN_ABS_DIFF,
        "auc_reproduces": auc_drift <= MAX_AUC_DRIFT,
        "verdict_matches": gpu_res.verdict == cpu_res.verdict,
    }
    verdict = "ADOPT" if all(checks.values()) else "DO_NOT_ADOPT"

    out = {
        "n_compared": len(df), "elapsed_s": elapsed,
        "pearson": pearson, "spearman": spearman,
        "mean_abs_diff_kcal": mad, "mean_bias_kcal": bias,
        "cpu": cpu_res.to_dict(), "gpu": gpu_res.to_dict(),
        "auc_drift": auc_drift, "checks": checks, "verdict": verdict,
    }
    (OUT_DIR / "comparison.json").write_text(json.dumps(out, indent=2) + "\n")
    df.to_csv(OUT_DIR / "per_ligand.csv", index=False)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8000)
    ap.add_argument("--search-depth", type=int, default=10)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    # Separate output per search_depth: the whole point of sweeping it is to
    # compare runs, which needs their artifacts kept apart.
    global OUT_DIR
    OUT_DIR = OUT_ROOT / f"sd{args.search_depth}"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    lig_dir, out_dir = OUT_DIR / "ligands", OUT_DIR / "poses"
    meta = stage_ligands(lig_dir)
    elapsed = run_vina_gpu(lig_dir, out_dir, threads=args.threads,
                           search_depth=args.search_depth)
    res = compare(meta, collect_scores(out_dir), elapsed)

    (Manifest(stage="validate_vina_gpu", approach="shared",
              params={"threads": args.threads, "search_depth": args.search_depth})
     .add_input("m3_results", M3_DIR / "results.jsonl")
     .add_output("comparison", OUT_DIR / "comparison.json")
     .note(f"verdict {res['verdict']}; spearman {res['spearman']:.3f}; "
           f"AUC drift {res['auc_drift']:.3f}")
     .write(OUT_DIR, filename="validation_manifest.json"))

    print(f"\n{'='*66}\nVina-GPU validation: {res['verdict']}\n{'='*66}")
    print(f"  ligands compared    {res['n_compared']}")
    print(f"  wall-clock          {elapsed:.1f} s  (CPU Vina took ~20 min on 20 cores)")
    print(f"  Pearson r           {res['pearson']:.3f}")
    print(f"  Spearman rho        {res['spearman']:.3f}   (need >= {MIN_SPEARMAN})")
    print(f"  mean |diff|         {res['mean_abs_diff_kcal']:.3f} kcal/mol  (need <= {MAX_MEAN_ABS_DIFF})")
    print(f"  mean bias           {res['mean_bias_kcal']:+.3f} kcal/mol")
    print(f"  ROC-AUC  cpu {res['cpu']['roc_auc']:.3f} -> gpu {res['gpu']['roc_auc']:.3f} "
          f"(drift {res['auc_drift']:.3f}, need <= {MAX_AUC_DRIFT})")
    print(f"  verdict  cpu {res['cpu']['verdict']} -> gpu {res['gpu']['verdict']}")
    for k, v in res["checks"].items():
        print(f"    {'PASS' if v else 'FAIL'}  {k}")


if __name__ == "__main__":
    main()
