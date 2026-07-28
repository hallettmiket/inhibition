"""
Purpose: T_1 step 1 — de novo generation into the Pin1 pocket with DiffSBDD.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: config/approaches/t1_de_novo.yaml, the prepared receptor, the DiffSBDD checkpoint
Output: the D1 frame — generated molecules with their 3D poses retained

WHAT T_1 IS FOR. It is the only approach with no seed. T_2 inherits ATRA's
chemotype, T_3 and T_4 inherit sulfopin's — so if all four converge on similar
chemistry, three of them were always going to. T_1 is what makes convergence
between approaches evidence rather than an artifact of shared starting material.
Its value is orthogonality, and that is worth more here than its hit rate.

POCKET-CONDITIONED, NOT LIGAND-CONDITIONED. DiffSBDD can define the pocket from
a reference ligand, which would be the easy route: point it at sulfopin in
6VAJ. That would hand T_1 the chemotype it exists to avoid, so the pocket is
defined by residue list instead. The site is the same; the prior is not.

THE POSES ARE KEPT. DiffSBDD emits 3D coordinates, not just SMILES, and those
coordinates are the model's own hypothesis about how the molecule sits. They are
written alongside the frame so docking can be compared against them later — a
generated pose that docking reproduces is a different quality of evidence from
one it does not.

SIZE BOUNDS ARE A REAL FILTER HERE. A seeded approach cannot emit a 6-atom
fragment or 90-atom sprawl; a diffusion model does both. The bounds in config
stamp rather than delete, per the usual rule.

RANKING WILL BE WEAK, AND THAT IS RECORDED UP FRONT. T_1 is non-covalent, and
the enrichment gate found non-covalent docking barely enriches on this receptor
(D0016: ROC-AUC 0.535, CI [0.215, 0.855], EF1% 0.0). T_1's eventual ranking
carries that caveat into the GUI rather than presenting a score as if it were
calibrated.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                      # noqa: E402
from shared import smiles as smi                  # noqa: E402

log = logging.getLogger("t1-generate")

EXPERIMENT = "01_t1_de_novo"
APPROACH = "t1"
CONFIG = REPO / "config" / "approaches" / "t1_de_novo.yaml"
DIFFSBDD_ENV = Path("/data/lab_vm/envs/dwi_diffsbdd")
RECEPTOR_PDB = Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdb")
WORK_ROOT = Path("/data/lab_vm/append_only/inhibition") / EXPERIMENT / "generation"
NICE = 19


def load_config() -> dict:
    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    g, f = cfg["generation"], cfg.get("filters") or {}
    return {
        "repo": Path(g["repo"]),
        "checkpoint": Path(g["checkpoint"]),
        "resi_list": list(g["resi_list"]),
        "n_samples": int(g["n_samples"]),
        "batch_size": int(g.get("batch_size", 64)),
        "sanitize": bool(g.get("sanitize", True)),
        "relax": bool(g.get("relax", False)),
        "timesteps": g.get("timesteps"),
        "min_heavy": int(f.get("min_heavy_atoms", 0)),
        "max_heavy": int(f.get("max_heavy_atoms", 10_000)),
        "size_classes": cfg.get("size_classes") or {},
        "mechanism": cfg["approach"]["mechanism"],
        "enrichment_verdict": cfg.get("docking", {}).get("enrichment_verdict"),
    }


def run_diffsbdd(cfg: dict, workdir: Path, n_samples: int,
                 timeout_s: int = 21600) -> Path:
    """Invoke DiffSBDD's own generate_ligands.py, unmodified.

    Run through the upstream entry point rather than importing its internals, so
    a checkpoint that expects a particular version of that script keeps working.
    """
    out_sdf = workdir / "generated.sdf"
    script = cfg["repo"] / "generate_ligands.py"

    # DiffSBDD asserts n_samples % batch_size == 0 and dies with a bare
    # AssertionError naming neither value. Reconcile them here and say so,
    # rounding UP so a requested count is never silently reduced.
    batch = min(cfg["batch_size"], n_samples)
    if n_samples % batch:
        adjusted = ((n_samples // batch) + 1) * batch
        log.info("n_samples %d is not a multiple of batch %d; generating %d",
                 n_samples, batch, adjusted)
        n_samples = adjusted
    if not script.is_file():
        raise SystemExit(f"DiffSBDD not staged at {script}\n"
                         "  python -m shared.sources stage --only diffsbdd_repo")
    if not cfg["checkpoint"].is_file():
        raise SystemExit(f"checkpoint not staged: {cfg['checkpoint']}")

    cmd = [str(DIFFSBDD_ENV / "bin" / "python"), str(script),
           str(cfg["checkpoint"]),
           "--pdbfile", str(RECEPTOR_PDB),
           "--outfile", str(out_sdf),
           "--n_samples", str(n_samples),
           "--batch_size", str(batch),
           "--resi_list", *cfg["resi_list"]]
    if cfg["sanitize"]:
        cmd.append("--sanitize")
    if cfg["relax"]:
        cmd.append("--relax")
    if cfg["timesteps"]:
        cmd += ["--timesteps", str(cfg["timesteps"])]

    log.info("running DiffSBDD for %d samples (the slow part)", n_samples)
    proc = subprocess.run(cmd, cwd=cfg["repo"], capture_output=True, text=True,
                          timeout=timeout_s)
    (workdir / "diffsbdd.log").write_text(proc.stdout + "\n" + proc.stderr,
                                          encoding="utf-8")
    if proc.returncode != 0:
        raise SystemExit(
            f"DiffSBDD failed ({proc.returncode}); see {workdir/'diffsbdd.log'}\n"
            + proc.stderr[-900:])
    if not out_sdf.is_file():
        raise SystemExit(f"DiffSBDD exited 0 but wrote no {out_sdf.name}")
    return out_sdf


def _size_class(hac: int, bounds: dict) -> str:
    """Label a molecule fragment / lead_like / out_of_range by heavy-atom count.

    A LABEL, not a verdict. T_1 is the only approach that spans fragment and
    lead space, because it is the only one with no seed dictating size. Pooling
    a 12-heavy-atom fragment with a T_4 lead on a physchem plot without saying
    so would make T_1 look like an outlier when it is really a different kind of
    object.
    """
    for name, (lo, hi) in bounds.items():
        if lo <= hac <= hi:
            return name
    return "out_of_range"


def parse_generated(sdf: Path, cfg: dict) -> pd.DataFrame:
    """Read the generated poses, dedup on InChIKey, stamp size bounds."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    rows: list[dict] = []
    seen: set[str] = set()
    n_read = n_unparseable = 0
    for idx, mol in enumerate(Chem.SDMolSupplier(str(sdf), sanitize=True)):
        n_read += 1
        if mol is None:
            n_unparseable += 1
            continue
        canon = smi.canonical(Chem.MolToSmiles(mol))
        if canon is None:
            n_unparseable += 1
            continue
        key = smi.inchikey(canon)
        if not key or key in seen:
            continue
        seen.add(key)
        hac = mol.GetNumHeavyAtoms()
        reject = None
        if hac < cfg["min_heavy"]:
            reject = "degenerate_too_small"
        elif hac > cfg["max_heavy"]:
            reject = "too_large"
        size_class = _size_class(hac, cfg["size_classes"])
        rows.append({"canonical_smiles": canon,
                     "candidate_id": smi.candidate_id(canon, prefix=APPROACH),
                     "approach": APPROACH,
                     "mechanism": cfg["mechanism"],
                     "generated_pose_index": idx,
                     "generated_heavy_atoms": hac,
                     "size_class": size_class,
                     "rejected_at": reject})
    log.info("read %d generated molecules, %d unparseable, %d unique kept",
             n_read, n_unparseable, len(rows))
    if n_unparseable:
        log.warning("%d/%d generated structures did not sanitize — normal for a "
                    "diffusion model, recorded rather than hidden",
                    n_unparseable, n_read)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="T_1: de novo generation into Pin1.")
    ap.add_argument("--n", type=int, default=None, help="override n_samples")
    ap.add_argument("--reuse", action="store_true",
                    help="parse an existing generated.sdf instead of re-running")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    WORK_ROOT.mkdir(parents=True, exist_ok=True)

    cfg = load_config()
    n = args.n or cfg["n_samples"]
    out_sdf = WORK_ROOT / "generated.sdf"
    if not (args.reuse and out_sdf.is_file()):
        out_sdf = run_diffsbdd(cfg, WORK_ROOT, n)

    df = parse_generated(out_sdf, cfg)
    if df.empty:
        raise SystemExit("DiffSBDD produced nothing parseable")

    out = dio.write_full_frame(
        df, approach=APPROACH, experiment=EXPERIMENT, stage="t1_generate",
        params={"engine": "diffsbdd",
                "checkpoint": str(cfg["checkpoint"]),
                "pocket_mode": "resi_list (NOT ref_ligand — no chemotype prior)",
                "resi_list": cfg["resi_list"],
                "n_samples_requested": n,
                "min_heavy_atoms": cfg["min_heavy"],
                "max_heavy_atoms": cfg["max_heavy"],
                "enrichment_verdict": cfg["enrichment_verdict"],
                "poses_retained": str(out_sdf)},
        inputs={"checkpoint": cfg["checkpoint"], "receptor": RECEPTOR_PDB})

    n_ok = int(df["rejected_at"].isna().sum())
    print(f"\nT_1 generation -> {out}")
    print(f"  {len(df)} unique molecules, {n_ok} inside the size bounds "
          f"[{cfg['min_heavy']}, {cfg['max_heavy']}] heavy atoms")
    print(f"  poses retained at {out_sdf}")
    if "generated_heavy_atoms" in df:
        s = df["generated_heavy_atoms"]
        print(f"  heavy atoms: median {s.median():.0f}, range {s.min()}-{s.max()}")
    for reason, k in df["rejected_at"].value_counts().items():
        print(f"    stamped {reason}: {k}")
    if "size_class" in df:
        print("  size classes (label, not verdict):")
        for cls, k in df["size_class"].value_counts().items():
            print(f"    {cls:14s} {k}")
    print(f"\n  NOTE: non-covalent docking is {cfg['enrichment_verdict']} on this "
          "receptor (D0016); T_1's ranking carries that caveat.")


if __name__ == "__main__":
    main()
