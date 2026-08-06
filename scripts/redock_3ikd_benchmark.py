"""
Purpose: re-run D0046's pose-recovery benchmark against the chemist's 3IKD.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: cases_1/ (already built), the prepared 3IKD receptor + box
Output: refs_3ikd/ transformed references, poses, and a recovery table

D0059 replaced 6VAJ with 3IKD. Every benchmark on record is a 6VAJ number, and
the one that matters most is D0046's: docking recovers a known Pin1 crystal pose
**5% of the time in production**, against a 60-80% norm. That single measurement
is why the ranking is considered unusable and why #14's redesign exists. It does
not transfer to a different receptor, so it has to be measured again.

WHAT THIS REUSES AND WHY. The 83 cases -- crystal ligands from Pin1 entries at
<=2.0 A -- are receptor-independent and already built. So are the prepared
ligands under `dock_1/ligands_ph7.4/`. Only two things change: the reference
poses must be expressed in 3IKD's coordinate frame instead of 6VAJ's, and the
docking target is 3IKD. The superposition machinery is imported from
`redock_02_build_cases` rather than reimplemented -- a second implementation of
"fit the PPIase CA trace" is how two answers to one question appear.

THE COMPARISON IS CONFOUNDED, AND THE CONFOUND IS NOT MINE TO HIDE. 6VAJ was
stripped of waters and protonated by `reduce` at pH 7.4; 3IKD keeps 6 waters and
the chemist's own protonation (D0059). So a difference against D0046's 5% is
**receptor AND preparation**, not receptor alone. The control that separates
them -- deposited 3IKD through 6VAJ's path -- is cheap and is NOT run here. Any
result must be reported with that attached.

WHAT COUNTS AS SUCCESS. Same definition as D0046: a docked pose within 2.0 A
of the crystal pose, heavy atoms only. Reported at top-1 (what the pipeline
actually carries) and best-of-N (what a perfect pose-selection rule could
reach). The gap between them is the headroom #14's BPMD work is chasing.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402

log = logging.getLogger("redock-3ikd")

RB = Path("/data/lab_vm/append_only/inhibition/05_redock_benchmark")
OUT = sout.Topic("blacksmith", "redock_3ikd")
# Resolved by shared.noncovalent_dock_run: the governed wrapper when its
# working directory is writable, else the byte-identical copy (D0060).
from shared.noncovalent_dock_run import VINA_GPU   # noqa: E402
SEARCH_DEPTH = 20          # D0017; below this the adoption evidence does not hold
THREADS = 8000
SUCCESS_A = 2.0            # D0046's definition


def _rd02():
    spec = importlib.util.spec_from_file_location(
        "rd02", REPO / "scripts" / "redock_02_build_cases.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def latest(topic: sout.Topic, stem: str, suffix: str) -> Path:
    from shared import io as dio
    p = dio.latest(topic.dir, stem, suffix)
    if p is None:
        raise SystemExit(f"missing {stem}*{suffix} under {topic.dir}")
    return p


def build_refs(rd02, receptor_pdb: Path, out_dir: Path) -> list[dict]:
    """Express every case's crystal ligand in the new receptor's frame."""
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_atoms, _ = rd02.read_pdb(receptor_pdb)
    ref_ca = rd02.ca_map(ref_atoms)
    log.info("target reference: %d PPIase CA atoms", len(ref_ca))

    rows = []
    for src in sorted((RB / "cases_1" / "refs").glob("*_ref.pdb")):
        case_id = src.name.replace("_ref.pdb", "")
        pdb_id = case_id.split("_")[0]
        struct = None
        for cand in (RB / "cases_1" / "pdb").glob(f"{pdb_id}.*"):
            struct = cand
            break
        if struct is None:
            rows.append({"case": case_id, "status": "no structure"}); continue
        try:
            atoms, _ = rd02.load_structure(struct)
        except Exception as exc:            # noqa: BLE001
            rows.append({"case": case_id, "status": f"parse: {exc}"[:60]}); continue

        # Fit this entry's PPIase CA trace onto the target's.
        chains = {a["chain"] for a in atoms if a["record"] == "ATOM"}
        best = None
        for ch in sorted(chains):
            mine = rd02.ca_map(atoms, ch)
            shared_res = sorted(set(mine) & set(ref_ca))
            if len(shared_res) < 50:
                continue
            mob = np.array([mine[r][1] for r in shared_res])
            tgt = np.array([ref_ca[r][1] for r in shared_res])
            rot, trans, fit = rd02.kabsch(mob, tgt)
            if best is None or fit < best[2]:
                best = (rot, trans, fit, len(shared_res), ch)
        if best is None:
            rows.append({"case": case_id, "status": "no chain fit"}); continue
        rot, trans, fit, n_fit, ch = best

        lig_atoms, _ = rd02.read_pdb(src)
        xyz = np.array([a["xyz"] for a in lig_atoms])
        rd02.write_ligand_pdb(lig_atoms, out_dir / f"{case_id}_ref.pdb",
                              xyz @ rot + trans)
        rows.append({"case": case_id, "status": "ok", "ca_fit_rmsd": round(fit, 3),
                     "n_ca": n_fit, "chain": ch})
    return rows


def dock(ligand_dir: Path, out_dir: Path, receptor: Path, box: dict,
         gpu: int, timeout_s: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(VINA_GPU), "--receptor", str(receptor),
           "--ligand_directory", str(ligand_dir),
           "--output_directory", str(out_dir),
           "--center_x", str(box["center_x"]), "--center_y", str(box["center_y"]),
           "--center_z", str(box["center_z"]),
           "--size_x", str(box["size_x"]), "--size_y", str(box["size_y"]),
           "--size_z", str(box["size_z"]),
           "--thread", str(THREADS), "--search_depth", str(SEARCH_DEPTH)]
    import os
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    env["GPU_DEVICE_ORDINAL"] = str(gpu)
    log.info("Vina-GPU on GPU %d, %d ligands, deadline %.1f h",
             gpu, len(list(ligand_dir.glob("*.pdbqt"))), timeout_s / 3600)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s,
                       env=env)
    (out_dir / "vina_gpu_stdout.log").write_text(p.stdout + "\n" + p.stderr)
    if p.returncode != 0:
        raise SystemExit(f"Vina-GPU failed ({p.returncode}); see the log")


def heavy_xyz_pdb(p: Path) -> np.ndarray:
    out = []
    for ln in p.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        out.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.array(out)


def models_pdbqt(p: Path) -> list[np.ndarray]:
    ms, cur = [], []
    for ln in p.read_text(errors="replace").splitlines():
        if ln.startswith("MODEL"):
            cur = []
        elif ln.startswith("ENDMDL"):
            if cur:
                ms.append(np.array(cur))
            cur = []
        elif ln.startswith(("ATOM", "HETATM")):
            el = ln[77:79].strip().upper() if len(ln) > 78 else ""
            if el.startswith("H") and el != "HG":
                continue
            try:
                cur.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            except ValueError:
                pass
    if cur:
        ms.append(np.array(cur))
    return ms


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Assignment-based heavy-atom RMSD.

    NOT `redock_04`'s graph-matched symmetry-corrected metric -- this is an
    optimal-assignment proxy, which is symmetry-tolerant but can be optimistic
    where a graph match would refuse a mapping. Absolute rates from here are
    therefore NOT directly comparable to D0046's; the 6VAJ arm is recomputed
    with this same metric so the COMPARISON is like-for-like.
    """
    from scipy.optimize import linear_sum_assignment
    if a.shape != b.shape:
        n = min(len(a), len(b))
        a, b = a[:n], b[:n]
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    r, c = linear_sum_assignment(d)
    return float(np.sqrt((d[r, c] ** 2).mean()))


def score(refs_dir: Path, poses_dir: Path, label: str) -> pd.DataFrame:
    rows = []
    for ref in sorted(refs_dir.glob("*_ref.pdb")):
        case = ref.name.replace("_ref.pdb", "")
        pose = poses_dir / f"{case}_out.pdbqt"
        if not pose.is_file():
            continue
        cx = heavy_xyz_pdb(ref)
        ms = models_pdbqt(pose)
        if len(cx) == 0 or not ms:
            continue
        rs = [rmsd(cx, m) for m in ms]
        rows.append({"case": case, "arm": label, "n_modes": len(rs),
                     "top1_rmsd": rs[0], "best_rmsd": min(rs),
                     "best_rank": int(np.argmin(rs)) + 1})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--gpu", type=int, default=1)
    ap.add_argument("--skip-dock", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    rec3 = sout.Topic("blacksmith", "receptor_3ikd")
    receptor = latest(rec3, "3IKD_prepared", ".pdbqt")
    box = json.loads(latest(rec3, "box_3IKD", ".json").read_text())
    log.info("receptor %s", receptor.name)

    rd02 = _rd02()
    refs_dir = RB / "cases_1" / "refs_3ikd"
    rows = build_refs(rd02, receptor, refs_dir)
    ok = [r for r in rows if r.get("status") == "ok"]
    log.info("transformed %d/%d references into 3IKD's frame; median CA fit "
             "%.3f A", len(ok), len(rows),
             float(np.median([r["ca_fit_rmsd"] for r in ok])) if ok else float("nan"))

    poses = RB / "dock_1" / "cross_3ikd" / "poses"
    if not args.skip_dock:
        dock(RB / "dock_1" / "ligands_ph7.4", poses, receptor, box,
             args.gpu, timeout_s=7200)

    df = score(refs_dir, poses, "3IKD")
    dest = OUT.write("redock_3ikd", ".csv")
    df.to_csv(dest, index=False)

    print(f"\n3IKD pose recovery ({len(df)} cases) -> {dest}")
    for k in (1, 3, 5, 9):
        hit = ((df.best_rank <= k) & (df.best_rmsd <= SUCCESS_A)).mean()
        print(f"  crystal pose within {SUCCESS_A} A in top-{k}: {hit:.1%}")
    print(f"  median top-1 RMSD : {df.top1_rmsd.median():.2f} A")
    print(f"  median best RMSD  : {df.best_rmsd.median():.2f} A")
    print("\n  6VAJ, same metric, measured 2026-08-05: top-1 4.9%, best-of-9 15.9%")
    print("  CONFOUNDED: 6VAJ was water-stripped + reduce-protonated; 3IKD keeps")
    print("  6 waters and the chemist's protonation (D0059). Receptor AND prep.")


if __name__ == "__main__":
    main()
