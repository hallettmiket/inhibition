"""
Purpose: measure how much of the redocking ranking gap a re-scorer closes, on OUR receptor.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: the 3IKD redocking benchmark's saved poses + transformed references
Output: 00_outputs/blacksmith/rescore_benchmark/rescore_benchmark_<N>.csv + a table

THE MEASUREMENT THAT MOTIVATES THIS. On the chemist's 3IKD, over 82 Pin1 crystal
ligands re-docked into their own receptor:

    sampling FINDS a <=2 A pose      41.5%   (34/82)
    scoring  RANKS it first          18.3%   (15/82)
    ------------------------------------------------
    pure ranking failure             23.2 points

And when a good pose exists it is in the **top 10 in 100% of cases** (median rank
4). So the right pose is nearly always present and nearly always visible in a
short window -- the scoring function simply does not put it first. That is 23
points available to a re-ranker at ZERO extra sampling cost, and it is the single
best-evidenced place to spend the compute freed by taking library-scale GROMACS
off the critical path.

This script tests whether gnina's CNN closes it, measured here rather than
inherited from a benchmark on other targets. Published numbers (GNINA vs Vina:
top-1 within 2 A, 27% -> 37% cross-docking) are an average over many proteins;
Pin1's shallow, solvent-exposed proline pocket is not the average protein, and
D0046/D0059 exist because this target has repeatedly behaved worse than norms.

WHAT IS COMPARED, on identical poses -- no re-docking, so sampling is held fixed
and only the ORDER changes:

    autodock      the order AutoDock-GPU returned          (the 18.3% baseline)
    vina_affinity gnina's Vina-like affinity, rescored
    cnn_score     gnina CNN pose score
    cnn_affinity  gnina CNN affinity
    oracle        the best pose present                    (the 41.5% ceiling)

POSE ORDER IS THE THING MOST LIKELY TO BREAK THIS, so it is asserted rather than
assumed. Poses travel PDBQT -> SDF (obabel) -> gnina, and a silent reorder or
drop anywhere would pair each RMSD with the wrong score and produce a confident
wrong answer. Counts are checked at every hop AND the SDF coordinates are
verified against the PDBQT models they came from. This project's entire defect
catalogue is values taken by position instead of identity; here position IS the
join key, so it has to be earned.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                # noqa: E402
from shared import covalent_protocol as cp        # noqa: E402
from shared import receptors as R                 # noqa: E402
import redock_3ikd_benchmark as rd                # noqa: E402

log = logging.getLogger("rescore-bench")
OUT = sout.Topic("blacksmith", "rescore_benchmark")
RB = Path("/data/lab_vm/append_only/inhibition/05_redock_benchmark")
POSES = RB / "dock_1" / "cross_3ikd" / "poses"
REFS = RB / "cases_1" / "refs_3ikd"
OBABEL = Path("/data/lab_vm/envs/dwi_cheminf/bin/obabel")

GOOD = 2.0          # A, the conventional pose-recovery threshold
_COORD_TOL = 0.05   # A; obabel round-trips coordinates, it does not move them


def sdf_coords(path: Path) -> list[np.ndarray]:
    """Heavy-atom coordinates per SDF record, in file order."""
    out, cur, n, i = [], [], 0, 0
    lines = path.read_text().splitlines()
    while i < len(lines):
        if i + 3 < len(lines) and len(lines[i + 3]) >= 6 and lines[i + 3][:3].strip().isdigit():
            try:
                n = int(lines[i + 3][:3])
            except ValueError:
                i += 1
                continue
            cur = []
            for l in lines[i + 4: i + 4 + n]:
                p = l.split()
                if len(p) >= 4 and p[3] != "H":
                    cur.append([float(p[0]), float(p[1]), float(p[2])])
            out.append(np.array(cur))
            while i < len(lines) and lines[i].strip() != "$$$$":
                i += 1
        i += 1
    return out


def to_sdf(pdbqt: Path, dest: Path) -> None:
    r = subprocess.run([str(OBABEL), "-ipdbqt", str(pdbqt), "-osdf", "-O", str(dest)],
                       capture_output=True, text=True)
    if not dest.is_file() or dest.stat().st_size == 0:
        raise RuntimeError(f"obabel produced nothing for {pdbqt.name}: {r.stderr[:200]}")


def gnina_rescore(receptor: Path, ligands: Path, gpu: int) -> pd.DataFrame:
    """Score every record in `ligands`, in file order.

    `--score_only` evaluates the pose as given; it does not minimise or re-dock,
    which is the whole point -- sampling is held fixed so any change in top-1
    accuracy is attributable to the ordering alone.
    """
    env = cp.gnina_env()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    r = subprocess.run(
        [str(cp.GNINA_BIN), "--receptor", str(receptor), "--ligand", str(ligands),
         "--score_only", "--cnn_scoring", "rescore", "--seed", "42"],
        capture_output=True, text=True, env=env, timeout=3600)
    if r.returncode != 0:
        raise RuntimeError(f"gnina rc={r.returncode}: {r.stderr[-300:]}")

    rows, cur = [], {}
    for line in r.stdout.splitlines():
        m = re.match(r"^(Affinity|CNNscore|CNNaffinity|CNNvariance):\s+(-?[\d.]+)", line)
        if not m:
            continue
        key, val = m.group(1), float(m.group(2))
        if key == "Affinity" and cur:
            rows.append(cur)
            cur = {}
        cur[key] = val
    if cur:
        rows.append(cur)
    return pd.DataFrame(rows)


def one_case(case: str, receptor: Path, gpu: int, work: Path) -> pd.DataFrame | None:
    ref = REFS / f"{case}_ref.pdb"
    pose = POSES / f"{case}_out.pdbqt"
    if not (ref.is_file() and pose.is_file()):
        return None

    cx = rd.heavy_xyz_pdb(ref)
    models = rd.models_pdbqt(pose)
    if len(cx) == 0 or not models:
        return None

    sdf = work / f"{case}.sdf"
    to_sdf(pose, sdf)
    coords = sdf_coords(sdf)

    # --- the join key is position; earn it -----------------------------------
    if len(coords) != len(models):
        raise RuntimeError(f"{case}: {len(models)} PDBQT models -> {len(coords)} SDF "
                           f"records. Conversion dropped or added poses; the "
                           f"RMSD/score pairing would be wrong.")
    for i, (a, b) in enumerate(zip(models, coords)):
        if a.shape != b.shape:
            raise RuntimeError(f"{case} pose {i}: {len(a)} heavy atoms in PDBQT vs "
                               f"{len(b)} in SDF")
        d = float(np.abs(a - b).max())
        if d > _COORD_TOL:
            raise RuntimeError(f"{case} pose {i}: coordinates moved {d:.3f} A through "
                               f"obabel. The SDF is not the pose that was scored.")

    scores = gnina_rescore(receptor, sdf, gpu)
    if len(scores) != len(models):
        raise RuntimeError(f"{case}: {len(models)} poses -> {len(scores)} gnina score "
                           f"blocks")

    rmsds = [rd.rmsd(cx, m) for m in models]
    df = scores.copy()
    df["case"] = case
    df["autodock_rank"] = np.arange(1, len(df) + 1)
    df["rmsd"] = rmsds
    return df


def summarise(all_poses: pd.DataFrame) -> pd.DataFrame:
    """Top-1 accuracy under each ordering, plus the oracle ceiling."""
    orders = {
        "autodock (baseline)":  ("autodock_rank", True),
        "gnina vina affinity":  ("Affinity", True),      # more negative is better
        "gnina CNNscore":       ("CNNscore", False),
        "gnina CNNaffinity":    ("CNNaffinity", False),
    }
    rows = []
    n = all_poses.case.nunique()
    for name, (col, asc) in orders.items():
        top = (all_poses.sort_values(col, ascending=asc)
               .groupby("case", sort=False).head(1))
        hit = (top.rmsd <= GOOD).sum()
        rows.append({"ordering": name, "top1_le_2A": hit / n * 100,
                     "n_hit": int(hit), "median_top1_rmsd": top.rmsd.median()})
    best = all_poses.groupby("case").rmsd.min()
    rows.append({"ordering": "ORACLE (best present)", "top1_le_2A": (best <= GOOD).mean() * 100,
                 "n_hit": int((best <= GOOD).sum()), "median_top1_rmsd": best.median()})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--gpu", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    os.nice(19)

    rec = sout.latest_path("blacksmith", "receptor_3ikd", "3IKD_prepared", ".pdbqt")
    R.resolve_3ikd_ian()          # refuse to benchmark against the wrong 3IKD
    log.info("receptor %s (3IKD_ian verified)", Path(rec).name)

    cases = sorted(p.name.replace("_ref.pdb", "") for p in REFS.glob("*_ref.pdb"))
    if args.limit:
        cases = cases[:args.limit]
    log.info("%d cases", len(cases))

    work = Path(tempfile.mkdtemp(prefix="rescore_"))
    frames, failed = [], []
    for i, case in enumerate(cases, 1):
        try:
            df = one_case(case, Path(rec), args.gpu, work)
            if df is None:
                failed.append((case, "missing ref or pose"))
                continue
            frames.append(df)
            if i % 10 == 0:
                log.info("[%d/%d]", i, len(cases))
        except Exception as exc:                       # noqa: BLE001
            failed.append((case, str(exc)[:160]))
            log.warning("%s: %s", case, str(exc)[:160])

    if not frames:
        log.error("nothing measured")
        return
    allp = pd.concat(frames, ignore_index=True)
    dest = OUT.write("rescore_benchmark", ".csv")
    allp.to_csv(dest, index=False)

    print(f"\n=== re-ranking {allp.case.nunique()} cases, {len(allp)} poses "
          f"(sampling held FIXED — no re-docking) ===\n")
    s = summarise(allp)
    print(f"  {'ordering':<24}{'top-1 <=2A':>12}{'n':>7}{'median RMSD':>14}")
    for r in s.itertuples():
        print(f"  {r.ordering:<24}{r.top1_le_2A:>11.1f}%{r.n_hit:>7}{r.median_top1_rmsd:>14.2f}")
    base = s[s.ordering == "autodock (baseline)"].top1_le_2A.iloc[0]
    ceil = s[s.ordering == "ORACLE (best present)"].top1_le_2A.iloc[0]
    best = s[~s.ordering.isin(["autodock (baseline)", "ORACLE (best present)"])]
    win = best.loc[best.top1_le_2A.idxmax()]
    print(f"\n  gap available to a re-ranker : {ceil - base:.1f} points")
    print(f"  best re-scorer               : {win.ordering} at {win.top1_le_2A:.1f}%")
    print(f"  gap closed                   : {(win.top1_le_2A - base) / (ceil - base) * 100:.0f}%")
    if failed:
        print(f"\n  {len(failed)} cases failed (NOT silently dropped):")
        for c, why in failed[:8]:
            print(f"    {c}: {why}")
    print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
