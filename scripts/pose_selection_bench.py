"""
Purpose: does any cheap pose-selection rule beat picking at random? Measured on crystal truth.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: the redock benchmark's crystal references + docked pose ensembles
Output: 00_outputs/blacksmith/pose_selection/pose_selection_<N>.csv + a report

STAGE 2 OF `docs/ranking_rationale.md`, and the experiment the rest of the
pipeline waits on.

THE BAR IS RANDOM, NOT THE SCORE. Measured on 3IKD over 82 crystal cases: a
correct pose (<=2 A) is FINDABLE in the nine 41.5% of the time, Vina's score
picks one 18.3% of the time, and picking uniformly at random gets 19.8%. The
score is indistinguishable from a coin flip, so "beat the docking score" is not
a meaningful target -- a rule must beat **random selection** to be worth
anything, and that is a slightly HARDER bar than the status quo.

WHY THIS RUNS BEFORE ANYTHING EXPENSIVE. If a free rule clears 19.8%, it is an
immediate improvement and it calibrates what BPMD must beat. If nothing free
clears it, that is the strongest available argument that poses inside one Vina
ensemble are not separable by cheap geometry -- which would make BPMD the only
remaining option and justify a week of GPU rather than assuming it.

WHAT A "RULE" IS HERE. A function from an ensemble of poses to ONE INDEX. Every
rule returns a pose that was actually generated; none may synthesise a
conformation (`shared/pose_vector.representative` exists for exactly this
reason). A rule may not look at the crystal reference -- that is the answer.

THE RANDOM BASELINE IS COMPUTED EXACTLY, NOT SAMPLED. For each case it is
(number of poses within 2 A) / (number of poses), averaged over cases. Sampling
it would add noise to the one number every other rule is compared against.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "integration" / "app"))

from shared import outputs as sout                  # noqa: E402
from shared import pose_vector as pv                # noqa: E402
import pose3d as p3d                                # noqa: E402

log = logging.getLogger("pose-selection")

RB = Path("/data/lab_vm/append_only/inhibition/05_redock_benchmark")
OUT = sout.Topic("blacksmith", "pose_selection")
SUCCESS_A = 2.0
CLUSTER_T = 2.0          # contact-profile distance; see cluster_poses.py


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------

def heavy_xyz(path: Path) -> np.ndarray:
    out = []
    for ln in path.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        out.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.array(out)


def pose_models(path: Path) -> list[np.ndarray]:
    ms, cur = [], []
    for ln in path.read_text(errors="replace").splitlines():
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
    """Optimal-assignment heavy-atom RMSD.

    Symmetry-tolerant, and NOT `redock_04`'s graph-matched metric -- it can be
    optimistic where a graph match would refuse a mapping. Every rule and the
    baseline are scored with the SAME function, so the comparison between them
    is sound even though the absolute rates are not comparable to D0046.
    """
    from scipy.optimize import linear_sum_assignment
    if a.shape != b.shape:
        n = min(len(a), len(b))
        a, b = a[:n], b[:n]
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    r, c = linear_sum_assignment(d)
    return float(np.sqrt((d[r, c] ** 2).mean()))


# --------------------------------------------------------------------------
# the receptor basis, shared by every contact-profile rule
# --------------------------------------------------------------------------

def pocket_basis(receptor: Path) -> tuple[tuple[int, ...], dict[int, np.ndarray]]:
    resi = tuple(p3d.pocket_resi())
    want, out = set(resi), {}
    for ln in receptor.read_text(errors="replace").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        try:
            r = int(ln[22:26])
        except ValueError:
            continue
        if r in want:
            out.setdefault(r, []).append(
                (float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return resi, {k: np.array(v) for k, v in out.items()}


# --------------------------------------------------------------------------
# the rules — each maps an ensemble to ONE index
# --------------------------------------------------------------------------

def rule_score(models, vectors) -> int:
    """What the pipeline does today: Vina's best-scoring mode is written first."""
    return 0


def rule_medoid(models, vectors) -> int:
    """Contact-profile medoid — the pose most typical of the ensemble."""
    return pv.representative(vectors)


def rule_largest_cluster_medoid(models, vectors) -> int:
    """Medoid of the most populated binding mode.

    Different from `rule_medoid`: the plain medoid can sit between two clusters,
    which is the mean's failure one level down. This one commits to a mode
    first, then takes a real member of it.
    """
    labels = pv.cluster(vectors, CLUSTER_T)
    sizes = pd.Series(labels).value_counts()
    biggest = int(sizes.idxmax())
    members = [i for i, l in enumerate(labels) if l == biggest]
    return pv.representative(vectors, members)


def rule_centroid_closest(models, vectors) -> int:
    """The pose whose HEAVY-ATOM CENTROID is nearest the ensemble's centre.

    A crude spatial-consensus rule, included as a control on the contact-profile
    ones: if it does just as well, the contact vector is adding nothing over
    "where is the middle of the cloud".
    """
    cents = np.array([m.mean(axis=0) for m in models])
    mid = cents.mean(axis=0)
    return int(np.argmin(np.linalg.norm(cents - mid, axis=1)))


RULES = {
    "vina_score_top1": rule_score,
    "contact_medoid": rule_medoid,
    "largest_cluster_medoid": rule_largest_cluster_medoid,
    "centroid_closest": rule_centroid_closest,
}


# --------------------------------------------------------------------------

def evaluate(refs_dir: Path, refs_suffix: str, poses_dir: Path,
             receptor: Path, label: str) -> pd.DataFrame:
    resi, rec = pocket_basis(receptor)
    log.info("%s: pocket basis %d residues", label, len(rec))

    rows = []
    for ref in sorted(refs_dir.glob(f"*{refs_suffix}")):
        case = ref.name.replace(refs_suffix, "")
        pose = poses_dir / f"{case}_out.pdbqt"
        if not pose.is_file():
            continue
        cx = heavy_xyz(ref)
        models = pose_models(pose)
        if len(cx) == 0 or len(models) < 2:
            continue

        rmsds = [rmsd(cx, m) for m in models]
        vectors = [pv.contact_vector(m, rec, resi) for m in models]

        row = {"case": case, "arm": label, "n_modes": len(models),
               "n_good": sum(1 for r in rmsds if r <= SUCCESS_A),
               "best_rmsd": min(rmsds)}
        # EXACT random baseline, not sampled.
        row["p_random"] = row["n_good"] / len(models)
        for name, fn in RULES.items():
            try:
                idx = fn(models, vectors)
            except Exception as exc:            # noqa: BLE001
                log.warning("%s on %s: %s", name, case, str(exc)[:60])
                continue
            row[f"{name}_rmsd"] = rmsds[idx]
            row[f"{name}_hit"] = bool(rmsds[idx] <= SUCCESS_A)
            row[f"{name}_idx"] = idx
        rows.append(row)
    return pd.DataFrame(rows)


def significance(df: pd.DataFrame, name: str) -> tuple[int, float, float, float]:
    """(hits, rate, z, p) against the EXACT random expectation.

    Each case is a Bernoulli with its OWN success probability p_i = (good poses)
    / (poses), so the null is a Poisson-binomial rather than a single binomial.
    Using a pooled rate would understate the variance and manufacture
    significance.

    THIS EXISTS BECAUSE A 4.6-POINT LEAD LOOKED LIKE A RESULT. `centroid_closest`
    beat random by +4.6% on 3IKD, which reads as a finding until you notice it is
    z = +1.57 over 82 cases -- fewer than four cases of difference -- and that
    the same rule is NEGATIVE on 6VAJ. No rule may be reported as beating random
    without this number attached.
    """
    from math import erfc, sqrt
    col = f"{name}_hit"
    p = df.p_random.values
    mu, sd = p.sum(), float(np.sqrt((p * (1 - p)).sum()))
    k = int(df[col].sum())
    z = (k - mu) / sd if sd else 0.0
    return k, k / len(df), z, erfc(abs(z) / sqrt(2))


def report(df: pd.DataFrame, label: str) -> None:
    n = len(df)
    findable = (df.n_good > 0).mean()
    rand = df.p_random.mean()
    print(f"\n=== {label}  ({n} cases) ===")
    print(f"  ceiling  best-of-N (a correct pose EXISTS) : {findable:.1%}")
    print(f"  floor    random pick among the N           : {rand:.1%}")
    print()
    res = []
    for name in RULES:
        if f"{name}_hit" not in df:
            continue
        k, rate, z, pval = significance(df, name)
        res.append((name, rate, rate - rand, z, pval))
    for name, rate, delta, z, pval in sorted(res, key=lambda r: -r[1]):
        verdict = ("BEATS random" if (pval < 0.05 and delta > 0)
                   else "not distinguishable from random")
        print(f"  {name:<24} {rate:>6.1%}  {delta:+.1%}  z={z:+.2f} "
              f"p={pval:.2f}   {verdict}")
    print()
    print("  A rule must beat the FLOOR, and beat it by more than noise. "
          "Beating\n  `vina_score_top1` is not sufficient — the score is itself "
          "at chance.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--arm", choices=["3ikd", "6vaj", "both"], default="both")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    from shared import io as dio
    rec3 = sout.Topic("blacksmith", "receptor_3ikd")
    arms = {
        "3ikd": (RB / "cases_1" / "refs_3ikd", "_ref.pdb",
                 RB / "dock_1" / "cross_3ikd" / "poses",
                 dio.latest(rec3.dir, "3IKD_prepared", ".pdbqt")),
        "6vaj": (RB / "cases_1" / "refs_6vaj", "_ref6vaj.pdb",
                 RB / "dock_1" / "cross" / "poses",
                 Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdb")),
    }
    want = ["3ikd", "6vaj"] if args.arm == "both" else [args.arm]

    frames = []
    for a in want:
        refs, suf, poses, receptor = arms[a]
        if not refs.is_dir() or not poses.is_dir():
            log.warning("%s: missing %s or %s", a, refs, poses)
            continue
        df = evaluate(refs, suf, poses, receptor, a.upper())
        if df.empty:
            continue
        frames.append(df)
        report(df, a.upper())

    if not frames:
        raise SystemExit("nothing evaluated")
    out = pd.concat(frames, ignore_index=True)
    dest = OUT.write("pose_selection", ".csv")
    out.to_csv(dest, index=False)
    print(f"\nwritten -> {dest}")


if __name__ == "__main__":
    main()
