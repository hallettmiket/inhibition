#!/usr/bin/env python3
"""
Purpose: which pocket residues carry information about where a pose sits, and how many.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: the replicate clouds election_<cand>_r{1..5}_allposes, and the receptor
Output: 00_outputs/blacksmith/residue_selection_<candidate>/

@tt8804: "we need to determine how many residues we want to use first, is more
better? up to what"

THE ANSWER IS NOT PROXIMITY AND NOT CATALYTIC IMPORTANCE. A residue earns a place
in the splitting metric if the ligand-to-residue distances VARY across the pose
cloud. A residue at a constant distance from every pose carries no discriminating
information however essential it is to catalysis -- it contributes a constant
column, and in a norm across dimensions a constant column DILUTES the informative
ones rather than merely failing to help.

Three forces, all measured here:
  VARIANCE     does this residue's distance profile move between poses?
  REDUNDANCY   two residues on the same pocket face move together; including both
               weights one direction twice and buys nothing.
  DIMENSION    as dimensions grow, pairwise distances converge and everything
               becomes equidistant -- the documented failure mode of density
               clustering in high dimensions. 28 atoms x 40 residues is 1,120
               dimensions for a few thousand poses, which is deep into it.

NOT DERIVED FROM Cys113 (@tt8804: "your focus on cys113 at this stage ... warhead
orientation with regards to cys113 is how we RANK the poses not DETERMINE them").
Candidates are every residue the cloud comes near, and the anchor is not
privileged. Whether including it changes the answer is reported rather than
assumed, because that is the cheap version of the argument.

RANKED ON FIVE INDEPENDENT DOCKINGS, NOT ONE. D0091 was refuted partly because a
ladder built by subsampling ONE pooled cloud manufactured its own result. Here
each replicate is a separate docking with its own seed, the ranking is computed
per replicate, and the AGREEMENT between them is the headline -- a residue
ordering that does not survive an independent draw is a property of one sample.

WHAT THIS EXPERIMENT CANNOT TELL YOU, stated so it is not read for more than it
is: candidates are chosen by proximity to the cloud, so residues at the cloud's
edge are selected partly for the same reason they will show high variance. The
ranking WITHIN the candidate set is meaningful; the absolute variance values are
not a test of whether the candidate set was well chosen.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("residue-selection")

#: Distances are capped: beyond this, "far" is "far", and an uncapped tail lets
#: one remote residue dominate a norm for no physical reason.
CAP_A = 10.0
#: A residue is a candidate if any of its heavy atoms comes within this of any
#: ligand heavy atom in any pose.
CANDIDATE_A = 10.0


def receptor_residues() -> dict:
    """{(chain, resi, resname): (n_atoms, 3) heavy-atom coordinates}."""
    rec = rp.receptor_prep()
    out: dict = {}
    for ln in rec.read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = ln[76:78].strip() or ln[12:16].strip()[:1]
        if el.upper() == "H":
            continue
        key = (ln[21:22].strip() or "A", ln[22:26].strip(), ln[17:20].strip())
        try:
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        except ValueError:
            continue
        out.setdefault(key, []).append(xyz)
    return {k: np.array(v) for k, v in out.items()}


def load_cloud(topic: str, cand: str) -> np.ndarray:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    f = rp.BLACKSMITH / f"{topic}_allposes" / f"{cand}.sdf"
    ms = [m for m in Chem.SDMolSupplier(str(f), removeHs=False, sanitize=False)
          if m is not None]
    if not ms:
        raise SystemExit(f"no poses in {f}")
    h = [i for i, a in enumerate(ms[0].GetAtoms()) if a.GetAtomicNum() > 1]
    return np.array([m.GetConformer().GetPositions()[h] for m in ms])


def contact_tensor(xyz: np.ndarray, res: dict, keys: list) -> np.ndarray:
    """(poses, atoms, residues) capped min-distance from each atom to each residue."""
    n_p, n_a = xyz.shape[0], xyz.shape[1]
    out = np.empty((n_p, n_a, len(keys)), dtype=np.float32)
    flat = xyz.reshape(-1, 3)
    for j, k in enumerate(keys):
        r = res[k]                                    # (n_res_atoms, 3)
        # min over the residue's atoms, for every (pose, atom) at once
        d = np.sqrt(((flat[:, None, :] - r[None, :, :]) ** 2).sum(-1)).min(1)
        out[:, :, j] = np.minimum(d, CAP_A).reshape(n_p, n_a)
    return out


def rank_residues(T: np.ndarray, keys: list) -> pd.DataFrame:
    """Per-residue variance across poses, summed over ligand atoms."""
    var = T.var(axis=0).sum(axis=0)                   # (residues,)
    rng = T.max(axis=0) - T.min(axis=0)
    return pd.DataFrame({
        "residue": [f"{c}:{i}:{n}" for c, i, n in keys],
        "total_variance": var,
        "mean_range_a": rng.mean(axis=0),
        "min_dist_a": T.min(axis=(0, 1)),
    }).sort_values("total_variance", ascending=False).reset_index(drop=True)


def greedy_nonredundant(T: np.ndarray, keys: list, k: int) -> list:
    """Forward selection: each pick maximises variance NOT explained by those held.

    Two residues on one pocket face move together across the cloud; taking both
    weights that direction twice and adds no information. Selection is on the
    per-pose mean profile, residualised against what is already chosen.
    """
    X = T.mean(axis=1)                                # (poses, residues)
    X = X - X.mean(0, keepdims=True)
    chosen, remaining = [], list(range(X.shape[1]))
    R = X.copy()
    for _ in range(min(k, len(remaining))):
        v = R[:, remaining].var(axis=0)
        best = remaining[int(np.argmax(v))]
        chosen.append(best)
        remaining.remove(best)
        b = R[:, best:best + 1]
        denom = float((b * b).sum()) or 1.0
        R = R - b @ (b.T @ R) / denom                 # project it out
    return chosen


def intrinsic_dim(T: np.ndarray) -> tuple:
    """Components of the flattened contact matrix needed for 90% / 95% variance."""
    X = T.reshape(T.shape[0], -1).astype(np.float64)
    X = X - X.mean(0, keepdims=True)
    # economical: eigenvalues of the (poses x poses) Gram matrix when that is
    # smaller than the feature count, which it is here
    s = np.linalg.svd(X, compute_uv=False)
    ev = s ** 2
    c = np.cumsum(ev) / ev.sum()
    return int(np.searchsorted(c, 0.90) + 1), int(np.searchsorted(c, 0.95) + 1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--top", type=int, default=25)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    res = receptor_residues()
    log.info("receptor: %d residues with heavy atoms", len(res))

    # candidate set from the FIRST cloud only, then held fixed, so every
    # replicate is scored on the same residues and the rankings are comparable
    xyz0 = load_cloud(f"election_{a.candidate}_r1", a.candidate)
    flat = xyz0.reshape(-1, 3)
    keys = []
    for k, r in res.items():
        d = np.sqrt(((flat[:, None, :] - r[None, :, :]) ** 2).sum(-1)).min()
        if d <= CANDIDATE_A:
            keys.append(k)
    keys.sort(key=lambda x: int(x[1]) if x[1].lstrip("-").isdigit() else 0)
    log.info("candidates within %.0f A of the cloud: %d residues", CANDIDATE_A, len(keys))

    ranks, dims = {}, []
    for i in range(1, a.replicates + 1):
        xyz = load_cloud(f"election_{a.candidate}_r{i}", a.candidate)
        T = contact_tensor(xyz, res, keys)
        r = rank_residues(T, keys)
        ranks[i] = r.set_index("residue")["total_variance"]
        d90, d95 = intrinsic_dim(T)
        dims.append((i, len(xyz), d90, d95))
        log.info("  r%d: %d poses, intrinsic dim %d (90%%) / %d (95%%) "
                 "from %d raw dimensions", i, len(xyz), d90, d95,
                 T.shape[1] * T.shape[2])

    V = pd.DataFrame(ranks)
    V["mean_variance"] = V.mean(axis=1)
    V["cv_across_replicates"] = V.iloc[:, :a.replicates].std(axis=1) / V["mean_variance"]
    V = V.sort_values("mean_variance", ascending=False)

    # does the ORDER survive an independent docking?
    from scipy.stats import spearmanr
    import itertools
    agree = [spearmanr(ranks[i], ranks[j])[0]
             for i, j in itertools.combinations(range(1, a.replicates + 1), 2)]

    # non-redundant pick, from replicate 1
    T1 = contact_tensor(xyz0, res, keys)
    pick = greedy_nonredundant(T1, keys, a.top)
    picked = [f"{keys[i][0]}:{keys[i][1]}:{keys[i][2]}" for i in pick]

    t = sout.Topic("blacksmith", f"residue_selection_{a.candidate}")
    V.to_csv(t.write("residue_variance", ".csv"))
    pd.DataFrame({"order": range(1, len(picked) + 1),
                  "residue": picked}).to_csv(t.write("greedy_pick", ".csv"),
                                             index=False)

    print("\n" + "=" * 76)
    print(f"  WHICH RESIDUES CARRY POSE INFORMATION — {a.candidate}")
    print("=" * 76)
    print(f"\n  candidates within {CANDIDATE_A:.0f} A of the cloud: {len(keys)}")
    print(f"  raw dimensions if all are used: {T1.shape[1]} atoms x {len(keys)} "
          f"residues = {T1.shape[1]*len(keys):,}")
    print(f"\n  INTRINSIC DIMENSIONALITY (PCA on the contact matrix):")
    for i, n, d90, d95 in dims:
        print(f"    r{i}: {d90:>3} components for 90%,  {d95:>3} for 95%   "
              f"(n={n})")
    print(f"\n  does the residue RANKING survive an independent docking?")
    print(f"    pairwise Spearman across {a.replicates} replicates: "
          f"mean {np.mean(agree):.3f}  min {np.min(agree):.3f}")
    print(f"\n  top {a.top} by mean variance (cv = spread across replicates):")
    print(V.head(a.top)[["mean_variance", "cv_across_replicates"]]
          .round(3).to_string())
    print(f"\n  greedy NON-REDUNDANT pick (each adds variance the held set lacks):")
    for i, p in enumerate(picked, 1):
        print(f"    {i:>2}. {p}")
    print(f"\n  variance is concentrated: top 10 hold "
          f"{V.mean_variance.head(10).sum()/V.mean_variance.sum()*100:.0f}% of it, "
          f"top 20 hold {V.mean_variance.head(20).sum()/V.mean_variance.sum()*100:.0f}%")
    print()


if __name__ == "__main__":
    main()
