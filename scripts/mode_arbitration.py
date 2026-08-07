"""
Purpose: can a co-folding model pick WHICH binding mode is the real one?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: the 15 crystal complexes' 500-pose mode decompositions + their Boltz-2 predictions
Output: 00_outputs/blacksmith/pose_split/mode_arbitration_<N>.csv

@tt8804's proposal: between pose splitting and MD, turn each mode into a single
high-confidence pose, validate it against an independent structure predictor, and
if the mode-to-pose transition is weak, move down the ranked list rather than
spend four GPU-hours of MD on it.

THE CLAIM THIS TESTS is the load-bearing one: **does agreement with Boltz-2
identify the correct mode better than mode population does?** If it does, it is
an extremely cheap arbiter -- seconds per molecule against ~1 GPU-hour for BPMD
-- and it is ORTHOGONAL, because Boltz-2 never sees our poses. If it does not,
the mode ranking stays on consensus and this leg of the proposal is dropped.

WHY THIS IS ANSWERABLE TODAY. Both halves already exist for the same 15
molecules: 500 poses each with mode labels and RMSD to the deposited ligand, and
a Boltz-2 prediction. Nothing needs docking or folding again.

THREE ARBITERS ARE COMPARED, and the baseline is the one to beat:
  * consensus       -- the most populated mode (what the pipeline does now)
  * boltz           -- the mode whose poses come closest to Boltz-2's prediction
  * boltz x consensus -- population weighted by agreement

WHAT THIS CANNOT SETTLE. Boltz-2 predicts where a ligand SITS, not whether it
reacts, and the deposited ligand is a post-reaction adduct whose leaving group is
absent -- so "correct mode" here means "the mode containing the pose that best
overlays the adduct". That is the right target for placing a molecule in the
pocket and the wrong one for choosing between reaction trajectories. Stated here
because the same conflation already cost us one wrong conclusion today.
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                  # noqa: E402
from shared import pose_modes as pmod               # noqa: E402

log = logging.getLogger("mode-arb")
OUT = sout.Topic("blacksmith", "pose_split")
DUMP = Path("/data/lab_vm/modifiable/inhibition/keeprule/nrun500")
BOLTZ = Path("/data/lab_vm/modifiable/inhibition/cofold/t1t2")
SUCCESS_A = 2.0


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def _assign_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Optimal-assignment RMSD between two conformers of one molecule."""
    from scipy.optimize import linear_sum_assignment
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    r, c = linear_sum_assignment(d)
    return float(np.sqrt((d[r, c] ** 2).mean()))


def _one_letter(name: str) -> str:
    import gemmi
    info = gemmi.find_tabulated_residue(name)
    return (info.one_letter_code.upper() if info else "X") or "X"


def boltz_ligand_in_docking_frame(ident: str, cb, rb, rd02, receptor) -> np.ndarray | None:
    """Boltz-2's predicted ligand, moved into the receptor frame the poses live in.

    The prediction is in its own frame and the poses are in 3IKD's, so nothing can
    be compared until the proteins are superposed. Correspondence is by SEQUENCE
    alignment, never residue number: Boltz numbers its output 1..N and the
    prepared receptor uses author numbering.
    """
    hits = list((BOLTZ / ident).rglob("*model_0.pdb"))
    if not hits:
        return None
    p = cb.parse(hits[0])
    if not len(p["lig"]):
        return None
    rec_atoms, _ = rd02.read_pdb(receptor)
    ca = rd02.ca_map(rec_atoms)
    # ca_map yields THREE-letter residue names. Joining them produced a 339-char
    # "sequence" for a 113-residue protein, so the alignment walked off the end
    # of the coordinate array. Converted explicitly rather than sliced.
    rec_seq = "".join(_one_letter(v[0]) for _, v in sorted(ca.items()))
    rec_xyz = np.array([v[1] for _, v in sorted(ca.items())])
    aln = cb.gemmi.align_string_sequences(list(p["ca_seq"]), list(rec_seq), [])
    ps, ts = aln.add_gaps(p["ca_seq"], 1), aln.add_gaps(rec_seq, 2)
    pi = ti = 0
    P, Q = [], []
    for a, b, m in zip(ps, ts, aln.match_string):
        if m == "|":
            P.append(p["ca_xyz"][pi]); Q.append(rec_xyz[ti])
        pi += a != "-"; ti += b != "-"
    if len(P) < 30:
        return None
    R, t = cb.kabsch(np.array(P), np.array(Q))
    return (R @ p["lig"].T).T + t


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    cb = _load("cofold_bench")
    rb = _load("redock_3ikd_benchmark")
    rd02 = rb._rd02()
    receptor = rb.latest(sout.Topic("blacksmith", "receptor_3ikd"),
                         "3IKD_prepared", ".pdbqt")

    rows = []
    for f in sorted(DUMP.glob("*_modefeat.npy")):
        ident = f.name.replace("_modefeat.npy", "")
        feat = np.load(f)
        coords = np.load(DUMP / f"{ident}_coords.npy")       # (n_poses, n_atoms, 3)
        d = pd.read_csv(DUMP / f"{ident}_poses.csv")
        rmsd_x = d.rmsd_to_crystal.values
        labels = pmod.split(feat)
        modes = sorted(set(int(l) for l in labels) - {-1})
        if not modes:
            continue

        bl = boltz_ligand_in_docking_frame(ident, cb, rb, rd02, receptor)


        per = []
        for k in modes:
            sel = labels == k
            # FULL-POSE RMSD, NOT CENTROID DISTANCE.
            #
            # The first version of this compared ligand centroids and concluded
            # Boltz-2 arbitrates worse than consensus. That comparison could not
            # have worked: whole-molecule centroids barely separate modes at all
            # -- sulfopin's two modes differ by ~1 A in centroid while their
            # sulfur atoms sit 7 A apart -- so the argmin was choosing between
            # near-identical numbers. Both disagreements had the TRUE mode within
            # 0.3 A of the prediction, which is the signature of a metric with no
            # resolution rather than of a bad predictor.
            #
            # Optimal-assignment RMSD over all heavy atoms is element-blind here
            # (the docked poses' elements are not persisted) but both structures
            # are the same molecule with the same atom count, so an assignment
            # can only be permissive, never wrong in a way that favours a mode.
            dboltz = (float(np.min([_assign_rmsd(c, bl) for c in coords[sel]]))
                      if bl is not None and coords[sel].shape[1] == len(bl)
                      else np.nan)
            per.append({
                "mode": k, "consensus": float(sel.mean()),
                "n": int(sel.sum()),
                "dist_to_boltz_a": dboltz,
                "best_rmsd_to_crystal": float(np.nanmin(rmsd_x[sel])),
                "contains_crystal": bool(np.nanmin(rmsd_x[sel]) <= SUCCESS_A),
                "spread_a": pmod.identity(feat, labels, k)["spread_a"],
                "dir_coherence": pmod.identity(feat, labels, k)["dir_coherence"],
            })
        pm = pd.DataFrame(per)

        # the three arbiters
        pick_cons = int(pm.loc[pm.consensus.idxmax(), "mode"])
        pick_boltz = (int(pm.loc[pm.dist_to_boltz_a.idxmin(), "mode"])
                      if pm.dist_to_boltz_a.notna().any() else None)
        w = pm.consensus / (1.0 + pm.dist_to_boltz_a)
        pick_both = int(pm.loc[w.idxmax(), "mode"]) if w.notna().any() else None
        truth = int(pm.loc[pm.best_rmsd_to_crystal.idxmin(), "mode"])

        rows.append({
            "ident": ident, "n_modes": len(modes),
            "true_mode": truth,
            "pick_consensus": pick_cons, "pick_boltz": pick_boltz,
            "pick_both": pick_both,
            "hit_consensus": pick_cons == truth,
            "hit_boltz": pick_boltz == truth,
            "hit_both": pick_both == truth,
            "any_mode_contains_crystal": bool(pm.contains_crystal.any()),
            "boltz_dist_of_true_mode": float(
                pm.loc[pm["mode"] == truth, "dist_to_boltz_a"].iloc[0]),
        })
        log.info("%s: %d modes, true=%d, consensus=%s boltz=%s",
                 ident, len(modes), truth, pick_cons, pick_boltz)

    t = pd.DataFrame(rows)
    dest = OUT.write("mode_arbitration", ".csv")
    t.to_csv(dest, index=False)

    print(f"\nWhich arbiter picks the mode containing the crystal pose?  (n={len(t)})\n")
    for c, lab in (("hit_consensus", "most populated mode (current)"),
                   ("hit_boltz", "closest to Boltz-2"),
                   ("hit_both", "consensus / (1 + distance to Boltz-2)")):
        print(f"  {lab:<38} {t[c].mean()*100:5.1f}%")
    print(f"\n  ceiling (some mode contains it): {t.any_mode_contains_crystal.mean()*100:.1f}%")
    print(f"  mean modes per molecule: {t.n_modes.mean():.1f}")
    print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
