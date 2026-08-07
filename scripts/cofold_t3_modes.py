"""
Purpose: T3 — does co-folding put sulfopin's sulfolane in the crystal mode or the decoy mode?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: the Boltz prediction for 6VAJ, the 20 docked sulfopin poses, the crystal ligand
Output: 00_outputs/blacksmith/cofold/cofold_t3_modes_<N>.csv

Issue #26: sulfopin's docked poses split in two. The sulfolane either sits in the
PROLINE pocket -- where the crystal puts it -- or against the BASIC cluster.
Docking energy prefers the wrong one, and by 0.16 kcal/mol, which is noise. If
something else could arbitrate that split, pose selection would stop being a
coin flip on this molecule and on every molecule like it.

THIS TEST IS ASYMMETRIC AND THE PRE-REGISTRATION SAYS SO. 6VAJ is IN Boltz-2's
training data (deposited 2019, cutoff 2023-06-01). A correct answer is therefore
nearly uninformative -- the model may simply be reciting 6VAJ. A WRONG answer is
strong: a model that cannot reproduce a pose it was trained on will not arbitrate
modes for molecules it has never seen. So this can rule the idea out and cannot
rule it in.

THE MODE MARKER IS THE SULFUR ATOM. Sulfopin has exactly one, in the sulfolane
ring, so the marker is unambiguous in an SDF, a PDB and a Boltz output alike --
no atom-name or SMARTS matching that could quietly pick a different group in one
of the three formats. Whole-molecule centroids do NOT separate the modes (they
differ by ~1 A); the sulfur separates them by 7.0 A.

EVERYTHING IS COMPARED IN ONE FRAME. The docked poses live in 3IKD's frame, the
crystal in 6VAJ's, and the prediction in its own. The crystal ligand is brought
into 3IKD's frame by the docking comparator, and the prediction is brought there
by superposing its protein. Nothing is compared across frames.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402

log = logging.getLogger("cofold-t3")
OUT = sout.Topic("blacksmith", "cofold")
WORK = Path("/data/lab_vm/modifiable/inhibition/cofold")
SULFOPIN_POSES = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
                      "nac_v2_poses/ref_Sulfopin__chloroacetamide.sdf")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(
        name, REPO / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sulfur_of_pdb(p: Path) -> np.ndarray:
    for ln in p.read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "S":
            return np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    raise SystemExit(f"no sulfur in {p}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--prediction", required=True, help="Boltz 6VAJ model_0.pdb")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    from rdkit import Chem
    cb = _load("cofold_bench")
    rb = _load("redock_3ikd_benchmark")
    rd02 = rb._rd02()
    receptor = rb.latest(sout.Topic("blacksmith", "receptor_3ikd"),
                         "3IKD_prepared", ".pdbqt")

    # ---- the two docked modes, in 3IKD's frame ---------------------------
    mols = [m for m in Chem.SDMolSupplier(str(SULFOPIN_POSES), removeHs=False) if m]
    S = np.array([m.GetConformer().GetPositions()[
        [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "S"][0]]
        for m in mols])
    ranks = [int(m.GetProp("pose_rank")) for m in mols]
    lab = fcluster(linkage(S, "average"), 2, "maxclust")
    modes = {k: {"n": int((lab == k).sum()),
                 "best_rank": min(r for r, l in zip(ranks, lab) if l == k),
                 "centroid": S[lab == k].mean(0)} for k in (1, 2)}

    # ---- the crystal, already transformed into 3IKD's frame ---------------
    crystal_ref = WORK / "docking" / "refs" / "6VAJ_ref.pdb"
    if not crystal_ref.is_file():
        raise SystemExit(f"run cofold_docking_comparator.py first ({crystal_ref})")
    s_crystal = sulfur_of_pdb(crystal_ref)

    # The crystal DEFINES which mode is correct. Naming a mode "proline" from
    # residue contacts would be an interpretation; distance to the deposited
    # sulfolane is a measurement.
    d_cry = {k: float(np.linalg.norm(s_crystal - v["centroid"])) for k, v in modes.items()}
    crystal_mode = min(d_cry, key=d_cry.get)
    decoy_mode = 3 - crystal_mode

    # ---- the prediction, brought into the same frame ----------------------
    pred = Path(args.prediction)
    p = cb.parse(pred)
    rec_atoms, _ = rd02.read_pdb(receptor)
    rec_ca = rd02.ca_map(rec_atoms)
    # rd02.ca_map is keyed by residue number; the prediction is numbered 1..N,
    # so pair by sequence alignment exactly as score_against does.
    # THREE-letter names from ca_map; see mode_arbitration._one_letter.
    def _ol(n):
        i = cb.gemmi.find_tabulated_residue(n)
        return (i.one_letter_code.upper() if i else "X") or "X"
    rec_seq = "".join(_ol(v[0]) for _, v in sorted(rec_ca.items()))
    rec_xyz = np.array([v[1] for _, v in sorted(rec_ca.items())])
    aln = cb.gemmi.align_string_sequences(list(p["ca_seq"]), list(rec_seq), [])
    ps, ts = aln.add_gaps(p["ca_seq"], 1), aln.add_gaps(rec_seq, 2)
    pi = ti = 0
    P, Q = [], []
    for a, b, m in zip(ps, ts, aln.match_string):
        if m == "|":
            P.append(p["ca_xyz"][pi]); Q.append(rec_xyz[ti])
        pi += a != "-"; ti += b != "-"
    if len(P) < 30:
        raise SystemExit(f"only {len(P)} CA pairs to 3IKD — cannot place the prediction")
    R, t = cb.kabsch(np.array(P), np.array(Q))
    s_pred_idx = [i for i, e in enumerate(p["el"]) if e == "S"]
    if not s_pred_idx:
        raise SystemExit("no sulfur in the predicted ligand")
    s_pred = (R @ p["lig"][s_pred_idx[0]]) + t

    d_pred = {k: float(np.linalg.norm(s_pred - v["centroid"])) for k, v in modes.items()}
    picked = min(d_pred, key=d_pred.get)
    d_pred_crystal = float(np.linalg.norm(s_pred - s_crystal))

    rows = [{"mode": k, "n_docked_poses": v["n"], "best_energy_rank": v["best_rank"],
             "is_crystal_mode": k == crystal_mode,
             "dist_to_crystal_S_A": round(d_cry[k], 2),
             "dist_to_prediction_S_A": round(d_pred[k], 2)} for k, v in modes.items()]
    df = pd.DataFrame(rows)
    dest = OUT.write("cofold_t3_modes", ".csv")
    df.to_csv(dest, index=False)

    sep = float(np.linalg.norm(modes[1]["centroid"] - modes[2]["centroid"]))
    print("\nT3 — sulfopin mode arbitration")
    print(f"  two docked modes, sulfur centroids {sep:.2f} A apart\n")
    print(df.to_string(index=False))
    print(f"\n  crystal mode  : {crystal_mode} "
          f"({modes[crystal_mode]['n']} poses, best energy rank "
          f"{modes[crystal_mode]['best_rank']})")
    print(f"  decoy mode    : {decoy_mode} "
          f"({modes[decoy_mode]['n']} poses, best energy rank "
          f"{modes[decoy_mode]['best_rank']}) — what docking energy prefers")
    print(f"\n  Boltz-2 lands in mode {picked} — "
          f"{'THE CRYSTAL MODE' if picked == crystal_mode else 'THE DECOY MODE'}")
    print(f"  predicted sulfur is {d_pred_crystal:.2f} A from the crystal sulfur")
    if picked == crystal_mode:
        print("\n  READING (fixed in advance): weak. 6VAJ is IN training, so this is")
        print("  consistent with recall rather than arbitration. It does NOT license")
        print("  using co-folding to split modes on unseen molecules.")
    else:
        print("\n  READING (fixed in advance): decisive against. A model that misses a")
        print("  pose it was TRAINED on will not arbitrate modes for molecules it has")
        print("  never seen. Mode splitting must be settled by physics.")
    print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
