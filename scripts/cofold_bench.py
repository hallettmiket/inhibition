"""
Purpose: run Boltz-2 on Pin1 complexes and score its poses against ground truth we hold.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: a set of (ident, SMILES) + optionally a deposited structure to score against
Output: 00_outputs/blacksmith/cofold/cofold_<TAG>_<N>.csv + Boltz's own predictions

Runs the experiment pre-registered in `docs/prereg_cofolding.md`. Nothing here
decides anything; it produces the numbers that document's readings are attached
to.

WHY A CO-FOLDING MODEL AT ALL. Every signal this project has tested -- docking
energy, enrichment, consensus, top-N viability, anchor quality -- is computed
from the same AutoDock pose set and inherits its failure modes (#23). Boltz-2
builds the complex from sequence and ligand and never sees those poses, so
whatever it says is ORTHOGONAL. That is the property being tested, not accuracy
for its own sake.

THE MSA IS COMPUTED ONCE. The protein is always Pin1, so its MSA is generated a
single time and reused for every ligand, collapsing per-candidate cost to one
forward pass. Without this the MSA dominates and the entire cheap-filter argument
fails -- so the cache is not an optimisation, it is the premise.

CONTAMINATION IS TRACKED PER ROW, NOT ASSUMED. Co-folding models train on the
PDB, so reproducing a pose from a structure in training measures memorisation.
Every row carries `held_out`, set from the deposition era, and the two sets are
never pooled into one accuracy number.

SCORING IS SYMMETRY-CORRECTED AND SUPERPOSITION-BASED. The prediction is in its
own frame, so protein CA atoms are superposed onto the deposited structure and
the transform applied to the predicted ligand before any RMSD is taken. A
symmetry-naive RMSD would penalise a correct pose of a molecule with equivalent
atoms -- the same trap `redock_3ikd_benchmark` documents.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                # noqa: E402

log = logging.getLogger("cofold")
OUT = sout.Topic("blacksmith", "cofold")
BOLTZ = Path.home() / ".micromamba/envs/dwi_boltz/bin/boltz"
WORK = Path("/data/lab_vm/modifiable/inhibition/cofold")
MSA_CACHE = WORK / "pin1_msa"

#: Pin1 PPIase domain as it appears in 3IKD_ian, residues 51-163.
PIN1 = ("EPARVRCSHLLVKHSQSRRPSSWRQEQITRTQEEALELINGYIQKIKSGEEDFESLASQFSDCSSAKARG"
        "DLGAFSRGQMQKPFEDASFALRTGEMSGPVFTDSGIHIILRTE")


def ensure_msa(seq: str) -> Path | None:
    """The MSA for one protein sequence, computed at most once ever.

    This is the premise of the whole cost argument, not an optimisation. MSA
    generation dominates a co-folding run; the forward pass does not. Because the
    protein is fixed and only the ligand varies, one MSA serves every candidate
    and a 300-molecule triage costs 300 forward passes plus ONE search. Cached on
    a hash of the sequence, so a changed construct gets its own MSA rather than
    silently reusing the wrong one.
    """
    MSA_CACHE.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256(seq.encode()).hexdigest()[:16]
    p = MSA_CACHE / f"{h}.csv"
    if p.is_file() and p.stat().st_size > 100:
        return p
    from boltz.main import compute_msa
    log.info("computing MSA for %d-residue construct %s (once)", len(seq), h)
    hits: list[Path] = []
    for attempt in range(1, 4):
        tmp = MSA_CACHE / f"tmp_{h}_{attempt}"
        tmp.mkdir(parents=True, exist_ok=True)
        try:
            compute_msa(data={h: seq}, target_id=h, msa_dir=tmp,
                        msa_server_url="https://api.colabfold.com",
                        msa_pairing_strategy="greedy")
        except Exception as e:                         # noqa: BLE001
            # The public ColabFold server occasionally returns a truncated
            # archive. That is transient, so retry — but never fall through to
            # single-sequence mode silently, because an MSA-less prediction is a
            # much weaker model and would be read as a Boltz-2 result.
            log.warning("MSA attempt %d/3 failed for %s: %s", attempt, h, e)
            continue
        hits = list(tmp.glob("*.csv"))
        if hits:
            break
    if not hits:
        log.error("no MSA for construct %s after 3 attempts", h)
        return None
    p.write_bytes(hits[0].read_bytes())
    log.info("cached MSA %s (%d lines)", p.name, len(p.read_text().splitlines()))
    return p


def yaml_for(seq: str, smiles: str, msa: Path | None) -> str:
    """Boltz input: the entry's own protein construct + its ligand from SMILES."""
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {seq}\n"
        f"      msa: {msa if msa else 'empty'}\n"
        "  - ligand:\n"
        "      id: B\n"
        f"      smiles: '{smiles}'\n"
    )


def run_boltz(ident: str, seq: str, smiles: str, outdir: Path, gpu: int,
              msa: Path | None, steps: int) -> Path | None:
    outdir.mkdir(parents=True, exist_ok=True)
    y = outdir / f"{ident}.yaml"
    y.write_text(yaml_for(seq, smiles, msa))
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(gpu))
    cmd = [str(BOLTZ), "predict", str(y), "--out_dir", str(outdir),
           "--output_format", "pdb", "--diffusion_samples", "1",
           "--sampling_steps", str(steps), "--override", "--no_kernels"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=7200)
    if r.returncode != 0:
        log.warning("%s: boltz rc=%d %s", ident, r.returncode, r.stderr[-400:])
        return None
    hits = sorted(outdir.rglob("*model_0.pdb"))
    return hits[0] if hits else None


def confidence(pred: Path) -> dict:
    """Boltz's own confidence for this prediction, from the JSON beside it."""
    out = {}
    for j in list(pred.parent.glob("confidence*.json")):
        try:
            d = json.loads(j.read_text())
        except Exception:                          # noqa: BLE001
            continue
        for k in ("confidence_score", "ptm", "iptm", "ligand_iptm",
                  "protein_iptm", "complex_plddt", "complex_iplddt",
                  "complex_pde", "complex_ipde"):
            if k in d:
                out[k] = d[k]
    return out


def parse(path: Path) -> dict:
    """Ligand atoms (coords + elements) and the CA trace with its sequence.

    The CA trace carries the one-letter sequence because residue NUMBERS cannot
    be trusted to correspond: Boltz numbers its output 1..N from the sequence it
    was given, while a deposited structure uses author numbering (Pin1's is
    51..163). Matching on number would silently superpose the wrong residues.
    """
    lig_xyz, lig_el, ca = [], [], []
    for l in path.read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        name, resn = l[12:16].strip(), l[17:20].strip()
        el = (l[76:78].strip() or name[0]).upper()
        if el == "H":
            continue
        xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
        if name == "CA" and gemmi_is_aa(resn):
            ca.append((int(l[22:26]), gemmi.find_tabulated_residue(resn).one_letter_code.upper(), xyz))
        elif resn not in ("HOH", "WAT") and not gemmi_is_aa(resn):
            lig_xyz.append(xyz)
            lig_el.append(el)
    return {"lig": np.array(lig_xyz), "el": np.array(lig_el),
            "ca_seq": "".join(c[1] for c in ca),
            "ca_xyz": np.array([c[2] for c in ca])}


def gemmi_is_aa(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return bool(info and info.is_amino_acid())


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Rotation+translation taking P onto Q (both n x 3, paired)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def sym_rmsd(a: np.ndarray, ael: np.ndarray,
             b: np.ndarray, bel: np.ndarray) -> float:
    """Element-aware assignment RMSD over the atoms the two molecules share.

    Two corrections, both of which change the number materially here:

    SYMMETRY. An atom-order RMSD punishes a correct pose of any molecule with
    equivalent atoms — a flipped phenyl scores as a miss. Optimal assignment
    removes that, as `redock_3ikd_benchmark` already does for docking.

    UNEQUAL ATOM COUNTS. These ligands are deposited as ADDUCTS and have lost
    their leaving group (sulfopin's Cl, and so on), while the SMILES Boltz is
    given is the intact free compound. The cost matrix is therefore rectangular
    and every truth atom is matched to a distinct predicted atom; the surplus
    predicted atoms are simply unmatched. Truncating the longer list by index
    instead would drop atoms by file order, which is arbitrary.

    Assignment is restricted WITHIN element, so a carbon can never be scored
    against an oxygen to flatter the result.
    """
    from scipy.optimize import linear_sum_assignment
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    d = np.where(ael[:, None] == bel[None, :], d, 1e6)
    r, c = linear_sum_assignment(d)
    keep = d[r, c] < 1e5
    if keep.sum() < 3:
        return float("nan")
    return float(np.sqrt((d[r, c][keep] ** 2).mean()))


def score_against(pred: Path, truth: Path) -> dict:
    """Superpose predicted protein onto the deposited one, then compare ligands.

    The prediction is in its own frame, so nothing can be compared until the
    proteins are on top of each other. Correspondence comes from a SEQUENCE
    alignment, not from residue numbers, because the two files number the same
    protein differently.
    """
    p, t = parse(pred), parse(truth)
    if not len(p["lig"]) or not len(t["lig"]):
        return {"rmsd_A": np.nan, "note": "no ligand atoms parsed"}

    aln = gemmi.align_string_sequences(list(p["ca_seq"]), list(t["ca_seq"]), [])
    ps, ts = aln.add_gaps(p["ca_seq"], 1), aln.add_gaps(t["ca_seq"], 2)
    pi = ti = 0
    P, Q = [], []
    for a, b, m in zip(ps, ts, aln.match_string):
        if m == "|":                      # identical residue in both — pair the CAs
            P.append(p["ca_xyz"][pi])
            Q.append(t["ca_xyz"][ti])
        pi += a != "-"
        ti += b != "-"
    if len(P) < 30:
        return {"rmsd_A": np.nan, "note": f"only {len(P)} aligned CA"}

    R, tr = kabsch(np.array(P), np.array(Q))
    moved = (R @ p["lig"].T).T + tr
    ca_rms = float(np.sqrt((((R @ np.array(P).T).T + tr - np.array(Q)) ** 2)
                           .sum(1).mean()))
    return {"rmsd_A": sym_rmsd(moved, p["el"], t["lig"], t["el"]),
            "n_ca_superposed": len(P), "ca_rmsd_A": ca_rms,
            "n_lig_atoms_pred": len(p["lig"]), "n_lig_atoms_truth": len(t["lig"])}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--set", required=True, help="CSV: ident,smiles[,truth,held_out]")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--gpu", type=int, default=3)
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--msa", default=None, help="cached Pin1 MSA (.a3m)")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if not BOLTZ.is_file():
        raise SystemExit(f"boltz not found at {BOLTZ}")
    df = pd.read_csv(args.set)
    if args.limit:
        df = df.head(args.limit)
    msa = Path(args.msa) if args.msa else None
    log.info("%d molecules, GPU %d, msa=%s", len(df), args.gpu, msa or "generated per-run")

    rows = []
    for i, r in enumerate(df.itertuples(), 1):
        wd = WORK / args.tag / str(r.ident).replace(":", "_")
        seq = getattr(r, "sequence", None) or PIN1
        rec = {"ident": r.ident, "smiles": r.smiles, "seq_len": len(seq),
               "held_out": getattr(r, "held_out", None)}
        rmsa = msa if msa else ensure_msa(seq)
        rec["msa"] = str(rmsa) if rmsa else "NONE — single-sequence mode"
        pred = run_boltz(str(r.ident).replace(":", "_"), seq, r.smiles, wd,
                         args.gpu, rmsa, args.steps)
        if pred is None:
            rec["status"] = "boltz failed"
        else:
            rec["status"] = "ok"
            rec["prediction"] = str(pred)
            rec.update(confidence(pred))
            truth = getattr(r, "truth", None)
            if isinstance(truth, str) and Path(truth).is_file():
                rec.update(score_against(pred, Path(truth)))
        rows.append(rec)
        log.info("[%d/%d] %s: %s  rmsd=%s", i, len(df), r.ident, rec["status"],
                 f"{rec.get('rmsd_A', float('nan')):.2f}" if rec.get("rmsd_A") == rec.get("rmsd_A") else "-")
        pd.DataFrame(rows).to_csv(OUT.write(f"cofold_{args.tag}", ".csv"), index=False)

    d = pd.DataFrame(rows)
    ok = d[d.status == "ok"]
    print(f"\n=== {len(ok)}/{len(d)} predicted ===")
    if "rmsd_A" in ok.columns and ok.rmsd_A.notna().any():
        for grp, g in ok.dropna(subset=["rmsd_A"]).groupby("held_out", dropna=False):
            lab = {True: "HELD OUT", False: "in training"}.get(grp, str(grp))
            print(f"  {lab:<12} n={len(g):<3} within 2 A: {(g.rmsd_A<=2).mean()*100:5.1f}%"
                  f"   median {g.rmsd_A.median():.2f} A")
        print("\n  The two sets are NOT pooled: reproducing a pose from a structure in\n"
              "  training measures memorisation, not prediction.")


if __name__ == "__main__":
    main()
