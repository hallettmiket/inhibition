"""
Purpose: re-run the docking-enrichment gate on 3IKD. D0041's verdict was measured on 6VAJ and is invalid.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: covalent actives + the property-matched covalent decoy set
Output: 00_outputs/blacksmith/gate_3ikd/gate_3ikd_<N>.csv + the gate's own graded verdict

WHY THIS EXISTS. `config/gates.yaml` still says `receptor: 6VAJ`. D0041 measured
docking enrichment there and returned ROC-AUC 0.599, CI [0.311, 0.874], EF1% 0.0
— and that verdict is the entire justification for
`docs/ranking_rationale.md` discarding affinity and ranking on geometry instead.

**D0059 replaced the receptor and invalidated that measurement. Nobody re-ran
it.** Everything built afterwards inherited a conclusion whose evidence had been
retracted, which is D0069's subject.

WHAT THIS RUNS. The protocol the covalent arm actually uses on this branch —
plain AutoDock-GPU on the chemist-prepared 3IKD, no reactive potential — against
the property-matched decoy set the gate specifies. `config/gates.yaml` insists on
the EXACT downstream protocol rather than a proxy, and on this branch that is
AutoDock, not the Vina the config names.

THE VERDICT IS COMPUTED BY THE PROJECT'S OWN GATE, not here.
`shared.enrichment_gate.evaluate` applies the thresholds, counts independent
chemotypes, and grades. This script only supplies the scores. That matters:
a gate whose caller decides the verdict is not a gate.

EXPECT UNDERPOWERED. The config requires **6 independent chemotypes** for any
verdict above UNDERPOWERED and issue #12 §A establishes there are **3**. That is
a fact about Pin1's covalent literature, not a defect here, and the gate is
explicitly designed to report evidence strength rather than veto — so the point
estimates are the deliverable, read beside D0041's 0.599 on the old receptor.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import enrichment_gate as eg        # noqa: E402
from shared import outputs as sout              # noqa: E402
import nac_screen as ns                         # noqa: E402

log = logging.getLogger("gate-3ikd")

OUT = sout.Topic("blacksmith", "gate_3ikd")
PLAIN = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_plain")
DECOYS = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate/decoys")


def latest(pattern: str) -> Path:
    fs = [f for f in glob.glob(str(DECOYS / pattern)) if "report" not in f]
    if not fs:
        raise SystemExit(f"no decoys matching {pattern}")
    return Path(max(fs, key=lambda p: int(re.search(r"_(\d+)\.csv$", p).group(1))))


def build_set() -> pd.DataFrame:
    """Covalent actives (label 1) + property-matched covalent decoys (label 0).

    Two active sets are carried, and kept separate:

      `anchors` — what config/gates.yaml names: verified entries of
        pin1_covalent_cys113_anchors_2.csv. This is the gate as specified.
      `xtal` — the 17 ligands with an observed covalent bond to Cys113. Stronger
        positives (structural fact rather than curation), and the same set the
        geometric framework was validated against, so the two are comparable.
    """
    anchors = pd.read_csv(REPO / "data/reference/pin1_covalent_cys113_anchors_2.csv")
    anchors = anchors[anchors.smiles_status == "verified"]
    a1 = pd.DataFrame({"canonical_smiles": anchors.canonical_smiles,
                       "name": anchors.name, "label": 1, "active_set": "anchors"})

    links = pd.read_csv("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
                        "pdb_covalent/covalent_links_3.csv").drop_duplicates("comp_id")
    a2 = pd.DataFrame({"canonical_smiles": links.smiles,
                       "name": links.pdb_id + ":" + links.comp_id,
                       "label": 1, "active_set": "xtal"})

    dec = pd.read_csv(latest("decoys_covalent_*.csv"))
    smi = "canonical_smiles" if "canonical_smiles" in dec.columns else "smiles"
    d = pd.DataFrame({"canonical_smiles": dec[smi],
                      "name": dec.get("chembl_id", pd.Series(range(len(dec)))).astype(str),
                      "label": 0, "active_set": "decoy"})

    out = pd.concat([a1, a2, d], ignore_index=True).dropna(subset=["canonical_smiles"])
    out = out.drop_duplicates("canonical_smiles", keep="first").reset_index(drop=True)
    log.info("actives: %d anchors, %d crystallographic; decoys: %d",
             (out.active_set == "anchors").sum(), (out.active_set == "xtal").sum(),
             (out.active_set == "decoy").sum())
    return out


def dock_plain(smiles: str, nrun: int, gpu: str) -> float:
    """Best AutoDock energy over `nrun` independent runs. No reactive anything."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    RDLogger.DisableLog("rdApp.*")

    w = Path(tempfile.mkdtemp(prefix="gate_"))
    try:
        mol = ns.largest_fragment(smiles)
        if mol is None:
            raise ValueError("unparseable")
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
            raise ValueError("embed failed")
        AllChem.MMFFOptimizeMolecule(mol)
        txt, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation()(mol)[0])
        if not ok:
            raise ValueError(err)
        lig = w / "lig.pdbqt"
        lig.write_text(txt)
        subprocess.run(
            [str(ns.AUTODOCK), "-M", "rec.maps.fld", "-L", str(lig),
             "--nrun", str(nrun), "--resnam", str((w / "out").resolve())],
            cwd=PLAIN, check=True, capture_output=True,
            env=dict(os.environ, CUDA_VISIBLE_DEVICES=gpu))
        es = [e for e in ns.pose_energies(w / "out.dlg") if not np.isnan(e)]
        if not es:
            raise ValueError("no poses")
        return min(es)
    finally:
        shutil.rmtree(w, ignore_errors=True)


def run_shard(shard: int, n_shards: int, nrun: int, gpu: str) -> None:
    df = build_set()
    mine = df.iloc[shard::n_shards]
    rows = []
    for r in mine.itertuples():
        try:
            rows.append({"canonical_smiles": r.canonical_smiles, "name": r.name,
                         "label": r.label, "active_set": r.active_set,
                         "best_dg": dock_plain(r.canonical_smiles, nrun, gpu),
                         "status": "ok"})
        except Exception as exc:                       # noqa: BLE001
            rows.append({"canonical_smiles": r.canonical_smiles, "name": r.name,
                         "label": r.label, "active_set": r.active_set,
                         "best_dg": np.nan, "status": f"failed: {str(exc)[:90]}"})
    dest = OUT.write(f"gate_3ikd_s{shard}", ".csv")
    pd.DataFrame(rows).to_csv(dest, index=False)
    log.info("shard %d: %d scored -> %s", shard, len(rows), dest.name)


def report() -> None:
    fs = glob.glob(str(OUT.dir / "gate_3ikd_s*.csv"))
    if not fs:
        raise SystemExit("nothing scored yet")
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    df = df[df.status == "ok"].drop_duplicates("canonical_smiles")
    dec = df[df.label == 0]
    print(f"\n=== enrichment gate on 3IKD — plain AutoDock, no reactive potential ===")
    print(f"  {len(df)} molecules scored ({len(dec)} decoys)")

    for aset in ("anchors", "xtal"):
        act = df[(df.label == 1) & (df.active_set == aset)]
        if len(act) < 2:
            print(f"\n  {aset}: {len(act)} actives — too few"); continue
        sub = pd.concat([act, dec], ignore_index=True)
        # LOWER dG is better, so the gate is told so explicitly rather than the
        # sign being flipped here — a flipped sign that reaches the gate as
        # "higher is better" would invert the verdict and still look plausible.
        res = eg.evaluate(sub, metric="best_dg", stratum=f"covalent_{aset}",
                          higher_is_better=False)
        print(f"\n  --- actives = {aset} (n={len(act)}) vs {len(dec)} decoys ---")
        for f in ("roc_auc", "roc_auc_ci", "ef_1pct", "bedroc",
                  "n_chemotypes", "verdict", "reasons"):
            if hasattr(res, f):
                print(f"    {f:<14} {getattr(res, f)}")
    print(f"\n  D0041 measured ROC-AUC 0.599, CI [0.311, 0.874], EF1% 0.0 on 6VAJ.")
    print(f"  That verdict is invalid under D0059 and has never been replaced.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--nrun", type=int, default=2000)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format=f"%(levelname)s [g{args.shard}] %(message)s")
    if args.report:
        report()
        return
    run_shard(args.shard, args.n_shards, args.nrun, args.gpu)


if __name__ == "__main__":
    main()
