"""
Purpose: consensus for EVERY candidate, so the composite ranking has two components everywhere.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: all scorable T_3/T_4 candidates + the reactive 3IKD receptor
Output: 00_outputs/blacksmith/nac_consensus/nac_consensus_s<shard>_<N>.csv

WHY A SECOND PASS OVER THE SAME MOLECULES.

`nac_rank` docked all 5,769 and kept only the summary — it discarded the poses,
so consensus (which is a statement about the poses) cannot be computed from its
output. Consensus currently exists for the 50 shortlisted candidates only,
because `export_nac_poses` re-docked those.

@tt8804 wants "a ranked list that we can query" over everything, with nothing
deleted. `shared/composite_rank.py` treats an unmeasured component as the
interval [0, 1] — total ignorance — so a candidate missing consensus is ranked
honestly rather than wrongly. But a component measured on 50 of 5,769 is barely a
component, and the coverage gap would show up as ~5,700 candidates all carrying
the same maximal uncertainty.

BOTH METRICS COME FROM ONE DOCKING RUN, not two. Enrichment is recomputed here
alongside consensus rather than joined from `nac_rank`, because a consensus from
this run and an enrichment from that one describe different dockings of the same
molecule — and D0068 established that two dockings of one molecule do not agree.
Joining them would silently pair a pose set with a frequency it did not produce.

RUN COUNT IS PART OF BOTH NUMBERS (D0068, D0070) and is written into every row.
`top_n` likewise, because `pose_consensus.require_same_n` refuses to compare
agreements measured at different N — and this is exactly the caller its docstring
warns about.

POSES ARE NOT KEPT. 5,769 multi-model SDFs would be a large write for data that
is only ever reduced to two numbers. The 50 shortlisted candidates keep their
poses (`export_nac_poses`) because those are the ones the GUI draws.
"""

from __future__ import annotations

import argparse
import glob
import logging
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import nac_criterion as nac          # noqa: E402
from shared import outputs as sout               # noqa: E402
from shared import pose_consensus as pc          # noqa: E402
import nac_screen as ns                          # noqa: E402
import nac_rank as nr                            # noqa: E402

log = logging.getLogger("nac-consensus-all")
OUT = sout.Topic("blacksmith", "nac_consensus")


def score_one(cand: ns.Candidate, rec_dir: Path, nrun: int, gpu: str,
              top_n: int) -> dict:
    """Enrichment AND consensus from one docking run of one molecule."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem

    work = Path(tempfile.mkdtemp(prefix="nacall_"))
    row = {"ident": cand.ident, "warhead_class": cand.warhead_class,
           "mechanism": cand.mechanism, "nrun": nrun, "top_n": top_n}
    try:
        best = None
        for j, lig in enumerate(ns.prepare_ligand(cand, work / "lig.pdbqt")):
            dlg = ns.dock(lig, rec_dir, work / f"c{j}", nrun, gpu)
            res = ns.measure_dlg(dlg, cand)
            energies = ns.pose_energies(dlg)
            if len(energies) != len(res):
                raise ValueError("energy/geometry length mismatch")
            frac = nac.viable_fraction(res)
            if best is not None and frac <= best[0]:
                continue
            pm = PDBQTMolecule.from_file(str(dlg), is_dlg=True, skip_typing=True)
            mols = [m for m in RDKitMolCreate.from_pdbqt_mol(pm) if m is not None]
            if not mols:
                continue
            mol = mols[0]
            match = mol.GetSubstructMatches(Chem.MolFromSmarts(cand.reactive_smarts))
            if not match:
                continue
            idx = list(match[0])
            poses = [pc.ReactivePose(
                        energy=float(e),
                        reactive_xyz=np.asarray(mol.GetConformer(k).GetPositions()[idx]),
                        atom_ids=tuple(idx))
                     for k, e in enumerate(energies) if not np.isnan(e)]
            best = (frac, poses, res)

        if best is None:
            raise ValueError("no usable poses")
        frac, poses, res = best
        row["viable_fraction"] = frac
        row["enrichment"] = frac / nac.isotropic_null(cand.mechanism)
        row["n_poses"] = len(poses)

        if len(poses) >= pc.MIN_POSES_FOR_CONSENSUS:
            c = pc.consensus(poses, top_n=min(top_n, len(poses)))
            row.update({"consensus": c.agreement,
                        "consensus_n": c.n_poses,
                        "consensus_lo": c.agreement_jackknife[0],
                        "consensus_hi": c.agreement_jackknife[1],
                        "consensus_median_rmsd": c.median_rmsd,
                        "consensus_modes": c.n_modes})
        else:
            # Recorded as unmeasurable, never as zero agreement.
            row["consensus"] = np.nan
            row["consensus_note"] = f"only {len(poses)} poses"
        row["status"] = "ok"
    except Exception as exc:                            # noqa: BLE001
        row["status"] = f"failed: {str(exc)[:130]}"
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return row


def done_idents() -> set[str]:
    ids = set()
    for f in glob.glob(str(OUT.dir / "nac_consensus_s*.csv")):
        try:
            ids.update(pd.read_csv(f).ident.astype(str))
        except Exception:                               # noqa: BLE001
            pass
    return ids


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format=f"%(levelname)s [c{args.shard}] %(message)s")

    rec = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    cands = [c for i, c in enumerate(nr.load_candidates())
             if i % args.n_shards == args.shard]
    if args.limit:
        cands = cands[:args.limit]
    done = done_idents()
    todo = [c for c in cands if c.ident not in done]
    log.info("shard %d/%d: %d assigned, %d already done, %d to do",
             args.shard, args.n_shards, len(cands), len(cands) - len(todo), len(todo))

    buf = []
    for k, c in enumerate(todo, 1):
        buf.append(score_one(c, rec, args.nrun, args.gpu, args.top_n))
        if len(buf) >= args.chunk or k == len(todo):
            dest = OUT.write(f"nac_consensus_s{args.shard}", ".csv")
            pd.DataFrame(buf).to_csv(dest, index=False)
            ok = sum(r["status"] == "ok" for r in buf)
            log.info("shard %d: %d/%d, wrote %d (%d ok) -> %s",
                     args.shard, k, len(todo), len(buf), ok, dest.name)
            buf = []


if __name__ == "__main__":
    main()
