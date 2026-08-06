"""
Purpose: does pose consensus converge where the viable-NAC frequency did not? The falsifier.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: the crystallographic Cys113 positives, docked at two search efforts
Output: 00_outputs/blacksmith/consensus_convergence/*.csv + a verdict

THE CLAIM UNDER TEST, STATED SO IT CAN FAIL.

`shared/pose_consensus.py` was built on the argument that agreement among the
top-N poses should be more stable than the viable-NAC *fraction*, because
frequency is computed over EVERY run and therefore inherits how hard the search
looked, while a top-N window does not.

D0068 measured the frequency failing exactly that way: at 10x the search effort
the same molecules fell from 2.91x to 0.96x, and the crystallographic positives
— never selected on score — fell identically. The argument for consensus is
plausible and, until this script runs, **untested**. Its own module docstring
says so.

    IF CONSENSUS ALSO MOVES between 200 and 2,000 runs, it inherits the search
    too, and it is not the fix for D0068. It may still be a useful component;
    it is not a more stable one, and the rationale must stop claiming otherwise.

    IF IT HOLDS while the frequency on the SAME dockings moves, that is direct
    evidence for the design, measured rather than argued.

WHY BOTH ARE MEASURED FROM THE SAME .dlg. Frequency and consensus are computed
from one docking run at each effort, not two. Comparing metrics across separate
runs would confound the metric's stability with the run-to-run scatter, and
run-to-run scatter is precisely the thing under examination.

TOP-N IS HELD FIXED ACROSS EFFORTS, because `require_same_n` exists for this
reason: a consensus at N=10 and one at N=50 are not the same quantity, and
letting N float with the run count would compare two different measurements and
call the difference convergence.
"""

from __future__ import annotations

import argparse
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

log = logging.getLogger("consensus-conv")
OUT = sout.Topic("blacksmith", "consensus_convergence")


def poses_and_fraction(cand: ns.Candidate, rec_dir: Path, nrun: int,
                       gpu: str) -> tuple[list[pc.ReactivePose], float]:
    """One docking run, reduced to both metrics — same poses, same file."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem

    work = Path(tempfile.mkdtemp(prefix="conv_"))
    try:
        best: tuple | None = None
        for j, lig in enumerate(ns.prepare_ligand(cand, work / "lig.pdbqt")):
            dlg = ns.dock(lig, rec_dir, work / f"c{j}", nrun, gpu)
            res = ns.measure_dlg(dlg, cand)
            energies = ns.pose_energies(dlg)
            if len(energies) != len(res):
                raise ValueError("energy/geometry length mismatch")
            pm = PDBQTMolecule.from_file(str(dlg), is_dlg=True, skip_typing=True)
            mols = [m for m in RDKitMolCreate.from_pdbqt_mol(pm) if m is not None]
            if not mols:
                continue
            mol = mols[0]
            match = mol.GetSubstructMatches(
                Chem.MolFromSmarts(cand.reactive_smarts))
            if not match:
                continue
            idx = list(match[0])
            poses = []
            for k, e in enumerate(energies):
                if np.isnan(e):
                    continue
                xyz = mol.GetConformer(k).GetPositions()[idx]
                poses.append(pc.ReactivePose(energy=float(e),
                                             reactive_xyz=np.asarray(xyz),
                                             atom_ids=tuple(idx)))
            frac = nac.viable_fraction(res)
            if best is None or frac > best[1]:
                best = (poses, frac)
        if best is None:
            raise ValueError("no usable poses")
        return best
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--top-n", type=int, default=10,
                    help="held FIXED across efforts; see the module docstring")
    ap.add_argument("--efforts", type=int, nargs=2, default=[200, 2000])
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    meta = {r.class_id: (r.mechanism, r.reactive_atom_smarts) for r in wh.itertuples()}
    cands = ns.crystal_positives(meta, None)
    if args.limit:
        cands = cands[:args.limit]
    rec = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    lo, hi = args.efforts
    log.info("%d crystallographic positives at %d and %d runs, top_n=%d",
             len(cands), lo, hi, args.top_n)

    rows = []
    for i, c in enumerate(cands, 1):
        row = {"ident": c.ident, "warhead_class": c.warhead_class}
        try:
            for nrun in (lo, hi):
                poses, frac = poses_and_fraction(c, rec, nrun, args.gpu)
                row[f"fraction_{nrun}"] = frac
                if len(poses) >= pc.MIN_POSES_FOR_CONSENSUS:
                    r = pc.consensus(poses, top_n=args.top_n)
                    row[f"consensus_{nrun}"] = r.agreement
                    row[f"median_rmsd_{nrun}"] = r.median_rmsd
                else:
                    row[f"consensus_{nrun}"] = np.nan
            row["status"] = "ok"
            log.info("[%d/%d] %-22s frac %.3f->%.3f   consensus %.3f->%.3f",
                     i, len(cands), c.ident[:22], row[f"fraction_{lo}"],
                     row[f"fraction_{hi}"], row[f"consensus_{lo}"],
                     row[f"consensus_{hi}"])
        except Exception as exc:                       # noqa: BLE001
            row["status"] = f"failed: {str(exc)[:110]}"
            log.warning("[%d/%d] %s: %s", i, len(cands), c.ident, row["status"])
        rows.append(row)

    df = pd.DataFrame(rows)
    dest = OUT.write("consensus_convergence", ".csv")
    df.to_csv(dest, index=False)

    ok = df[df.status == "ok"]
    if ok.empty:
        print("\n  nothing measured"); return
    fl, fh = ok[f"fraction_{lo}"], ok[f"fraction_{hi}"]
    cl, ch = ok[f"consensus_{lo}"], ok[f"consensus_{hi}"]
    print(f"\n=== consensus vs frequency, {len(ok)} crystallographic positives ===")
    print(f"  {'':<12}{lo:>10} runs{hi:>10} runs{'median |change|':>18}")
    print(f"  {'frequency':<12}{fl.median():>15.3f}{fh.median():>15.3f}"
          f"{(fh - fl).abs().median():>18.3f}")
    print(f"  {'consensus':<12}{cl.median():>15.3f}{ch.median():>15.3f}"
          f"{(ch - cl).abs().median():>18.3f}")
    from scipy.stats import spearmanr, wilcoxon
    print(f"\n  rank agreement across efforts (higher = more stable):")
    print(f"    frequency  Spearman rho = {spearmanr(fl, fh).statistic:+.3f}")
    print(f"    consensus  Spearman rho = {spearmanr(cl, ch).statistic:+.3f}")
    try:
        print(f"\n  does each shift systematically? (Wilcoxon, paired)")
        print(f"    frequency p = {wilcoxon(fl, fh).pvalue:.4f}")
        print(f"    consensus p = {wilcoxon(cl, ch).pvalue:.4f}")
    except ValueError as exc:
        print(f"    unavailable: {exc}")
    print("\n  READING, fixed in advance: consensus is the fix for D0068 only if it")
    print("  moves LESS than the frequency on these same dockings. If it moves as")
    print("  much, it inherits the search too and the rationale must say so.")


if __name__ == "__main__":
    main()
