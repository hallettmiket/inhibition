"""
Purpose: score the TOP-N poses instead of the whole run population, and see which survives more searching.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: crystallographic positives + siblings, docked at two search efforts
Output: 00_outputs/blacksmith/topn_vs_fraction/*.csv + a verdict

@tt8804's diagnosis, which is better than the one it replaces:

    "the tool makes the best poses in order right? we just want to see the top
     few poses, you are over diluting it by running something that can be rerun
     infinitely."

THE DESIGN ERROR THIS TESTS. `nac_criterion.viable_fraction` divides by EVERY
run. AutoDock returns poses ranked by energy, so the question a chemist asks —
"do this molecule's BEST poses present the warhead?" — is a question about the
top of that ranking. Asking instead "what fraction of all 200 attempts landed in
position" puts every mediocre pose in the denominator, and every additional run
adds more of them. In the limit the fraction approaches the pocket's background
rate for anything, independent of the molecule. That is not a property of the
molecule, and D0068 measured exactly that collapse without identifying the cause.

THE HINT ALREADY IN THE DATA. `pose_consensus` scores the TOP-N by energy, and
D0070 measured it preserving rank across efforts (Spearman +0.568) where the
whole-population frequency did not (-0.047). Same dockings, same molecules — the
only difference is the window. That should have been read as evidence about the
WINDOW rather than about consensus.

WHAT IS COMPARED, on identical poses from one docking run at each effort:

    fraction_all      viable / all runs                 (the current metric)
    topN_viable       viable among the N best by energy (the proposed metric)

If `topN_viable` holds where `fraction_all` collapses, the metric was diluting a
real signal and the fix is a redefinition, not 10x the compute.

TOP-N IS SCORED ON THE MOLECULE'S OWN BEST POSES, so it is not a per-run rate and
adding runs cannot dilute it — more searching can only IMPROVE which poses occupy
the top N. That is the property the ranking needed and never had.
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
import nac_screen as ns                          # noqa: E402
import nac_rank as nr                            # noqa: E402

log = logging.getLogger("topn-vs-fraction")
OUT = sout.Topic("blacksmith", "topn_vs_fraction")
TOP_NS = (1, 3, 5, 10, 20)


def measure(cand: ns.Candidate, rec_dir: Path, nrun: int, gpu: str) -> dict:
    """Both metrics from ONE docking run, so the comparison isolates the window."""
    work = Path(tempfile.mkdtemp(prefix="topn_"))
    row = {"ident": cand.ident, "warhead_class": cand.warhead_class, "nrun": nrun}
    try:
        best = None
        for j, lig in enumerate(ns.prepare_ligand(cand, work / "lig.pdbqt")):
            dlg = ns.dock(lig, rec_dir, work / f"c{j}", nrun, gpu)
            res = ns.measure_dlg(dlg, cand)
            energies = ns.pose_energies(dlg)
            if len(energies) != len(res):
                raise ValueError("length mismatch")
            frac = nac.viable_fraction(res)
            if best is None or frac > best[0]:
                best = (frac, res, energies)
        if best is None:
            raise ValueError("no usable poses")
        frac, res, energies = best

        row["fraction_all"] = frac
        row["enrichment_all"] = frac / nac.isotropic_null(cand.mechanism)
        # Rank poses by energy — LOWER is better — then ask about the top of that
        # ranking, which is the list a docking program actually hands a chemist.
        order = sorted((i for i, e in enumerate(energies) if not np.isnan(e)),
                       key=lambda i: energies[i])
        row["n_scored"] = len(order)
        for n in TOP_NS:
            top = order[:n]
            if not top:
                continue
            row[f"top{n}_viable"] = sum(res[i].viable for i in top) / len(top)
            row[f"top{n}_best_dist"] = min(res[i].distance for i in top)
        row["status"] = "ok"
    except Exception as exc:                            # noqa: BLE001
        row["status"] = f"failed: {str(exc)[:120]}"
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--efforts", type=int, nargs="+", default=[200, 2000])
    ap.add_argument("--gpu", default="7")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--extra", nargs="*", default=[],
                    help="extra candidate_ids to include (e.g. the BDHI lead)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    meta = {r.class_id: (r.mechanism, r.reactive_atom_smarts) for r in wh.itertuples()}
    cands = ns.crystal_positives(meta, None)
    if args.limit:
        cands = cands[:args.limit]
    if args.extra:
        by_id = {c.ident: c for c in nr.load_candidates()}
        cands += [by_id[i] for i in args.extra if i in by_id]

    rec = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    log.info("%d molecules at efforts %s", len(cands), args.efforts)

    rows = []
    for i, c in enumerate(cands, 1):
        for nrun in args.efforts:
            r = measure(c, rec, nrun, args.gpu)
            rows.append(r)
            if r["status"] == "ok":
                log.info("[%d/%d] %-22s %5d runs  all %.3f  top1 %.2f  top5 %.2f  top10 %.2f",
                         i, len(cands), c.ident[:22], nrun, r["fraction_all"],
                         r.get("top1_viable", float("nan")),
                         r.get("top5_viable", float("nan")),
                         r.get("top10_viable", float("nan")))
            else:
                log.warning("[%d/%d] %s @%d: %s", i, len(cands), c.ident, nrun, r["status"])

    df = pd.DataFrame(rows)
    dest = OUT.write("topn_vs_fraction", ".csv")
    df.to_csv(dest, index=False)

    ok = df[df.status == "ok"]
    lo, hi = min(args.efforts), max(args.efforts)
    a, b = ok[ok.nrun == lo].set_index("ident"), ok[ok.nrun == hi].set_index("ident")
    common = a.index.intersection(b.index)
    if len(common) < 3:
        print("\n  too few molecules measured at both efforts"); return

    from scipy.stats import spearmanr
    print(f"\n=== does the metric survive more searching? "
          f"({len(common)} molecules, {lo} vs {hi} runs) ===\n")
    print(f"  {'metric':<18}{'median @'+str(lo):>14}{'median @'+str(hi):>14}"
          f"{'median |change|':>17}{'rank rho':>10}")
    for name in ["fraction_all"] + [f"top{n}_viable" for n in TOP_NS]:
        if name not in a.columns or name not in b.columns:
            continue
        x, y = a.loc[common, name], b.loc[common, name]
        m = x.notna() & y.notna()
        if m.sum() < 3:
            continue
        rho = spearmanr(x[m], y[m]).statistic
        print(f"  {name:<18}{x[m].median():>14.3f}{y[m].median():>14.3f}"
              f"{(y[m]-x[m]).abs().median():>17.3f}{rho:>10.3f}")
    print("\n  A metric defined on the molecule's OWN BEST poses cannot be diluted by")
    print("  adding runs — more searching can only improve which poses occupy the top N.")
    print("  A whole-population fraction can, and D0068 measured it doing so.")


if __name__ == "__main__":
    main()
