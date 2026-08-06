"""
Purpose: score every real T_3/T_4 candidate through the validated near-attack gate.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: the latest T_3 + T_4 candidate frames, and the reactive 3IKD receptor
Output: 00_outputs/blacksmith/nac_rank/nac_rank_s<shard>_<N>.csv (one per chunk)

THE PRODUCTION RUN. `scripts/nac_screen.py` established that the criterion
separates crystallographic covalent binders from warhead-matched measured
inactives -- chloroacetamide AUC 0.822 (p=0.0020), Michael AUC 0.800 (p=0.0288),
pooled on enrichment AUC 0.722 (p=0.0031). This applies the same gate to the
5,769 scorable candidates the four approaches actually generated.

WHAT IS RANKED, AND ON WHAT. `enrichment` -- the fraction of independent docking
runs reaching a mechanism-appropriate near-attack conformation, divided by the
fraction an isotropically-approaching nucleophile would reach by chance. The
division is not cosmetic: raw viable fractions are NOT comparable across
mechanisms (the windows differ in solid angle), and T_3 is entirely acrylamide
while T_4 spans eight classes. Ranking the two on raw fractions would rank the
mechanisms, not the molecules.

SCOPE. T_1 and T_2 are absent by design -- they carry no warhead, so the
question "can it present its warhead" is not false for them, it is undefined.
They are reported as not-applicable elsewhere and never ranked last (D0043).

RESUMABLE AND SHARDED. Each worker owns candidates where `index %% n_shards ==
shard`, writes a fresh chunk file every `--chunk` molecules, and on restart skips
any candidate_id already present in a chunk file. Nothing is appended to and
nothing is overwritten, so the append-only guarantee holds without the run having
to survive in one piece -- a 6-hour job that dies at hour 5 must not cost 5
hours.

A DRIVER'S EXIT CODE IS NOT EVIDENCE THAT THE WORK HAPPENED. Every chunk records
per-candidate status, and `--report` counts what actually landed rather than
trusting that the workers returned 0.
"""

from __future__ import annotations

import argparse
import glob
import logging
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import nac_criterion as nac         # noqa: E402
from shared import outputs as sout              # noqa: E402
import nac_screen as ns                         # noqa: E402

log = logging.getLogger("nac-rank")

OUT = sout.Topic("blacksmith", "nac_rank")
DATA = Path("/data/lab_vm/append_only/inhibition")
KNOWN = set(nac.MECHANISMS)


def latest_frame(subdir: str) -> Path:
    fs = glob.glob(str(DATA / subdir / "D*_*.parquet"))
    if not fs:
        raise SystemExit(f"no candidate frames under {subdir}")
    return Path(max(fs, key=lambda p: int(re.search(r"_(\d+)\.parquet", p).group(1))))


def load_candidates() -> list[ns.Candidate]:
    """Every generated candidate carrying a warhead of known mechanism.

    `rejected_at` is honoured where present: a candidate the generation stage
    already threw out must not reappear in a ranking. The mechanism comes from
    the warhead library keyed on class_id, never from a per-frame column, so
    T_3 and T_4 are scored against the same definitions.
    """
    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    meta = {r.class_id: (r.mechanism, r.reactive_atom_smarts) for r in wh.itertuples()}

    out = []
    for approach, subdir in (("T_3", "03_t3_reinvent"), ("T_4", "04_t4_combinatorial")):
        f = latest_frame(subdir)
        df = pd.read_parquet(f)
        if "rejected_at" in df.columns:
            df = df[df.rejected_at.isna()]
        log.info("%s: %s -> %d candidates", approach, f.name, len(df))
        for r in df.itertuples():
            cls = getattr(r, "warhead_class", None)
            if cls not in meta:
                continue
            mech, smarts = meta[cls]
            if mech not in KNOWN:
                continue
            out.append(ns.Candidate(str(r.candidate_id), r.canonical_smiles,
                                    cls, mech, smarts, approach))
    # Deterministic order, so shard membership does not depend on frame order.
    out.sort(key=lambda c: c.ident)
    log.info("scorable: %d", len(out))
    return out


def already_done() -> set[str]:
    ids = set()
    for f in glob.glob(str(OUT.dir / "nac_rank_s*.csv")):
        try:
            ids.update(pd.read_csv(f).ident.astype(str))
        except Exception as exc:                       # noqa: BLE001
            log.warning("unreadable chunk %s: %s", Path(f).name, exc)
    return ids


def score_one(c: ns.Candidate, rec_dir: Path, nrun: int, gpu: str) -> dict:
    work = Path(tempfile.mkdtemp(prefix="nacrank_"))
    try:
        per_centre = []
        for j, lig in enumerate(ns.prepare_ligand(c, work / "lig.pdbqt")):
            dlg = ns.dock(lig, rec_dir, work / f"c{j}", nrun, gpu)
            res_j = ns.measure_dlg(dlg, c)
            energies = ns.pose_energies(dlg)
            # Both parse MODEL blocks in file order, so they align — asserted
            # rather than assumed, because a silent off-by-one would attach each
            # pose's energy to a different pose's geometry and still look
            # perfectly plausible.
            if len(energies) != len(res_j):
                raise ValueError(f"{c.ident}: {len(energies)} energies for "
                                 f"{len(res_j)} poses")
            per_centre.append((res_j, energies))
        res, energies = max(per_centre, key=lambda t: nac.viable_fraction(t[0]))
        angles = np.array([r.angle for r in res])
        vf = nac.viable_fraction(res)

        # STAGE 4's raw material: the energy of poses that PASSED the geometric
        # gate, kept separate from the energy of all poses. The docking score is
        # not trusted to rank molecules (D0041, D0046, D0061); it is only being
        # asked to compare poses of ONE molecule that already cleared geometry.
        viable_e = [e for r, e in zip(res, energies) if r.viable and not np.isnan(e)]
        all_e = [e for e in energies if not np.isnan(e)]
        return {
            "ident": c.ident, "approach": c.label, "warhead_class": c.warhead_class,
            "mechanism": c.mechanism, "n_centres": len(per_centre), "n_poses": len(res),
            "viable_fraction": vf,
            "enrichment": vf / nac.isotropic_null(c.mechanism),
            "median_angle": float(np.median(angles)),
            "median_dist": float(np.median([r.distance for r in res])),
            "best_viable_dg": min(viable_e) if viable_e else np.nan,
            "median_viable_dg": float(np.median(viable_e)) if viable_e else np.nan,
            "best_dg": min(all_e) if all_e else np.nan,
            "status": "ok", "smiles": c.smiles,
        }
    except Exception as exc:                           # noqa: BLE001
        # Recorded as its own row. A candidate that could not be measured is not
        # a candidate that scored zero, and must never rank as if it had.
        return {"ident": c.ident, "approach": c.label,
                "warhead_class": c.warhead_class, "mechanism": c.mechanism,
                "viable_fraction": np.nan, "enrichment": np.nan,
                "status": f"failed: {str(exc)[:150]}", "smiles": c.smiles}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def run_shard(shard: int, n_shards: int, nrun: int, gpu: str, chunk: int) -> None:
    rec = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    cands = [c for i, c in enumerate(load_candidates()) if i % n_shards == shard]
    done = already_done()
    todo = [c for c in cands if c.ident not in done]
    log.info("shard %d/%d: %d assigned, %d already done, %d to do",
             shard, n_shards, len(cands), len(cands) - len(todo), len(todo))

    buf = []
    for k, c in enumerate(todo, 1):
        buf.append(score_one(c, rec, nrun, gpu))
        if len(buf) >= chunk or k == len(todo):
            dest = OUT.write(f"nac_rank_s{shard}", ".csv")
            pd.DataFrame(buf).to_csv(dest, index=False)
            ok = sum(r["status"] == "ok" for r in buf)
            log.info("shard %d: %d/%d done, wrote %d rows (%d ok) -> %s",
                     shard, k, len(todo), len(buf), ok, dest.name)
            buf = []


def report() -> None:
    fs = sorted(glob.glob(str(OUT.dir / "nac_rank_s*.csv")))
    if not fs:
        raise SystemExit("nothing written yet")
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    df = df.drop_duplicates("ident", keep="first")
    ok = df[df.status == "ok"].copy()
    print(f"\n=== NAC ranking: {len(df)} scored, {len(ok)} ok, "
          f"{len(df) - len(ok)} failed ({len(fs)} chunks) ===")
    if df.status.ne("ok").any():
        print("  failures:")
        print(df[df.status != "ok"].status.str.slice(0, 60).value_counts().head().to_string())
    if ok.empty:
        return
    print(f"\n  enrichment over an isotropic approach (1.0 = no better than chance)")
    for a, g in ok.groupby("approach"):
        print(f"    {a}  n={len(g):>5}  median {g.enrichment.median():>5.2f}x  "
              f"max {g.enrichment.max():>5.2f}x  above chance: "
              f"{(g.enrichment > 1).mean():.0%}")
    print(f"\n  per warhead class")
    for c, g in ok.groupby("warhead_class"):
        print(f"    {c:<24} n={len(g):>5}  median {g.enrichment.median():>5.2f}x  "
              f"max {g.enrichment.max():>5.2f}x")
    print(f"\n  TOP 25 (crystallographic positives ran 1.6-4.3x on the same gate)")
    for r in ok.nlargest(25, "enrichment").itertuples():
        print(f"    {r.enrichment:>5.2f}x  {r.viable_fraction:>6.1%}  "
              f"{r.approach}  {r.warhead_class[:20]:<21}{r.ident[:26]}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--report", action="store_true",
                    help="summarise what has actually landed, and exit")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format=f"%(levelname)s [s{args.shard}] %(message)s")
    if args.report:
        report()
        return
    run_shard(args.shard, args.n_shards, args.nrun, args.gpu, args.chunk)


if __name__ == "__main__":
    main()
