#!/usr/bin/env python3
"""
Purpose: how does the number of binding modes grow with the number of docked poses?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-19
Input: --candidate, --nrun (one deep dock), --ladder (pose counts to subsample to)
Output: append_only/00_outputs/blacksmith/mode_saturation_<candidate>/

@tt8804: "run a test to see how number of modes is related to number of poses
generated go up to 100,000 maybe? id imagine it is logorithmic".

WHY IT MATTERS. `docking.n_runs: 500` is justified in config by a SAMPLING
argument -- 500 runs give >=95% probability of drawing at least one pose within
2 A of the true one. That is a statement about finding a good pose, not about
how many DISTINCT modes the cloud resolves. If mode count is still climbing at
500, then "how many ways can this molecule sit" is an artefact of how hard we
sampled, and every per-mode score (consensus = mode_size / n_poses, and
`conditional_eb` built on it) is measured against a denominator that moves.

Logarithmic growth would mean the modes saturate: new poses land in modes that
already exist, and 500 is enough. Linear growth would mean they do not, and the
mode abstraction is resolution-limited rather than chemistry-limited.

SUBSAMPLED, NOT RE-DOCKED. AutoDock-GPU's `--nrun` is N independent GA runs, so
a cloud of N poses is N i.i.d. draws and a random K-subset is distributed exactly
as a K-run cloud. One deep dock therefore gives the whole ladder, holds the
molecule and receptor fixed, and costs one run instead of a dozen. Each rung is
repeated `--reps` times on different subsets, because a single subset at K=100 is
one draw and would read as structure.

The split here is the FULL production recipe -- stage 1 then stage 2 -- so the
count is modes as the pipeline would report them, not clusters in the abstract.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import nac_screen as ns                                # noqa: E402
import nac_rank as nr                                  # noqa: E402
from shared import pose_modes as pmod                  # noqa: E402
from shared import pose_subsplit as psub               # noqa: E402
from shared import run_paths as rp                     # noqa: E402
from shared import target_config as tc                 # noqa: E402

log = logging.getLogger("saturation")

DEFAULT_LADDER = "50,100,200,500,1000,2000,5000,10000,20000,50000,100000"


def covering_number(coords, r: float) -> int:
    """Fewest poses such that every pose is within `r` A of one of them.

    GREEDY FARTHEST-POINT, the standard 2-approximation: start anywhere,
    repeatedly add the pose furthest from everything chosen so far, stop when
    nothing is further than `r`.

    WHY THIS AND NOT A CLUSTER COUNT (@tt8804). A cluster count moves for two
    opposing reasons as the cloud densifies: under-sampled regions consolidate
    from several singletons into one group (count falls), while genuinely new
    regions appear for the first time (count rises). Those cancel to an unknown
    degree, so a flat curve would not prove saturation.

    A covering number is DENSITY-INDEPENDENT by construction. It does not care
    how many times a spot was hit, only whether the spot is covered -- which is
    exactly the de-duplication question: the smallest set of poses that
    represents everything found. It plateaus if the reachable space is bounded
    and climbs forever if it is not.

    Computed incrementally (O(n) per centre, no n x n matrix), so it is usable
    at depths where the pairwise matrix would not fit.
    """
    n = len(coords)
    if n == 0:
        return 0
    dmin = np.full(n, np.inf)
    chosen, nxt = 0, 0
    while True:
        d = np.sqrt(((coords - coords[nxt]) ** 2).sum(axis=2).mean(axis=1))
        dmin = np.minimum(dmin, d)
        chosen += 1
        far = int(np.argmax(dmin))
        if dmin[far] <= r or chosen >= n:
            return chosen
        nxt = far


def pb_valid(mols) -> "np.ndarray":
    """Boolean mask over the concatenated cloud: does each pose pass PoseBusters?

    Run on the SAME rebuilt conformers the clustering sees, not on a re-read
    file, so a pose cannot be valid in one step and absent in the next.
    """
    import tempfile as _tf
    from pathlib import Path as _P
    from rdkit import Chem
    from posebusters import PoseBusters
    from shared import run_paths as _rp
    pb = PoseBusters(config="dock")
    rec = _rp.receptor_prep()
    out = []
    td = _P(_tf.mkdtemp(prefix="pb_"))
    for j, (mol, _m) in enumerate(mols):
        f = td / f"c{j}.sdf"
        w = Chem.SDWriter(str(f))
        for cid in range(mol.GetNumConformers()):
            w.write(mol, confId=cid)
        w.close()
        df = pb.bust([f], None, rec)
        cols = [c for c in df.columns if df[c].dtype == bool]
        out.append(df[cols].all(axis=1).to_numpy())
        log.info("  PoseBusters call %d/%d: %d/%d valid",
                 j + 1, len(mols), int(out[-1].sum()), len(out[-1]))
    return np.concatenate(out)


def deep_dock(cand, nrun: int, gpu: str, seed: int | None, calls: int = 1):
    """`calls` docks at `nrun` each, returning (features, heavy coords, mols).

    SEVERAL CALLS, BECAUSE ONE CANNOT GO DEEP ENOUGH. Measured on this build,
    `--nrun 5000` corrupts its own stack AND STILL WRITES A .dlg, and 10000
    exits 0 with no output (D0088). So depth past ~2,000 has to come from
    repeated calls, each with its OWN seed -- reusing one seed returns the
    identical cloud and would show a saturation curve that is an artefact of
    asking the same question twice.

    Not the screen's `--all-poses` file: that one drops DBSCAN noise, and a
    saturation curve needs the raw cloud -- the poses that fail to join a mode
    are exactly what decides whether the count is still climbing.
    """
    # THE REACTIVE receptor, built the same way the screen builds it. The
    # prepared-receptor directory does not hold `rec.reactive_config`, and
    # AutoDock fails with "Could not open dpf file" rather than anything that
    # names the real problem.
    rec_dir = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    work = Path(tempfile.mkdtemp(prefix="satur_"))
    ligs = list(ns.prepare_ligand(cand, work / "lig.pdbqt"))
    if not ligs:
        raise RuntimeError(f"{cand.ident}: ligand preparation produced nothing")
    feats, heavies, mols = [], [], []
    for i in range(calls):
        sd = None if seed is None else int(seed) + 1000 * i
        log.info("docking %s at nrun=%d, call %d/%d (seed %s)",
                 cand.ident, nrun, i + 1, calls, sd)
        dlg = ns.dock(ligs[0], rec_dir, work / f"c{i}", nrun, gpu, seed=sd)
        mol, match = ns.rebuild_and_match(dlg, cand)
        feats.append(pmod.features(mol, match))
        heavies.append(np.array([
            np.array([mol.GetConformer(c).GetPositions()[a.GetIdx()]
                      for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
            for c in range(mol.GetNumConformers())]))
        mols.append((mol, match))
    return (np.concatenate(feats), np.concatenate(heavies), mols)


#: Stage 2 builds a FULL pairwise RMSD matrix per stage-1 mode, so its cost and
#: memory are O(n^2) in that mode's population: 10,000 poses in one mode is an
#: 800 MB matrix, 100,000 is 80 GB. The ladder is capped accordingly, and the
#: ceiling is reported rather than worked around -- "we cannot measure past here
#: with the current splitter" is the honest answer to how deep this can go.
STAGE2_MAX_MODE_POSES = 12_000


def modes_at(feat, heavy, k: int, rng, recipe: str = "shipped",
             fixed_min_samples: int | None = None,
             eps: float = pmod.DEFAULT_EPS) -> dict:
    """Split a random K-subset.

    `shipped` is stage 1 + stage 2 as production runs them -- which CAPS the
    count at `max_sub` per stage-1 mode, so the curve it produces flattens by
    construction and says nothing about saturation. `fine` is the uncapped
    0.1 nm recipe (exp/4_election), where the number of modes is whatever the
    geometry supports and the question is answerable.
    """
    idx = rng.choice(len(feat), size=k, replace=False)
    # A FIXED DENSITY THRESHOLD, OR THE SHIPPED PROPORTIONAL ONE.
    #
    # `min_population_frac = 0.05` means a mode must hold 5% of ALL poses to
    # exist -- 25 poses at nrun=500, 100 at 2000, 500 at 10000. The bar rises
    # with the sample, so deeper docking never reveals a rarer mode and the
    # count cannot exceed 20 by arithmetic. That is why mode count is flat in
    # sampling depth rather than growing. A FIXED threshold asks the ordinary
    # discovery question instead: how many distinct modes exist at all, given
    # enough looks.
    frac = (fixed_min_samples / k) if fixed_min_samples else pmod.MIN_POPULATION_FRAC
    lab = pmod.split(feat[idx], eps=eps, min_population_frac=frac)
    n_stage1 = len({int(x) for x in lab if x >= 0})
    assigned = int((lab >= 0).sum())
    if assigned == 0:
        return {"poses": k, "stage1": 0, "modes": 0, "assigned_frac": 0.0,
                "largest": 0, "capped": False}
    biggest = max(int((lab == c).sum()) for c in {int(x) for x in lab if x >= 0})
    if biggest > STAGE2_MAX_MODE_POSES:
        return {"poses": k, "stage1": n_stage1, "modes": float("nan"),
                "assigned_frac": assigned / k, "largest": biggest, "capped": True}
    if recipe == "hdbscan":
        # ONE STEP, on pose similarity alone (D0088). `min_cluster_size` is
        # ABSOLUTE, unlike the shipped rule's 5%-of-sample threshold, so the bar
        # does not rise with depth -- which is the whole reason the shipped
        # curve is flat and says nothing about saturation.
        from shared import pose_cluster as pcl
        lab2 = pcl.cluster(heavy[idx])
        real = [int((lab2 == c).sum()) for c in sorted(set(lab2) - {-1})]
        n_noise = int((lab2 == -1).sum())
        cov = {f"cover_{r}a": covering_number(heavy[idx], r)
               for r in (1.0, 1.5, 2.0)}
        return {"poses": k, "stage1": len(real),
                "modes": len(real), **cov,
                # SINGLETONS COUNTED SEPARATELY (@tt8804): a noise pose is a
                # group of one, so "unique poses" is groups + singletons.
                "modes_with_singletons": len(real) + n_noise,
                "noise_frac": n_noise / k,
                "assigned_frac": 1.0 - n_noise / k,
                "largest": max(real) if real else 0, "capped": False}
    kw = (dict(max_sub=None, min_sub_size=3, cut_a=1.0) if recipe == "fine"
          else {})
    sub, _ = psub.subdivide(lab, heavy[idx], **kw)
    sizes = [int((sub == c).sum()) for c in {int(x) for x in sub if x >= 0}]
    return {"poses": k, "stage1": n_stage1, "modes": len(sizes),
            "assigned_frac": assigned / k,
            "largest": max(sizes) if sizes else 0, "capped": False}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--nrun", type=int, default=20000,
                    help="depth of the ONE dock; the ladder cannot exceed it")
    ap.add_argument("--ladder", default=DEFAULT_LADDER)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--calls", type=int, default=1,
                    help="repeat the dock this many times with distinct seeds; "
                         "depth past ~2,000 poses needs more than one call")
    ap.add_argument("--persist", action="store_true",
                    help="write the deep cloud, so covers can be recomputed "
                         "without paying for the dock again")
    ap.add_argument("--posebusters", action="store_true",
                    help="keep only PoseBusters-valid poses before clustering")
    ap.add_argument("--recipe", choices=("shipped", "fine", "hdbscan"), default="shipped",
                    help="shipped caps modes at max_sub per stage-1 mode; "
                         "fine is the uncapped 0.1 nm split")
    ap.add_argument("--eps", type=float, default=pmod.DEFAULT_EPS,
                    help="how close two poses must be to count as the same mode")
    ap.add_argument("--min-samples", type=int, default=None,
                    help="fixed poses-per-mode threshold; default is the "
                         "shipped 5%% OF THE SAMPLE, which rises with depth")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    cands = {c.ident: c for c in nr.load_candidates()}
    if args.candidate not in cands:
        raise SystemExit(f"{args.candidate} not in the candidate table")
    feat, heavy, mols = deep_dock(cands[args.candidate], args.nrun, args.gpu,
                                  args.seed, args.calls)
    if args.persist:
        from shared import outputs as _so
        from rdkit import Chem as _C
        tp = _so.Topic("blacksmith", f"deep_cloud_{args.candidate}")
        f = tp.write("cloud", ".sdf")
        w = _C.SDWriter(str(f))
        for mol, _m in mols:
            for cid in range(mol.GetNumConformers()):
                w.write(mol, confId=cid)
        w.close()
        log.info("persisted %d poses -> %s", len(feat), f)
    if args.posebusters:
        keep = pb_valid(mols)
        log.info("PoseBusters: %d of %d poses valid (%.1f%%)",
                 int(keep.sum()), len(keep), 100 * keep.mean())
        feat, heavy = feat[keep], heavy[keep]
    n = len(feat)
    log.info("cloud: %d poses", n)

    ladder = [int(x) for x in args.ladder.split(",") if int(x) <= n]
    if n not in ladder:
        ladder.append(n)
    rng = np.random.default_rng(args.seed)
    rows = []
    for k in sorted(set(ladder)):
        for r in range(args.reps if k < n else 1):
            rows.append(dict(rep=r, recipe=args.recipe,
                             **modes_at(feat, heavy, k, rng, args.recipe,
                                            args.min_samples, args.eps)))
        log.info("  %6d poses -> %.1f modes", k,
                 np.mean([x["modes"] for x in rows if x["poses"] == k]))

    d = pd.DataFrame(rows)
    _tag = (f"ms{args.min_samples}" if args.min_samples else "frac5pct") + f"_eps{args.eps}"
    out = rp.BLACKSMITH / f"mode_saturation_{args.candidate}_{args.recipe}_{_tag}"
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "saturation_1.csv", index=False)

    g = d.groupby("poses").agg(modes=("modes", "mean"), sd=("modes", "std"),
                               stage1=("stage1", "mean"),
                               assigned=("assigned_frac", "mean"),
                               largest=("largest", "mean")).reset_index()
    print(f"\n  {args.candidate}: one dock at nrun={args.nrun}, "
          f"{n} poses, {args.reps} subsets per rung\n")
    print(f"  {'poses':>8}{'stage-1':>9}{'modes':>8}{'sd':>6}{'assigned':>10}{'largest':>9}")
    for _, r in g.iterrows():
        print(f"  {int(r.poses):>8}{r.stage1:>9.1f}{r.modes:>8.1f}"
              f"{(0 if pd.isna(r.sd) else r.sd):>6.1f}{100*r.assigned:>9.0f}%{r.largest:>9.0f}")

    # LOG vs LINEAR, decided by fit rather than by eye.
    x, y = g.poses.to_numpy(float), g.modes.to_numpy(float)
    m = y > 0
    lg = np.polyfit(np.log(x[m]), y[m], 1)
    ln = np.polyfit(x[m], y[m], 1)
    r2 = lambda p, xx: 1 - (np.sum((y[m] - np.polyval(p, xx)) ** 2)
                            / np.sum((y[m] - y[m].mean()) ** 2))
    r2_log, r2_lin = r2(lg, np.log(x[m])), r2(ln, x[m])
    print(f"\n  fit R^2   logarithmic {r2_log:.3f}   linear {r2_lin:.3f}"
          f"   -> {'LOG (saturating)' if r2_log > r2_lin else 'LINEAR (not saturating)'}")
    print(f"  modes per doubling of poses: {lg[0] * np.log(2):.1f}")
    at500 = float(g[g.poses == 500].modes.iloc[0]) if (g.poses == 500).any() else float("nan")
    print(f"  at the configured n_runs=500: {at500:.1f} modes, "
          f"{100 * float(g[g.poses == 500].assigned.iloc[0]):.0f}% of poses assigned"
          if (g.poses == 500).any() else "")
    (out / "saturation_1.json").write_text(json.dumps(
        {"candidate": args.candidate, "nrun": args.nrun, "poses": int(n),
         "r2_log": r2_log, "r2_linear": r2_lin,
         "modes_per_doubling": float(lg[0] * np.log(2))}, indent=2))
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
