#!/usr/bin/env python3
"""
Purpose: does bounding the mode's diameter get the validated mode elected 5/5?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-18
Input: --replicates 5 (each an independent 500-pose screen, distinct seed)
Output: append_only/00_outputs/blacksmith/election_<candidate>/

@tt8804: "we need to get to 5/5 election".

WHERE THE 5/5 IS SUPPOSED TO COME FROM. exp/1_mode_stability established that
sampling is not the limit: the mode 3.0.0 elected and validated with a 100 ns run
(max RMSD 0.317 nm) is RECOVERED by all five replicates. It is elected 2/5 under
`conditional_eb`, 3/5 with the prior floor. The hypothesis this tests is that the
remainder is the MIXTURE problem (D0086), not a scoring problem: the validated
mode carries 29 tight viable poses alongside 58 in a halo, so its viable fraction
is diluted by poses that a correctly-cut mode would not contain, and a smaller
purer fragment elsewhere outscores it. If that is right, cutting modes properly
raises election without touching the score -- and no scoring tweak reaches 5/5.

THIS IS THE FAITHFUL COMPARISON exp/3_linkage could not make. It clusters on the
REAL stage-1 feature -- reactive-atom position plus the direction the warhead
faces, `pose_modes.features` -- by resolving the same reactive SMARTS the screen
used, from `data/reference/warhead_classes_10.csv`. exp/3 used whole-pose RMSD as
a stand-in, so its DBSCAN column was not what production produces.

Each replicate is screened with its own `--seed`, because `docking.seed` now pins
the cloud (#77) and five identical clouds would answer nothing.
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_modes as pmod                  # noqa: E402
from shared import pose_subsplit as psub               # noqa: E402
from shared import pose_cluster as pclust              # noqa: E402
from shared import run_paths as rp                     # noqa: E402

log = logging.getLogger("election")

#: Molecules with a VALIDATED 100 ns answer, and the mode that produced it.
#: A reference is a mode some run elected AND a 100 ns trajectory then confirmed;
#: anything else is an opinion. Resolved from the rank table that scored it, so
#: the geometry here is the geometry the screen recorded rather than a constant
#: retyped into this file.
REFERENCES = {
    # 3.0.0 elected mode 0; 100 ns held at max RMSD 0.317 nm.
    "t4_716800c125a7": {"topic": "nac_v5", "mode": 2, "pose_rank": 3,
                        "tier": "held", "validated": "max 0.710 nm at 100 ns"},
    # 3.1.0 elected mode 3; 100 ns held at mean 0.099 / max 0.466 nm -- the
    # tightest run in the screen (@tt8804).
    "t4_80fbed3bdf1e": {"topic": "nac_v5", "mode": 3, "pose_rank": 4,
                        "tier": "held", "validated": "mean 0.099 / max 0.466 nm"},
    # A DISCRIMINATING CASE. Both of this molecule's elevated modes were run for
    # 100 ns and they disagree: m2 HELD (mean 0.253 / max 0.812, residence
    # 1.000) while m1 LEFT at 2 ns (residence 0.025). So the benchmark is not
    # merely "find a mode" here -- it has to prefer m2 over a sibling that a
    # trajectory showed to be wrong.
    "t4_b306425b6a73": {"topic": "nac_v5", "mode": 2, "pose_rank": 3,
                        "tier": "held",
                        "validated": "mean 0.253 / max 0.812 nm; sibling m1 left at 2 ns"},
}
DEFAULT_CANDIDATE = "t4_716800c125a7"
MATCH_A, MATCH_DEG = 2.0, 45.0


def reference_for(cand: str) -> dict:
    """Warhead position and direction of the POSE that was actually validated.

    THE REPRESENTATIVE POSE, NOT THE MODE'S MEAN. This first read `dir_x/y/z`
    off the rank table -- the mean warhead direction over the mode's poses. But
    a mode is a mixture (this record's whole finding) and its poses span ~83 deg
    internally, so its mean direction is not a stable quantity to match against:
    on t4_80fbed3bdf1e the nearest mode sat 0.99-1.62 A away in every replicate
    while its mean angle wandered 24-79 deg, and a 45 deg tolerance scored two of
    five as "not recovered" when the mode was plainly there.

    The pose elevated to 100 ns is ONE geometry, exact and unambiguous. It is
    what the trajectory validated, so it is what a replicate has to find.
    """
    if cand not in REFERENCES:
        raise SystemExit(
            f"{cand} has no validated 100 ns mode to score against. A reference "
            f"must be a mode a run elected AND a trajectory confirmed; known: "
            f"{', '.join(sorted(REFERENCES))}")
    meta = REFERENCES[cand]
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    sdf = rp.BLACKSMITH / f"{meta['topic']}_poses" / f"{cand}.sdf"
    if not sdf.is_file():
        raise SystemExit(f"no representative-pose file at {sdf}")
    smarts = _smarts_for_class(_class_of(cand, meta["topic"]))
    patt = Chem.MolFromSmarts(smarts)
    want = int(meta["pose_rank"])
    for m in Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True):
        if m is None or not m.HasProp("pose_rank"):
            continue
        if int(m.GetProp("pose_rank")) != want:
            continue
        hit = m.GetSubstructMatches(patt)
        if not hit:
            raise SystemExit(f"{cand}: pose_rank {want} does not match {smarts!r}")
        f = pmod.features(m, hit[0])[0]
        return {"candidate": cand, "centroid": f[:3], "direction": f[3:],
                "note": (f"{meta['topic']} mode {meta['mode']}, pose_rank {want}"
                         f", {meta['validated']}")}
    raise SystemExit(f"{cand}: pose_rank {want} not in {sdf.name}")


def _class_of(cand: str, topic: str) -> str:
    fs = glob.glob(str(rp.BLACKSMITH / topic / "poses_s*_*.csv"))
    if not fs:
        fs = glob.glob(str(rp.BLACKSMITH / "nac_v5" / "poses_s*_*.csv"))
    p = pd.concat([pd.read_csv(f, usecols=["ident", "warhead_class"]) for f in fs],
                  ignore_index=True)
    return str(p[p.ident == cand].warhead_class.iloc[0])


def _smarts_for_class(cls: str) -> str:
    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    row = wh[wh.class_id == cls]
    if row.empty:
        raise SystemExit(f"class {cls} not in the warhead library")
    return str(row.reactive_atom_smarts.iloc[0])


def smarts_for(cand: str, topic: str) -> str:
    """The reactive SMARTS the screen used, keyed on this molecule's class."""
    fs = glob.glob(str(rp.BLACKSMITH / topic / "poses_s*_*.csv"))
    p = pd.concat([pd.read_csv(f, usecols=["ident", "warhead_class"]) for f in fs],
                  ignore_index=True)
    cls = p[p.ident == cand].warhead_class.iloc[0]
    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    row = wh[wh.class_id == cls]
    if row.empty:
        raise RuntimeError(f"{cand}: class {cls} not in the warhead library")
    return str(row.reactive_atom_smarts.iloc[0])


def screen(cand: str, topic: str, gpu: str, nrun: int, seed: int) -> None:
    only = Path(f"/tmp/election_{topic}.txt")
    only.write_text(cand + "\n")
    r = subprocess.run(
        ["nice", "-n", "19", sys.executable, str(REPO / "scripts/nac_screen_v2.py"),
         "--only", str(only), "--topic", topic, "--nrun", str(nrun),
         "--gpu", gpu, "--seed", str(seed), "--all-poses"],
        capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-1500:])


def cloud(topic: str, cand: str, smarts: str):
    """(feature matrix, heavy coords, per-pose table) for one replicate.

    Aligned on `pose_idx`, which the SDF now carries (#76) -- without it the
    cloud and its measurements are two lists at the same offset.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    fs = glob.glob(str(rp.BLACKSMITH / topic / "poses_s*_*.csv"))
    tab = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    tab = tab[tab.ident == cand].set_index("pose_idx")
    sdf = rp.BLACKSMITH / f"{topic}_allposes" / f"{cand}.sdf"
    patt = Chem.MolFromSmarts(smarts)
    idx, feats, heavy = [], [], []
    for m in Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True):
        if m is None or not m.HasProp("pose_idx"):
            continue
        i = int(m.GetProp("pose_idx"))
        if i not in tab.index:
            continue
        hit = m.GetSubstructMatches(patt)
        if not hit:
            continue
        f = pmod.features(m, hit[0])          # (1, 6) -- one conformer per record
        idx.append(i)
        feats.append(f[0])
        pos = m.GetConformer().GetPositions()
        heavy.append(np.array([pos[a.GetIdx()] for a in m.GetAtoms()
                               if a.GetAtomicNum() > 1]))
    if not idx:
        raise RuntimeError(f"{cand}/{topic}: no pose matched {smarts!r}")
    return np.array(idx), np.array(feats), np.array(heavy), tab.loc[idx]


#: @tt8804's design: "after pose splitting we have more modes that are smaller
#: more homogenous clusters of very similar poses (within 0.1 nm rmsd), there
#: should be no cap on the number of modes. next we just rank by attack angle
#: ... and we heavily decrease how much we rank by consensus."
FINE_CUT_A = 1.0          # 0.1 nm
FINE_MIN_SUB = 3


def full_split(feat, heavy, method: str) -> np.ndarray:
    """Stage 1, then stage 2. `method` names the whole splitting recipe.

    `fine` is the redesign: the same stage 1, then a 0.1 nm stage-2 cut with NO
    cap, so the number of modes is whatever the geometry says and each one is a
    tight cluster of near-identical poses rather than a bag.
    """
    if method == "hdbscan":
        # @tt8804's design, and the ORDER matters: ONE clustering step, on pose
        # similarity alone, with attack geometry used only to rank afterwards.
        # The shipped recipe clusters on the reactive atom's position and
        # direction -- which ARE the score's distance and angle terms -- so it
        # forms groups along the axis it later grades them on.
        return pclust.cluster(heavy)
    if method == "fine":
        lab = pmod.split(feat, method="dbscan")
        if (lab >= 0).sum() == 0:
            return lab
        sub, _ = psub.subdivide(lab, heavy, max_sub=None,
                                min_sub_size=FINE_MIN_SUB, cut_a=FINE_CUT_A)
        return sub
    lab = pmod.split(feat, method=method)
    if (lab >= 0).sum() == 0:
        return lab
    sub, _ = psub.subdivide(lab, heavy)
    return sub


def prior():
    """The library prior conditional_eb is fitted with, from the shipped run."""
    f = sorted(glob.glob(str(rp.BLACKSMITH / "rank_v2/rank_v2_T4_nac_v5_conditional_eb_*.csv")))[-1]
    d = pd.read_csv(f)
    n = d.n_in_range
    p = (d.n_viable_given_in_range / n.replace(0, np.nan))[n > 0]
    mu, var = float(np.nanmean(p)), float(np.nanvar(p))
    conc = max(mu * (1 - mu) / var - 1, 1e-6)
    return mu, conc


#: The off-normal acceptance window for these mechanisms
#: (`nac_criterion.PERPENDICULAR_MAX_OFF_NORMAL`): a pose is in attack geometry
#: at <= 30 deg, so LOWER is better.
ANGLE_MAX_DEG = 30.0


def score_by_angle(lab, feat, tab, ref, consensus_w: float = 0.0):
    """Rank a mode by its ATTACK ANGLE, not by how much of the cloud it holds.

    @tt8804: "we just rank by attack angle ... and we heavily decrease how much
    we rank by consensus."

    WHY THIS FOLLOWS FROM THE FINE SPLIT. `conditional_eb` scores the FRACTION of
    a mode's in-range poses that reach attack geometry -- a meaningful number
    only while a mode is a bag of dissimilar poses. Cut modes finely enough that
    every pose in one is the same pose, and that fraction goes to 0 or 1 and
    stops discriminating; what is left to compare is how good the geometry
    actually is. The angle is a property of the pose, so it neither rewards a
    mode for being large nor punishes it for being small -- which is the same
    small-n failure the empirical-Bayes floor was patching from the other side.

    `consensus_w` keeps a dial on the old behaviour: 0 ignores mode population
    entirely, higher values re-weight by it.
    """
    ang = tab.angle.to_numpy(float)
    inr = tab.in_range.to_numpy().astype(bool)
    n_tot = max(len(lab), 1)
    rows = []
    for c in sorted({int(x) for x in lab if x >= 0}):
        m = lab == c
        sel = m & inr
        if sel.sum() == 0:
            continue
        # Median over the poses that are at a reactable DISTANCE: the angle of a
        # pose 8 A away is not an attack angle at all.
        med = float(np.median(ang[sel]))
        quality = max(0.0, (ANGLE_MAX_DEG - med) / ANGLE_MAX_DEG)
        cons = float(m.sum()) / n_tot
        f = feat[m]
        cen = f[:, :3].mean(axis=0)
        d = f[:, 3:].mean(axis=0)
        d = d / (np.linalg.norm(d) or 1.0)
        rows.append({"mode": c, "n": int(m.sum()), "n_in": int(sel.sum()),
                     "median_angle": med,
                     "eb": quality * (cons ** consensus_w if consensus_w else 1.0),
                     "dist": float(np.linalg.norm(cen - ref["centroid"])),
                     "angle": float(np.degrees(np.arccos(
                         np.clip(float(d @ ref["direction"]), -1, 1))))})
    if not rows:
        return None
    d = pd.DataFrame(rows).sort_values("eb", ascending=False).reset_index(drop=True)
    d["matches_ref"] = (d.dist <= MATCH_A) & (d.angle <= MATCH_DEG)
    return d


def score_modes(lab, feat, tab, mu, conc, floor, ref):
    """conditional_eb per mode, and each mode's distance to the reference."""
    conc = max(conc, floor)
    a0, b0 = mu * conc, (1 - mu) * conc
    inr = tab.in_range.to_numpy().astype(bool)
    via = tab.viable.to_numpy().astype(bool)
    rows = []
    for c in sorted({int(x) for x in lab if x >= 0}):
        m = lab == c
        n_in = int((m & inr).sum())
        if n_in == 0:
            continue
        k = int((m & inr & via).sum())
        f = feat[m]
        cen = f[:, :3].mean(axis=0)
        d = f[:, 3:].mean(axis=0)
        d = d / (np.linalg.norm(d) or 1.0)
        ang = np.degrees(np.arccos(np.clip(float(d @ ref["direction"]), -1, 1)))
        rows.append({"mode": c, "n": int(m.sum()), "n_in": n_in,
                     "eb": ((k + a0) / (n_in + a0 + b0)),
                     "dist": float(np.linalg.norm(cen - ref["centroid"])),
                     "angle": ang})
    if not rows:
        return None
    d = pd.DataFrame(rows).sort_values("eb", ascending=False).reset_index(drop=True)
    d["matches_ref"] = (d.dist <= MATCH_A) & (d.angle <= MATCH_DEG)
    return d


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default=DEFAULT_CANDIDATE)
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--nrun", type=int, default=500)
    ap.add_argument("--floor", type=float, default=10.0)
    ap.add_argument("--seed0", type=int, default=1000,
                    help="replicate k uses seed seed0+k; change it for a fresh set")
    ap.add_argument("--consensus-w", type=float, default=0.0,
                    help="weight on mode population in the angle arm; 0 ignores it")
    ap.add_argument("--skip-screen", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    cand = args.candidate
    ref = reference_for(cand)
    log.info("reference: %s", ref["note"])
    mu, conc = prior()
    log.info("library prior: mean %.4f, fitted %.2f poses, floor %.0f", mu, conc, args.floor)

    out_rows = []
    for i in range(1, args.replicates + 1):
        topic = f"election_{cand}_s{args.seed0}_r{i}"
        if not args.skip_screen:
            log.info("replicate %d -> %s (seed %d)", i, topic, args.seed0 + i)
            screen(cand, topic, args.gpu, args.nrun, args.seed0 + i)
        smarts = smarts_for(cand, topic)
        _, feat, heavy, tab = cloud(topic, cand, smarts)
        for method in ("dbscan", "complete", "fine", "hdbscan"):
            lab = full_split(feat, heavy, method)
            d = (score_by_angle(lab, feat, tab, ref, args.consensus_w)
                 if method in ("fine", "hdbscan")
                 else score_modes(lab, feat, tab, mu, conc, args.floor, ref))
            if d is None or not d.matches_ref.any():
                out_rows.append({"replicate": i, "method": method,
                                 "modes": 0 if d is None else len(d),
                                 "ref_rank": np.nan, "elected": False})
                continue
            rank = int(d.index[d.matches_ref].min()) + 1
            out_rows.append({"replicate": i, "method": method, "modes": len(d),
                             "ref_rank": rank, "elected": rank == 1,
                             "ref_n": int(d.loc[d.matches_ref, "n"].iloc[0]),
                             "ref_eb": float(d.loc[d.matches_ref, "eb"].iloc[0]),
                             "top_eb": float(d.eb.iloc[0]),
                             "top_n": int(d.n.iloc[0])})

    t = pd.DataFrame(out_rows)
    out = rp.BLACKSMITH / f"election_{cand}_s{args.seed0}"
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "election_1.csv", index=False)
    print(f"\n  {cand}: {args.replicates} independent {args.nrun}-pose screens")
    print(f"  reference = {ref['note']}\n")
    print(t.to_string(index=False))
    print()
    for m, g in t.groupby("method"):
        ranks = sorted(int(x) for x in g.ref_rank.dropna())
        print(f"  {m:<9} elected {int(g.elected.sum())}/{len(g)}   "
              f"ranks {ranks}   modes/replicate {g.modes.mean():.1f}")
    (out / "election_1.json").write_text(json.dumps(
        {"candidate": cand, "replicates": args.replicates,
         "floor": args.floor,
         "elected": {m: int(g.elected.sum()) for m, g in t.groupby("method")}}, indent=2))
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
