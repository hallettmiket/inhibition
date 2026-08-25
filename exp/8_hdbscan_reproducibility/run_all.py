#!/usr/bin/env python3
"""
Purpose: do independent dockings of one molecule produce the SAME modes?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: the 5 independent replicate clouds under election_<candidate>_r{1..5}_allposes
Output: 00_outputs/blacksmith/hdbscan_reproducibility_<candidate>/

@tt8804: "hdbscan does generate modes that are tight and homogenous good. but we
first have a lot of noise and also is it reproducable. lets check the
reproducability first."

WHY THIS AND NOT ACCURACY. There is no experimental ground truth for the state
this pipeline models -- a TRANSIENT PRE-COVALENT complex against a NAIVE
receptor. Sulfopin's crystal is the covalent adduct in an induced-fit pocket,
which is the wrong state on both axes, so it cannot serve as an answer key
(@tt8804). What is left is internal consistency, and the sharpest form of it is:
a mode that does not survive an independent draw of the pose cloud is not a
binding mode, it is a partition of one sample.

THE MEASUREMENT. Each replicate is an independent 500-run docking with a
distinct seed. Cluster each one on its own, take each mode's medoid, and ask
whether a mode found in one replicate has a counterpart in another -- matched by
the SAME in-place heavy-atom RMSD the clustering uses, so "the same mode" means
the same thing here as it does there.

COMPARED AGAINST THE SHIPPED RULE, always. "Is HDBSCAN reproducible" has no
answer on its own; it only means something beside what we run today.

NOISE IS MEASURED, NOT ASSUMED AWAY. HDBSCAN labels ~29% of poses noise
(D0088). Two different things could be true: noise is a sparse fringe that any
draw would discard, or noise is where a rare transient pose lives. The
reproducibility of the noise SET distinguishes them, so it is reported.
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import pose_cluster as pc               # noqa: E402
from shared import pose_modes as pmod               # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("hdbscan-repro")

#: A mode in one replicate is "the same mode" as one in another when their
#: medoids sit within this heavy-atom RMSD. 2.0 A is the Astex redocking
#: convention and the same bar the rest of the project uses for "the same pose".
MATCH_A = 2.0


def load_cloud(topic: str, ident: str):
    """(heavy-atom coords (n,a,3), the mols) for one replicate's persisted cloud."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    sdf = rp.BLACKSMITH / f"{topic}_allposes" / f"{ident}.sdf"
    if not sdf.is_file():
        raise SystemExit(f"no cloud at {sdf}")
    ms = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=False)
          if m is not None]
    if not ms:
        raise SystemExit(f"{sdf} holds no readable pose")
    heavy = [i for i, a in enumerate(ms[0].GetAtoms()) if a.GetAtomicNum() > 1]
    xyz = np.array([m.GetConformer().GetPositions()[heavy] for m in ms])
    return xyz, ms


def medoids(xyz: np.ndarray, lab: np.ndarray) -> dict[int, np.ndarray]:
    """The most central pose of each mode, as coordinates.

    The medoid, not the mean: a mean of coordinates is not a pose any docking
    produced and can sit in a place no conformer can reach.
    """
    out = {}
    for k in sorted(set(lab) - {-1}):
        idx = np.flatnonzero(lab == k)
        d = pc.rmsd_matrix(xyz[idx])
        out[int(k)] = xyz[idx[int(np.argmin(d.sum(axis=1)))]]
    return out


def _rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def match(m_a: dict, m_b: dict, tol: float = MATCH_A) -> dict:
    """How many of A's modes have a counterpart in B, and vice versa."""
    if not m_a or not m_b:
        return {"a_matched": 0, "b_matched": 0, "a_n": len(m_a), "b_n": len(m_b),
                "a_frac": float("nan"), "b_frac": float("nan")}
    d = np.array([[_rmsd(va, vb) for vb in m_b.values()] for va in m_a.values()])
    a_hit = (d.min(axis=1) <= tol).sum()
    b_hit = (d.min(axis=0) <= tol).sum()
    return {"a_matched": int(a_hit), "b_matched": int(b_hit),
            "a_n": len(m_a), "b_n": len(m_b),
            "a_frac": a_hit / len(m_a), "b_frac": b_hit / len(m_b)}


def class_of(topic: str, ident: str) -> str:
    """The warhead class THIS RUN recorded, never a literal.

    The first version hardcoded `bdhi_c5`. It happened to be right for this
    molecule and would have silently matched the wrong reactive atom on any
    other -- giving a shipped-rule clustering built from a different atom's
    position and direction, which is a plausible answer to a different question.
    """
    import pandas as _pd
    f = sorted((rp.BLACKSMITH / topic).glob("agg_s*.csv"))
    if not f:
        raise SystemExit(f"no aggregate under {rp.BLACKSMITH / topic}")
    d = _pd.concat([_pd.read_csv(x) for x in f], ignore_index=True)
    col = "parent_ident" if "parent_ident" in d.columns else "ident"
    cls = sorted(set(d[d[col].astype(str).str.startswith(ident)].warhead_class))
    if len(cls) != 1:
        raise SystemExit(f"{ident}: {len(cls)} warhead classes recorded: {cls}")
    return str(cls[0])


def smarts_for(cls: str) -> str:
    """The class's reactive-atom pattern, from the library."""
    from shared import warhead_library as wl
    df = wl.load()
    row = df[df.class_id.astype(str) == cls]
    if len(row) != 1:
        raise SystemExit(f"{len(row)} library rows for warhead class {cls!r}")
    return str(row.iloc[0]["reactive_atom_smarts"])


def group(xyz: np.ndarray, mols, method: str, reactive_smarts: str) -> np.ndarray:
    """Mode label per pose, under one rule."""
    if method == "hdbscan":
        return pc.cluster(xyz)
    if method == "shipped":
        # The production rule: DBSCAN over reactive-atom position + warhead
        # direction. Needs the reactive-atom match, so it is rebuilt here from
        # the same SMARTS the screen uses.
        from rdkit import Chem
        patt = Chem.MolFromSmarts(reactive_smarts)
        hit = mols[0].GetSubstructMatches(patt)
        if not hit:
            raise SystemExit(f"shipped rule: {smarts!r} does not match this molecule")
        feat = np.array([pmod.features(m, hit[0])[0] for m in mols])
        return pmod.split(feat)
    raise SystemExit(f"unknown method {method!r}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--methods", default="hdbscan,shipped")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    rows, per_rep = [], []
    for method in a.methods.split(","):
        meds, labs, clouds = {}, {}, {}
        for r in range(1, a.replicates + 1):
            topic = f"election_{a.candidate}_r{r}"
            xyz, ms = load_cloud(topic, a.candidate)
            smarts = smarts_for(class_of(topic, a.candidate))
            lab = group(xyz, ms, method, smarts)
            meds[r], labs[r], clouds[r] = medoids(xyz, lab), lab, xyz
            n_mode = len(set(lab) - {-1})
            noise = float((lab == -1).mean())
            per_rep.append(dict(method=method, rep=r, poses=len(xyz),
                                modes=n_mode, noise_frac=round(noise, 4),
                                largest=int(pd.Series(lab[lab >= 0]).value_counts().max())
                                if n_mode else 0))
            log.info("%s r%d: %d poses, %d modes, %.1f%% noise",
                     method, r, len(xyz), n_mode, noise * 100)

        for i, j in itertools.combinations(range(1, a.replicates + 1), 2):
            m = match(meds[i], meds[j])
            rows.append(dict(method=method, a=i, b=j, **m))

        # A mode present in EVERY replicate -- the reproducible core.
        core = 0
        for k, v in meds[1].items():
            if all(any(_rmsd(v, w) <= MATCH_A for w in meds[r].values())
                   for r in range(2, a.replicates + 1)):
                core += 1
        rows.append(dict(method=method, a=0, b=0, a_n=len(meds[1]), b_n=core,
                         a_matched=core, b_matched=core,
                         a_frac=core / max(1, len(meds[1])), b_frac=float("nan")))
        log.info("%s: %d of r1's %d modes are found in ALL %d replicates",
                 method, core, len(meds[1]), a.replicates)

    t = sout.Topic("blacksmith", f"hdbscan_reproducibility_{a.candidate}")
    pd.DataFrame(rows).to_csv(t.write("pairs", ".csv"), index=False)
    pd.DataFrame(per_rep).to_csv(t.write("per_replicate", ".csv"), index=False)

    print("\n" + "=" * 74)
    print(f"  MODE REPRODUCIBILITY ACROSS INDEPENDENT DOCKINGS — {a.candidate}")
    print("=" * 74)
    pr = pd.DataFrame(per_rep)
    print("\nper replicate:")
    print(pr.to_string(index=False))
    pp = pd.DataFrame([r for r in rows if r["a"] > 0])
    print("\npairwise mode recovery (fraction of A's modes with a counterpart in B,"
          f" medoid RMSD <= {MATCH_A} A):")
    for meth, g in pp.groupby("method"):
        both = pd.concat([g.a_frac, g.b_frac])
        print(f"  {meth:9s} mean {both.mean()*100:5.1f}%   "
              f"min {both.min()*100:5.1f}%   max {both.max()*100:5.1f}%")
    print("\nmodes found in ALL replicates:")
    for r in rows:
        if r["a"] == 0:
            print(f"  {r['method']:9s} {r['a_matched']} of {r['a_n']} "
                  f"({r['a_frac']*100:.0f}%)")
    print()


if __name__ == "__main__":
    main()
