"""
Purpose: build the 2.1.0 ranking from the v2 screen — both score repairs, both pose orderings, filtered within warhead class.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: 00_outputs/blacksmith/nac_v2/{agg,poses}_*.csv + the persisted top-20 poses
Output: 00_outputs/blacksmith/rank_v2/rank_v2_<TIER>_<N>.csv

Computes every candidate score side by side rather than picking one in advance,
because which of them predicts anything is an open question and the elevation
cohort is the instrument that will answer it.

THE FOUR SCORES

  enrichment_joint        what 2.0.0 ranked on: P(distance AND angle) over the
                          orientational null. Carried forward unchanged so the
                          comparison has a baseline.

  enrichment_conditional  P(angle | distance in window) over the same null. The
                          reactive potential is a SAMPLER that biases distance
                          (D0064), so the joint form multiplies in how well the
                          sampler worked on each molecule. Conditioning cancels
                          it: the bias applies equally to the numerator and to
                          the conditioning set. Measured spread of that bias
                          across molecules is large -- poses landing in the
                          distance window ran 38% to 97% in the first two
                          molecules alone.

  consensus_autodock      agreement among the top 10 by AutoDock energy. What
                          2.0.0 measured.

  consensus_gnina         agreement among the top 10 by gnina CNNscore. On this
                          receptor, over 82 Pin1 crystal ligands with sampling
                          held fixed, gnina's CNN puts a sub-2 A pose first
                          26.8% of the time against AutoDock's 18.3%. A better
                          ordering changes WHICH ten poses consensus sees.

THE FILTER IS A QUOTA WITHIN WARHEAD CLASS, not a library-wide bar (D0073): a
single 0.90 cut promoted rigid chemistry and demoted flexible chemistry, because
pose agreement is partly a statement about how many ways a molecule can sit.
Within-class competition removes the cross-class part of that. It does NOT
remove the within-class part -- T_3 is 100% acrylamide and consensus still
correlates with rotatable bonds at rho = -0.312 there -- so rigidity is reported
alongside every class rather than claimed to be handled.

NO VALIDATION GATE. @tt8804: the 15 crystallographic depositions are too few and
too poor to decide which chemistry to pursue. Classes are ranked on their own
terms.
"""

from __future__ import annotations

import argparse
import glob
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import nac_criterion as nac           # noqa: E402
from shared import outputs as sout                # noqa: E402
from shared import pose_consensus as pc           # noqa: E402

log = logging.getLogger("rank-v2")
OUT = sout.Topic("blacksmith", "rank_v2")
DATA = Path("/data/lab_vm/append_only/inhibition")
V2 = DATA / "00_outputs/blacksmith/nac_v2"
POSES = DATA / "00_outputs/blacksmith/nac_v2_poses"

TOP_N = 10
DEFAULT_QUOTA = 0.20
CONSENSUS_FLOOR = 0.50


def _shards(pattern: str) -> pd.DataFrame:
    fs = sorted(glob.glob(str(V2 / pattern)))
    if not fs:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)


def load_v2() -> tuple[pd.DataFrame, pd.DataFrame]:
    agg = _shards("agg_s*_*.csv")
    if agg.empty:
        raise SystemExit(f"no v2 aggregates under {V2}")
    agg = agg.drop_duplicates("ident", keep="last")
    poses = _shards("poses_s*_*.csv")
    if not poses.empty:
        poses = poses.drop_duplicates(["ident", "energy_rank"], keep="last")
    return agg, poses


def conditional_enrichment(row) -> float:
    """P(angle viable | distance in window) / P(angle viable | isotropic).

    Returns NaN rather than 0.0 when no pose reached the distance window: the
    molecule's orientation ability was not measured, which is not the same as
    measured and found absent.
    """
    if not row.n_in_range:
        return np.nan
    return (row.n_viable_given_in_range / row.n_in_range) / row.isotropic_null


def reactive_poses(ident: str, smarts: str, energies: list[float]):
    """Rebuild ReactivePose objects from the persisted SDF.

    The poses were persisted precisely so consensus could be recomputed under a
    different ordering without docking again.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    sdf = POSES / f"{ident}.sdf"
    if not sdf.is_file():
        return None
    supp = Chem.SDMolSupplier(str(sdf), removeHs=False)
    mols = [m for m in supp if m is not None]
    if len(mols) < pc.MIN_POSES_FOR_CONSENSUS or len(mols) != len(energies):
        return None
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        return None
    match = mols[0].GetSubstructMatch(patt)
    if not match:
        return None
    out = []
    for m, e in zip(mols, energies):
        conf = m.GetConformer()
        xyz = np.array([list(conf.GetAtomPosition(i)) for i in match])
        out.append(pc.ReactivePose(energy=float(e), reactive_xyz=xyz,
                                   atom_ids=tuple(match)))
    return out


def topn_viable(poses: pd.DataFrame, n: int = TOP_N) -> pd.DataFrame:
    """Fraction of the top-N poses BY ENERGY that are reaction-competent.

    This is the metric D0068 argued for and 2.0.0 never implemented. Enrichment
    counts viable poses anywhere in the 200-run population; the top N by energy
    are a different set, and measured here they often do not overlap at all --
    three of the top five T_3 molecules by conditional enrichment have NO viable
    pose among their 20 lowest-energy poses, including the top-ranked one.

    Two properties the population fraction lacks:
      - it cannot be diluted by more searching, because it is defined on the
        molecule's own best poses rather than on a per-run rate (D0068);
      - it asks the question selection actually needs answered -- is there a
        pose worth starting a trajectory from -- rather than whether a
        reaction-competent geometry exists somewhere docking did not favour.
    """
    g = poses[poses.energy_rank <= n]
    out = g.groupby("ident").agg(
        topn_scanned=("viable", "size"),
        topn_viable_n=("viable", "sum"),
        topn_best_dist=("distance", "min"),
    ).reset_index()
    out["topn_viable_frac"] = out.topn_viable_n / out.topn_scanned
    return out


def consensus_both(agg: pd.DataFrame, poses: pd.DataFrame,
                   smarts: dict) -> pd.DataFrame:
    """Consensus under AutoDock's ordering and under gnina's.

    gnina's CNNscore is HIGHER-is-better, and `consensus` selects by energy
    ascending, so the CNN score is negated to put the best pose first. Getting
    that sign wrong would silently score the ten WORST poses, so it is done in
    one place with the reason attached.
    """
    rows = []
    have_cnn = "CNNscore" in poses.columns
    for ident, g in poses.groupby("ident", sort=False):
        cls = g.warhead_class.iloc[0]
        s = smarts.get(cls)
        if s is None:
            continue
        g = g.sort_values("energy_rank")
        rec = {"ident": ident}
        rp = reactive_poses(ident, s, g.energy.tolist())
        if rp:
            try:
                rec["consensus_autodock"] = pc.consensus(rp, top_n=TOP_N).agreement
            except pc.ConsensusError:
                pass
        if have_cnn and g.CNNscore.notna().all():
            rp2 = reactive_poses(ident, s, (-g.CNNscore).tolist())
            if rp2:
                try:
                    rec["consensus_gnina"] = pc.consensus(rp2, top_n=TOP_N).agreement
                except pc.ConsensusError:
                    pass
        if len(rec) > 1:
            rows.append(rec)
    return pd.DataFrame(rows)


def add_props(df: pd.DataFrame) -> pd.DataFrame:
    smi, extra = {}, {}
    for sub, stem in (("03_t3_reinvent", "D3"), ("04_t4_combinatorial", "D4")):
        g = list((DATA / sub).glob(f"{stem}_*.parquet"))
        if not g:
            continue
        d = pd.read_parquet(max(g, key=lambda p: int(re.search(r"_(\d+)\.parquet$", p.name).group(1))))
        smi.update(dict(zip(d.candidate_id.astype(str), d.canonical_smiles)))
        for col in ("QED",):
            if col in d.columns:
                extra.setdefault(col, {}).update(dict(zip(d.candidate_id.astype(str), d[col])))
    df["smiles"] = df.ident.map(smi)
    for c, m in extra.items():
        df[c] = df.ident.map(m)

    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors as rdd
    RDLogger.DisableLog("rdApp.*")

    def rotb(s):
        if not isinstance(s, str):
            return np.nan
        m = Chem.MolFromSmiles(s)
        return np.nan if m is None else rdd.CalcNumRotatableBonds(m)

    df["rotb"] = df.smiles.map(rotb)
    df["tier"] = df.ident.str.split("_").str[0].str.upper()
    return df


def filter_and_rank(df: pd.DataFrame, score: str, cons: str,
                    quota: float, floor: float) -> pd.DataFrame:
    out = []
    for cls, g in df.groupby("warhead_class"):
        g = g.copy()
        c = g[cons]
        k = max(1, int(round(len(g) * quota)))
        cut = c.nlargest(k).min() if c.notna().any() else np.nan
        g["consensus_cut"] = cut
        g["passes"] = (c >= cut) & (c >= floor)
        s = g[g.passes].copy()
        s["class_rank"] = s[score].rank(ascending=False, method="min")
        g = g.merge(s[["ident", "class_rank"]], on="ident", how="left")
        out.append(g)
    return pd.concat(out, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--score", default="enrichment_conditional",
                    choices=("enrichment_conditional", "enrichment_joint",
                             "topn_viable_frac"))
    ap.add_argument("--consensus", default="consensus_gnina",
                    choices=("consensus_gnina", "consensus_autodock"))
    ap.add_argument("--quota", type=float, default=DEFAULT_QUOTA)
    ap.add_argument("--floor", type=float, default=CONSENSUS_FLOOR)
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    agg, poses = load_v2()
    ok = agg[agg.status == "ok"].copy()
    log.info("v2: %d measured, %d ok, %d poses rows", len(agg), len(ok), len(poses))

    ok["enrichment_joint"] = ok.enrichment
    ok["enrichment_conditional"] = ok.apply(conditional_enrichment, axis=1)
    ok["frac_in_range"] = ok.n_in_range / ok.n_poses

    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    smarts = dict(zip(wh.class_id, wh.reactive_atom_smarts))
    cons = consensus_both(agg, poses, smarts) if not poses.empty else pd.DataFrame()
    if not cons.empty:
        ok = ok.merge(cons, on="ident", how="left")
    if not poses.empty:
        ok = ok.merge(topn_viable(poses), on="ident", how="left")
    for c in ("consensus_autodock", "consensus_gnina"):
        if c not in ok.columns:
            ok[c] = np.nan

    ok = add_props(ok)
    cons_col = args.consensus if ok[args.consensus].notna().any() else "consensus_autodock"
    if cons_col != args.consensus:
        log.warning("%s unavailable; falling back to %s", args.consensus, cons_col)

    from scipy.stats import spearmanr
    if "topn_viable_frac" in ok.columns:
        v = ok.dropna(subset=["topn_viable_frac", "enrichment_joint"])
        if len(v) > 10:
            r = spearmanr(v.topn_viable_frac, v.enrichment_joint)
            log.info("top-%d viable vs enrichment_joint: rho %+.3f (n=%d)",
                     TOP_N, r.statistic, len(v))
            log.info("molecules with ZERO viable poses in their top %d: %d of %d "
                     "(%.1f%%)", TOP_N, int((v.topn_viable_frac == 0).sum()),
                     len(v), (v.topn_viable_frac == 0).mean() * 100)
    m = ok.dropna(subset=["enrichment_joint", "enrichment_conditional"])
    if len(m) > 10:
        r = spearmanr(m.enrichment_joint, m.enrichment_conditional)
        log.info("joint vs conditional rank correlation: rho %+.3f (n=%d) — "
                 "how much the sampler bias was moving the ranking",
                 r.statistic, len(m))
    mc = ok.dropna(subset=["consensus_autodock", "consensus_gnina"])
    if len(mc) > 10:
        r = spearmanr(mc.consensus_autodock, mc.consensus_gnina)
        log.info("consensus autodock vs gnina ordering: rho %+.3f (n=%d)",
                 r.statistic, len(mc))

    for tier in ("T3", "T4"):
        g = ok[ok.tier == tier]
        if g.empty:
            continue
        ranked = filter_and_rank(g, args.score, cons_col, args.quota, args.floor)
        dest = OUT.write(f"rank_v2_{tier}", ".csv")
        ranked.sort_values(["warhead_class", "class_rank"]).to_csv(dest, index=False)
        surv = ranked[ranked.passes]
        print(f"\n{'='*76}\n{tier}: {len(ranked)} molecules, {len(surv)} survive "
              f"({args.quota:.0%} within class, floor {args.floor})")
        print(f"ranked on {args.score}, filtered on {cons_col}\n{'='*76}")
        for cls, s in surv.groupby("warhead_class"):
            t = s.nsmallest(args.top, "class_rank")
            print(f"\n  --- {cls}  ({len(s)} of {(ranked.warhead_class==cls).sum()}) ---")
            print(f"  {'rk':>3} {'ident':<20}{'cond':>7}{'joint':>7}{'cons':>6}"
                  f"{'inrng':>7}{'QED':>7}")
            for r in t.itertuples():
                print(f"  {int(r.class_rank):>3} {r.ident:<20}"
                      f"{getattr(r, args.score):>7.2f}{r.enrichment_joint:>7.2f}"
                      f"{getattr(r, cons_col):>6.2f}{r.frac_in_range:>7.2f}"
                      f"{getattr(r, 'QED', float('nan')):>7.3f}")
        print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
