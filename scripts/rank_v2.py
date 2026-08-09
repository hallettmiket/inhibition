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
#: 2.2.0 writes one row per BINDING MODE to its own topic. Selected with
#: --topic; the two must never be pooled, since their idents differ (`x` vs
#: `x_m0`) so a dedup would not collapse them.
V3 = DATA / "00_outputs/blacksmith/nac_v3"
POSES = DATA / "00_outputs/blacksmith/nac_v2_poses"

TOP_N = 10
DEFAULT_QUOTA = 0.20
CONSENSUS_FLOOR = 0.50


def load_references() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reference molecules' aggregates and poses, from the same v2 topic.

    Kept OUT of `load_v2` deliberately. References are a yardstick, not
    competitors: @tt8804 ruled the known Pin1 binders too few and too poor to
    decide which chemistry to pursue, so they must never enter a per-class quota
    or take a rank slot from a candidate. They are scored by the SAME functions
    so the comparison is like-for-like, and reported separately.
    """
    agg = _shards("agg_ref_*.csv")
    poses = _shards("poses_ref_*.csv")
    if not agg.empty:
        agg = agg[agg.status == "ok"].drop_duplicates("ident")
    if not poses.empty:
        poses = poses.drop_duplicates(["ident", "energy_rank"], keep="last")
    return agg, poses


SRC = V2          # rebound by --topic


def derived_scores(d: pd.DataFrame) -> pd.DataFrame:
    """Add `conditional_x_consensus` and `conditional_lcb` to a scored frame.

    A FUNCTION, AND NOT TWO COPIES, because the references go down a separate
    path from the candidates. When this was inline on the candidate frame only,
    the reference table came out without the very column the ranking was ordered
    by -- so the panel that exists to answer "what does the known binder score"
    could show every ingredient of the score and not the score. The references
    are a yardstick, and a yardstick has to be marked in the same units.

    WILSON LOWER BOUND — evidence enters as uncertainty, not as a weight.

    Conditional enrichment is a PROPORTION estimated from n_in_range trials, so
    its uncertainty is a solved problem rather than something to approximate
    with a multiplier. Measured on this library: modes with <10 poses in the
    window reach 12.25 while modes with >=100 top out at 7.82 -- the extremes
    are small-sample artefacts, not molecules.

    Multiplying by consensus fixes that but overcorrects: it penalises a mode
    for being a MINORITY rather than for being thinly evidenced, and no minority
    mode then reaches the top 40 -- which would bury sulfopin's case, where the
    reactive mode is the minority one.

    The Wilson bound penalises exactly the thing that should be penalised. A
    60-pose minority mode with a clean signal survives it; an 8-pose one does
    not. No weight is chosen, and the statistic is standard.
    """
    d = d.copy()
    if "consensus" in d.columns and "enrichment_conditional" in d.columns:
        d["conditional_x_consensus"] = d.enrichment_conditional * d.consensus
    need = {"n_in_range", "n_viable_given_in_range", "isotropic_null"}
    if not need.issubset(d.columns):
        return d
    n = d.n_in_range.astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        phat = d.n_viable_given_in_range / n.replace(0, np.nan)
        z = 1.96
        den = 1.0 + z * z / n
        centre = (phat + z * z / (2 * n)) / den
        half = z * np.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / den
        d["conditional_lcb"] = (centre - half).clip(lower=0.0) / d.isotropic_null

    # EMPIRICAL-BAYES SHRINKAGE — the same protection against thin evidence,
    # WITHOUT the size bias (#42).
    #
    # A lower confidence bound is conservative in one direction: it pushes an
    # uncertain estimate toward ZERO. Uncertainty scales with 1/sqrt(n) and n is
    # the mode's population, so the bound systematically demotes minority modes.
    # Measured on T4: rho(conditional_lcb, mode size) = +0.163, and mode 0 takes
    # 68.9% of the top ten against a 63.6% base rate among scorable modes.
    #
    # Shrinkage pulls an uncertain estimate toward the POPULATION MEAN instead of
    # toward zero, so a thinly-evidenced mode lands at average rather than at the
    # bottom -- there is no systematic direction, which is the property being
    # asked for. The prior is fitted to this library by method of moments, so it
    # is not a tuning knob: it is whatever the population says.
    #
    # Measured: rho(eb, size) = -0.142 and mode 0 takes 63.3% of the top ten
    # against the 63.6% base rate -- the dominant mode ends up represented at the
    # rate it actually occurs, which is what "unbiased" means here.
    ok = n > 0
    p = (d.n_viable_given_in_range / n.replace(0, np.nan))[ok]
    if len(p.dropna()) > 50:
        mu, var = float(np.nanmean(p)), float(np.nanvar(p))
        conc = max(mu * (1 - mu) / var - 1, 1e-6) if var > 0 else 1e-6
        a0, b0 = mu * conc, (1 - mu) * conc
        log.info("empirical-Bayes prior: alpha=%.2f beta=%.2f (mean %.4f, worth "
                 "~%.0f poses)", a0, b0, mu, a0 + b0)
        with np.errstate(divide="ignore", invalid="ignore"):
            d["conditional_eb"] = (((d.n_viable_given_in_range + a0)
                                    / (n + a0 + b0)) / d.isotropic_null)
        # A mode with NO in-range pose has no estimate to shrink; the prior alone
        # would hand it the population mean, which is a score invented from
        # nothing. Those modes are unscorable, not average.
        d.loc[~ok, "conditional_eb"] = np.nan
    return d


def _shards(pattern: str) -> pd.DataFrame:
    fs = sorted(glob.glob(str(SRC / pattern)))
    if not fs:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)


def load_v2() -> tuple[pd.DataFrame, pd.DataFrame]:
    agg = _shards("agg_s*_*.csv")
    if agg.empty:
        raise SystemExit(f"no aggregates under {SRC}")
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


#: Weights for the composite molecule score. EQUAL AND UNVALIDATED, on purpose:
#: @tt8804 asked to keep them equal for now and revisit. Nothing has yet shown
#: which component predicts anything, so a tuned weight would be fitted to
#: nothing. They are named and exported so the sensitivity sweep is one edit.
WEIGHTS = {"anchor_quality": 0.5, "topn_viable_frac": 0.5}


#: Moved to `shared.nac_criterion` so `nac_screen_v2` can pick a mode's
#: representative with the SAME function the ranking scores with. Re-exported
#: here under its original name; every existing caller is unaffected.
anchor_quality = nac.anchor_quality


def composite(poses: pd.DataFrame, n: int = TOP_N) -> pd.DataFrame:
    """Per-molecule weighted score over its top-N poses by energy.

    `anchor_quality` is averaged over the window rather than maximised: one
    excellent pose among ten poor ones is a different molecule from ten decent
    ones, and the max would call them equal.
    """
    g = poses[poses.energy_rank <= n].copy()
    g["anchor_quality"] = [anchor_quality(d, a, m) for d, a, m
                           in zip(g.distance, g.angle, g.mechanism)]
    out = g.groupby("ident").agg(
        anchor_quality=("anchor_quality", "mean"),
        anchor_quality_best=("anchor_quality", "max"),
    ).reset_index()
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
    ap.add_argument("--score", default="weighted_score",
                    choices=("weighted_score", "enrichment_conditional",
                             "enrichment_joint", "topn_viable_frac",
                             "conditional_x_consensus", "conditional_lcb",
                             "conditional_eb"))
    ap.add_argument("--consensus", default="consensus_gnina",
                    choices=("consensus_gnina", "consensus_autodock"))
    ap.add_argument("--quota", type=float, default=DEFAULT_QUOTA)
    ap.add_argument("--floor", type=float, default=CONSENSUS_FLOOR)
    ap.add_argument("--topic", choices=("nac_v2", "nac_v3"), default="nac_v2",
                    help="nac_v3 = 2.2.0 per-mode rows at nrun=500")
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    global SRC
    SRC = V3 if args.topic == "nac_v3" else V2
    log.info("reading %s", SRC)

    # THE GATE'S THRESHOLDS WERE CALIBRATED FOR A DIFFERENT QUANTITY.
    #
    # 2.1.0's `consensus` was agreement among the top-10 poses BY ENERGY, where
    # 0.50 means "half the favoured poses agree". 2.2.0's `consensus` is MODE
    # POPULATION over all 500 poses, where 0.50 means "this mode is the majority
    # of everything the molecule did" -- a much stronger claim about a different
    # measurement.
    #
    # Measured on the library: single-mode molecules sit at a median 0.860 and
    # 99.9% clear 0.50 trivially; multi-mode molecules sit at 0.142 and only 33%
    # of their modes clear it, removing 142 molecules outright and every MINORITY
    # mode of the rest. Those minority modes are exactly what pose splitting
    # exists to surface -- sulfopin's reactive mode is the minority one.
    #
    # Coherence is now guaranteed by CONSTRUCTION: clustering labels anything
    # below MIN_POPULATION_FRAC (5%) as noise, so a mode that exists is already a
    # real mode. Re-asking the question with a threshold built for another
    # measurement double-counts it and biases against multi-modal binders.
    if args.topic == "nac_v3":
        if args.consensus == "consensus_gnina":
            args.consensus = "consensus"          # the mode population itself
        if args.floor == CONSENSUS_FLOOR:
            args.floor = 0.05                     # "the mode is real"
        if args.quota == DEFAULT_QUOTA:
            args.quota = 1.0                      # rank every real mode
        log.info("nac_v3: gating on `%s` at floor %.2f, quota %.2f",
                 args.consensus, args.floor, args.quota)

    agg, poses = load_v2()
    ok = agg[agg.status == "ok"].copy()
    log.info("v2: %d measured, %d ok, %d poses rows", len(agg), len(ok), len(poses))

    ok["enrichment_joint"] = ok.enrichment
    ok["enrichment_conditional"] = ok.apply(conditional_enrichment, axis=1)
    # PER-MODE ROWS DIVIDE BY THE MODE, NOT THE MOLECULE.
    #
    # On a 2.2.0 row `n_in_range` counts poses that are in range AND in this
    # mode, while `n_poses` is everything the molecule produced -- so the ratio
    # would be neither fraction, and would shrink simply because a molecule had
    # more modes. 2.1.0 rows have no `n_poses_mode` and fall back to `n_poses`,
    # which is correct for them.
    denom = (ok.n_poses_mode if "n_poses_mode" in ok.columns else ok.n_poses)
    ok["frac_in_range"] = ok.n_in_range / denom.replace(0, np.nan)

    # THE SCORE @tt8804 ACTUALLY SPECIFIED: "rank the poses by how well anchored
    # they are WEIGHTED BY THE CONSENSUS SCORE of that pose."
    #
    # It is also the fix for a small-sample artefact. Conditional enrichment is a
    # ratio, so a mode with few poses in the distance window takes extreme
    # values: measured on this library, modes with <10 poses in window reach
    # 12.25 while modes with >=100 top out at 7.82. Ungated, a mode holding 6% of
    # a molecule's poses was ranking first in its class on 8 poses' worth of
    # evidence.
    #
    # Multiplying by consensus makes the weight of evidence part of the score
    # rather than a threshold bolted beside it: a 6%-population mode scoring 12
    # lands at 0.74, an 86%-population mode scoring 2 lands at 1.76. No cut-off
    # is invented, and a minority mode can still win if its enrichment is high
    # ENOUGH -- which is what sulfopin's case requires.
    ok = derived_scores(ok)

    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    smarts = dict(zip(wh.class_id, wh.reactive_atom_smarts))
    # 2.2.0 ALREADY HAS ITS CONSENSUS, AND IT IS A DIFFERENT QUANTITY.
    #
    # `consensus_both` recomputes agreement among the top-N poses BY ENERGY --
    # the ordering #23/#30 showed is uninformative and that this version exists
    # to remove. Under nac_v3 the screen has already written `consensus` as the
    # MODE POPULATION over all 500 poses, with no energy in it.
    #
    # It also cannot merge: the poses frame is keyed on the parent molecule while
    # these rows are keyed on <parent>_m<k>, so `on="ident"` matches nothing and
    # every consensus column comes back NaN. Measured, not guessed -- with the
    # recompute in place, 0 of 6,302 molecules passed the gate, because the gate
    # compares NaN and NaN >= NaN is False.
    per_mode = args.topic == "nac_v3"
    cons = (pd.DataFrame() if per_mode or poses.empty
            else consensus_both(agg, poses, smarts))
    if not cons.empty:
        ok = ok.merge(cons, on="ident", how="left")
    if not poses.empty and not per_mode:
        # Same keying problem, and the same objection: both read a top-N-by-energy
        # window. `topn_viable_frac` is the metric that scored Sulfopin 0.000.
        ok = ok.merge(topn_viable(poses), on="ident", how="left")
        ok = ok.merge(composite(poses), on="ident", how="left")
        # the weighted score @tt8804 asked for: enrichment's question, scored
        # continuously via anchoring, combined with how much of the favoured
        # pose window is reaction-competent at all.
        ok["weighted_score"] = sum(
            WEIGHTS[c] * ok[c].fillna(0.0) for c in WEIGHTS if c in ok.columns)
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

    # references, through the identical scoring functions
    ragg, rposes = load_references()
    if not ragg.empty:
        ragg["enrichment_joint"] = ragg.enrichment
        ragg["enrichment_conditional"] = ragg.apply(conditional_enrichment, axis=1)
        ragg["frac_in_range"] = ragg.n_in_range / ragg.n_poses
        if not rposes.empty:
            ragg = ragg.merge(topn_viable(rposes), on="ident", how="left")
            ragg = ragg.merge(composite(rposes), on="ident", how="left")
            rc = consensus_both(ragg, rposes, smarts)
            if not rc.empty:
                ragg = ragg.merge(rc, on="ident", how="left")
        for c in WEIGHTS:
            if c not in ragg.columns:
                ragg[c] = np.nan
        ragg["weighted_score"] = sum(
            WEIGHTS[c] * ragg[c].fillna(0.0) for c in WEIGHTS)
        # The SAME derived scores the candidates get, or the references cannot be
        # read on the axis the candidates were ranked on.
        ragg = derived_scores(ragg)
        ragg["reference_name"] = ragg.ident.str.replace("^ref_", "", regex=True) \
                                          .str.split("__").str[0]
        # Carry the SMILES and the provenance across. Without smiles the GUI can
        # draw a candidate's structure and not a reference's, which makes the one
        # comparison the panel exists for the one thing it cannot show.
        ref_csv = pd.read_csv(REPO / "data/reference/pin1_reference_binders_4.csv")
        meta = ref_csv.set_index("name")
        for col in ("canonical_smiles", "potency", "tier", "citation", "pdb"):
            if col in meta.columns:
                ragg[col if col != "canonical_smiles" else "smiles"] = \
                    ragg.reference_name.map(meta[col])
        dest = OUT.write(f"rank_v2_REF_{args.score}", ".csv")
        ragg.to_csv(dest, index=False)
        log.info("references: %d scored -> %s", len(ragg), Path(dest).name)

    for tier in ("T3", "T4"):
        g = ok[ok.tier == tier]
        if g.empty:
            continue
        ranked = filter_and_rank(g, args.score, cons_col, args.quota, args.floor)
        # the SCORE goes in the filename. Writing every ordering to
        # rank_v2_<tier>_<n> and letting consumers take the newest is how the
        # overnight selection silently used enrichment_joint -- the 2.0.0
        # quantity -- when the chain had ranked on topn_viable_frac first.
        dest = OUT.write(f"rank_v2_{tier}_{args.score}", ".csv")
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
