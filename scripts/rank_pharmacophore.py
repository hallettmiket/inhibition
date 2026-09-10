"""
Purpose: rank the T_2 degree-2 pools on pharmacophore similarity, with the gate's verdict attached.
Author: @tt8804 (with Claude Code)
Date: 2026-09-09
Input: each degree-2 pool's latest frame + the frozen reference binder set
Output: a new frame version per pool carrying the metric, the rank and the verdict

WHY THIS POOL. The T_2 degree-2 samples are 60,000 molecules that have never
been ranked by anything: the ATRA pool's docking run segfaulted on all six
chunks (D0060) and the Guo pool was never docked at all. They are also the one
place where a scorer needing no receptor is not merely convenient but the ONLY
thing that can run, because the docking route is blocked on a permissions fix
that has not landed in a month.

THE SEED IS HELD OUT, PER POOL. Both degree-2 seeds are themselves in the
reference binder set, and these molecules are CReM derivatives of them.
Scoring Guo's derivatives against a model set containing Guo asks how much they
still resemble their parent -- the same circularity control B4 forbids on the
novelty axis. Measured before this was written: Guo's pool sits at median 0.291
to its own seed against 0.269 to any other binder, with the seed nearest for
58.8% of molecules; ATRA's sits at 0.053 against 0.107, nearest for 16.8%. The
leak is real for one pool and not the other, which is exactly why the rule is
applied to both rather than to the one where it looked necessary.

RANKED WITHIN POOL, NEVER ACROSS. The two pools sit at different absolute
similarity levels for reasons that are properties of their seeds, not of their
chemistry -- ATRA is a polyene carrying almost no pharmacophore features. A
merged ranking would put Guo's derivatives on top by construction. This is the
same trap `config/seeds.yaml` already records for docking scores across seeds
("compare score DISTRIBUTIONS, or best-of-N at matched N"), one metric later.

SIZE DECORRELATION IS APPLIED, AND IT NEEDED A NEW COLUMN. D0049 established
that this project ranks on a size-decorrelated residual because Vina's score is
partly a heavy-atom sort. A Tanimoto to a fixed model set has its own size
dependence -- a bigger molecule has more feature pairs and more chances to
match -- so the same treatment applies. These frames carry no `HAC`, because
D2_7 was generation-only and never annotated, so it is computed here.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, RDLogger

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

RDLogger.DisableLog("rdApp.*")

from shared import compute as cmp                       # noqa: E402
from shared import io as dio                            # noqa: E402
from shared import pharmacophore as ph                  # noqa: E402
from shared import rank_shortlist as rs                 # noqa: E402

log = logging.getLogger("rank-pharmacophore")

DATA = Path("/data/lab_vm/append_only/inhibition")
STRATUM = "pharmacophore_aid504891"

#: The degree-2 pools, and the seed each one must hold out. Keyed by the seed
#: name in `config/seeds.yaml` so the SMILES is never transcribed here -- a
#: hand-copied seed SMILES that drifts from the config is a pin, and this
#: project has written that defect five times.
POOLS = {
    "atra": "02_t2_atra_crem_degree2",
    "guo_pfizer": "02e_t2_guo_crem_degree2",
}


def rank_one(seed_name: str, experiment: str, workers: int,
             quota: int) -> pd.DataFrame:
    from scripts.pharmacophore_gate import score_parallel

    seeds = yaml.safe_load(
        (REPO / "config" / "seeds.yaml").read_text())["seeds"]
    seed_smiles = seeds[seed_name]["canonical_smiles"]

    frame_path = dio.latest(DATA / experiment, "D2", ".parquet")
    if frame_path is None:
        raise SystemExit(f"no D2 frame for {experiment}")
    df = dio.read_frame(frame_path)
    log.info("[%s] %s: %d rows", seed_name, frame_path.name, len(df))

    model = ph.build_model_set(hold_out=seed_smiles)
    scored = score_parallel(list(df.canonical_smiles), model, workers)
    out = pd.concat([df.reset_index(drop=True), scored], axis=1)

    n_unscored = int(out[ph.METRIC].isna().sum())
    if n_unscored:
        log.warning("[%s] %d of %d molecules could not be fingerprinted",
                    seed_name, n_unscored, len(out))

    # HAC for the size decorrelation. Computed rather than assumed present:
    # D2_7 was generation-only and carries no descriptors.
    out["HAC"] = [
        m.GetNumHeavyAtoms() if (m := Chem.MolFromSmiles(s)) else pd.NA
        for s in out.canonical_smiles]

    # The metric is registered in RANK_DIRECTION as higher-is-better; `rank`
    # reads the direction rather than assuming it, which is what the
    # 2026-09-09 generalisation of that registry was for.
    out = rs.rank(out, metric=ph.METRIC, group_col=None,
                  min_docked=20, decorrelate_size=True)
    # GATE BEFORE SHORTLIST, NOT AFTER. `shortlist` reads `gate_verdict` off
    # the frame to write it into `shortlist_reason`; running it first records
    # every selection as "gate=UNGATED" while `rank_validated` says otherwise
    # two columns away. Two fields disagreeing about the same fact is how a
    # reader ends up trusting whichever one they happened to read.
    out = rs.attach_gate(out, STRATUM, ph.METRIC)
    out = rs.shortlist(out, quota=quota)

    # Provenance of the score, carried onto the frame rather than only logged:
    # which model molecule the top of the list is nearest to. If one binder is
    # the nearest neighbour of the whole shortlist, the ranking is a
    # single-molecule similarity search wearing a set's name.
    top = out.nsmallest(quota, "rank") if "rank" in out else out.head(quota)
    log.info("[%s] shortlist's nearest model molecule: %s", seed_name,
             top["ph2d_nearest_active"].value_counts().to_dict())

    # `write_full_frame` OWNS THE VERSIONING, THE SCHEMA CHECK AND THE
    # MANIFEST. Writing the parquet here and building a manifest beside it
    # would be a second frame-writing path free to drift from the one every
    # other stage uses -- the defect `merge_poses_onto_frame` was split out to
    # avoid. It validates against the D^i contract BEFORE writing, so a frame
    # that would violate the schema never reaches disk.
    dest = dio.write_full_frame(
        out, approach="t2", experiment=experiment,
        stage="t2_rank_pharmacophore",
        params={"metric": ph.METRIC, "higher_is_better": True,
                "seed_name": seed_name, "held_out_seed": seed_smiles,
                "model_set": list(model.name), "n_model": len(model),
                "model_rule": {"max_amide_bonds": ph.MAX_AMIDE_BONDS,
                               "max_heavy_atoms": ph.MAX_HEAVY_ATOMS},
                "gate_stratum": STRATUM,
                "gate_verdict": str(out["gate_verdict"].iloc[0]),
                "rank_validated": bool(out["rank_validated"].iloc[0]),
                "ranking_scope": "within this pool only — NOT comparable "
                                 "across seeds (different seed baselines)",
                "quota": quota, "n_unscored": n_unscored,
                "scorer": "Gobbi 2D pharmacophore fingerprint, max Tanimoto "
                          "to the held-out model set — SIMILARITY to known "
                          "actives, not a hypothesis match"},
        inputs={"frame": frame_path})
    log.info("[%s] wrote %s", seed_name, dest.name)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pools", default="atra,guo_pfizer")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--quota", type=int, default=25)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    workers = cmp.available_workers(args.workers)

    for seed_name in [s.strip() for s in args.pools.split(",") if s.strip()]:
        if seed_name not in POOLS:
            raise SystemExit(f"unknown pool {seed_name!r}; "
                             f"known: {sorted(POOLS)}")
        out = rank_one(seed_name, POOLS[seed_name], workers, args.quota)
        print(f"\n=== {seed_name} — top {args.quota} "
              f"(gate {out['gate_verdict'].iloc[0]}, "
              f"rank_validated={bool(out['rank_validated'].iloc[0])}) ===")
        cols = [c for c in ("rank", ph.METRIC, "ph2d_nearest_active", "HAC",
                            "candidate_id") if c in out.columns]
        print(out.nsmallest(args.quota, "rank")[cols].to_string(index=False))


if __name__ == "__main__":
    main()
